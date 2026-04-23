mod parallel;
mod spool;

use crate::align::competitive::resolve_bundle_inputs;
use crate::adjudicate::FragmentAdjudicator;
use crate::aggregate::CandidateAggregator;
use crate::align::{CompetitiveAligner, CompetitivePreset};
use crate::calibration::{
    load_profile, load_reference_bundle, LoadedProfile, LoadedReferenceBundle,
    ReleaseGate as CalibrationReleaseGate,
};
use crate::cli::{ScreenArgs, ScreenMode};
use crate::decision::{
    apply_negative_control, load_benchmark_gates, DecisionContext, DecisionEngine, DecisionOutcome,
    NegativeControlDecisionInput, ReleaseContext, ReleaseGate as DecisionReleaseGate,
};
use crate::error::{Result, RvScreenError};
use crate::io::{
    FragmentReaderFactory, FragmentRecord as IoFragmentRecord, ScreenInput,
};
use crate::pipeline::parallel::RoundParallelism;
use crate::pipeline::spool::RepresentativeRoundSpool;
use crate::qc::{FragmentRecord as QcFragmentRecord, QcFilter, QcStats};
use crate::report::write_report_bundle;
use crate::sampling::{RepresentativeSampler, Sampler, StreamingSampler};
use crate::types::{CandidateCall, RoundRecord, RunManifest, SampleSummary, StopReason};
use std::path::PathBuf;

const SUPPORTED_BACKEND: &str = "minimap2";
const SUPPORTED_PRESET: &str = "sr-conservative";
const MAX_ACTIVE_SAMPLING_ROUNDS: usize = 64;

#[derive(Debug, Clone, PartialEq)]
pub struct ScreenRunOutcome {
    pub summary: SampleSummary,
    pub candidates: Vec<CandidateCall>,
    pub rounds: Vec<RoundRecord>,
    pub manifest: RunManifest,
}

pub struct PreparedScreenRunner {
    loaded_bundle: LoadedReferenceBundle,
    loaded_profile: LoadedProfile,
    benchmark_gates: Option<CalibrationReleaseGate>,
    input: ScreenInput,
    round_targets: Vec<u64>,
    qc_filter: QcFilter,
    adjudicator: FragmentAdjudicator,
    decision_engine: DecisionEngine,
    aligner: CompetitiveAligner,
    thread_budget: usize,
}

pub fn run_screen(args: &ScreenArgs) -> Result<ScreenRunOutcome> {
    PreparedScreenRunner::prepare(args)?.run(args)
}

impl PreparedScreenRunner {
    pub fn prepare(args: &ScreenArgs) -> Result<Self> {
        let loaded_bundle = load_reference_bundle(&args.reference_bundle)?;
        let loaded_profile = load_profile(&args.calibration_profile)?;
        let benchmark_gates = load_benchmark_gates(&loaded_profile.profile_dir)?;
        validate_profile_compatibility(&loaded_bundle, &loaded_profile.profile)?;

        let input = ScreenInput::from_cli_inputs(&args.input, &loaded_bundle.bundle_dir)?;
        validate_supported_input(&loaded_profile.profile.supported_input, input.kind_label())?;

        let round_targets = active_round_targets(
            &loaded_profile.profile.sampling.rounds,
            loaded_profile.profile.sampling.max_rounds,
        )?;
        let qc_filter = QcFilter::default();
        let adjudicator = FragmentAdjudicator::new(&loaded_profile.profile.fragment_rules);
        let decision_engine = DecisionEngine::new(
            &loaded_profile.profile.decision_rules,
            &loaded_profile.profile.candidate_rules,
        );
        let (reference_input, manifest_path) = resolve_bundle_inputs(&loaded_bundle.bundle_dir)?;
        let aligner = CompetitiveAligner::from_resolved_inputs(
            reference_input,
            manifest_path,
            CompetitivePreset::SrConservative,
            args.threads,
        )?;
        let thread_budget = aligner.thread_budget();

        Ok(Self {
            loaded_bundle,
            loaded_profile,
            benchmark_gates,
            input,
            round_targets,
            qc_filter,
            adjudicator,
            decision_engine,
            aligner,
            thread_budget,
        })
    }

    pub fn run(&self, args: &ScreenArgs) -> Result<ScreenRunOutcome> {
        self.execute(args, true)
    }

    pub fn run_without_report(&self, args: &ScreenArgs) -> Result<ScreenRunOutcome> {
        self.execute(args, false)
    }

    pub fn aligner_reference_input(&self) -> &std::path::Path {
        self.aligner.reference_input()
    }

    pub fn aligner_uses_prebuilt_index(&self) -> bool {
        self.aligner.uses_prebuilt_index()
    }

    fn execute(&self, args: &ScreenArgs, write_report: bool) -> Result<ScreenRunOutcome> {
        let negative_control =
            NegativeControlDecisionInput::from_optional_path(args.negative_control.as_deref())?;
        let parallelism = RoundParallelism::new(args.threads)?;
        if parallelism.threads() != self.thread_budget {
            return Err(RvScreenError::validation(
                "threads",
                format!(
                    "prepared screen runner was initialized for {} worker threads, but execution requested {}",
                    self.thread_budget,
                    parallelism.threads()
                ),
            ));
        }
        let round_dependencies = RoundDependencies {
            qc_filter: &self.qc_filter,
            aligner: &self.aligner,
            adjudicator: &self.adjudicator,
            decision_engine: &self.decision_engine,
            negative_control: &negative_control,
            parallelism: &parallelism,
        };

        let execution = match args.mode {
            ScreenMode::Representative => execute_representative_rounds(
                &self.input,
                &self.loaded_profile.profile.seed,
                &self.round_targets,
                &round_dependencies,
            )?,
            ScreenMode::Streaming => {
                execute_streaming_rounds(&self.input, &self.round_targets, &round_dependencies)?
            }
        };

        let release_status = DecisionReleaseGate::evaluate(&ReleaseContext {
            sampling_mode: args.mode,
            negative_control: &negative_control,
            calibration_profile: &self.loaded_profile.profile,
            benchmark_gates: self.benchmark_gates.as_ref(),
        });
        let summary = SampleSummary {
            sample_id: infer_sample_id(&args.input),
            reference_bundle_version: self.loaded_bundle.bundle.version.clone(),
            calibration_profile_version: self.loaded_profile.profile.profile_id.clone(),
            backend: self.loaded_profile.profile.backend.clone(),
            seed: self.loaded_profile.profile.seed,
            input_fragments: execution.qc_stats.total_fragments,
            qc_passing_fragments: execution.qc_stats.passed_fragments,
            sampled_fragments: execution.final_round_sampled_fragments,
            rounds_run: execution.rounds.len() as u64,
            stop_reason: execution.stop_reason,
            decision_status: execution.final_outcome.status,
            release_status,
        };
        let manifest = RunManifest {
            reference_bundle_version: self.loaded_bundle.bundle.version.clone(),
            calibration_profile_version: self.loaded_profile.profile.profile_id.clone(),
            backend: self.loaded_profile.profile.backend.clone(),
            seed: self.loaded_profile.profile.seed,
            sampling_mode: match args.mode {
                ScreenMode::Representative => "representative".to_string(),
                ScreenMode::Streaming => "streaming".to_string(),
            },
            negative_control: negative_control
                .manifest(self.loaded_profile.profile.negative_control_required),
            input_files: args
                .input
                .iter()
                .map(|path| path.display().to_string())
                .collect(),
        };

        if write_report {
            write_report_bundle(
                &args.out,
                &summary,
                &execution.candidates,
                &execution.rounds,
                &manifest,
            )?;
        }

        Ok(ScreenRunOutcome {
            summary,
            candidates: execution.candidates,
            rounds: execution.rounds,
            manifest,
        })
    }
}

#[derive(Debug, Clone)]
struct PipelineExecution {
    qc_stats: QcStats,
    final_round_sampled_fragments: u64,
    final_outcome: DecisionOutcome,
    stop_reason: StopReason,
    candidates: Vec<CandidateCall>,
    rounds: Vec<RoundRecord>,
}

#[derive(Debug, Clone)]
struct RoundExecution {
    qc_stats: QcStats,
    sampled_fragments: u64,
    candidates: Vec<CandidateCall>,
}

struct RepresentativePassExecution {
    qc_stats: QcStats,
    round_spool: RepresentativeRoundSpool,
}

struct RoundDependencies<'a> {
    qc_filter: &'a QcFilter,
    aligner: &'a CompetitiveAligner,
    adjudicator: &'a FragmentAdjudicator,
    decision_engine: &'a DecisionEngine,
    negative_control: &'a NegativeControlDecisionInput,
    parallelism: &'a RoundParallelism,
}

fn execute_representative_rounds(
    input: &dyn FragmentReaderFactory,
    seed: &u64,
    round_targets: &[u64],
    dependencies: &RoundDependencies<'_>,
) -> Result<PipelineExecution> {
    let sampler = RepresentativeSampler::new(*seed, round_targets.to_vec())?;

    if sampler.round_count() == 1 {
        let execution = execute_single_round_representative_pass(input, &sampler, dependencies)?;
        let background_adjusted =
            apply_negative_control(&execution.candidates, dependencies.negative_control);
        let evaluated = dependencies
            .decision_engine
            .evaluate_candidates(&background_adjusted);
        let outcome = dependencies.decision_engine.decide(
            &evaluated,
            &DecisionContext::new(1, 1),
        );
        let rounds = vec![round_record(&evaluated, execution.sampled_fragments, &outcome)];

        return Ok(PipelineExecution {
            qc_stats: execution.qc_stats,
            final_round_sampled_fragments: execution.sampled_fragments,
            stop_reason: outcome
                .stop_reason
                .clone()
                .unwrap_or(StopReason::MaxRoundsReached),
            final_outcome: outcome,
            candidates: evaluated,
            rounds,
        });
    }

    let RepresentativePassExecution {
        qc_stats,
        round_spool,
    } = execute_representative_pass(input, &sampler, dependencies)?;
    let mut cumulative_aggregator = CandidateAggregator::new();
    let mut cumulative_sampled_fragments = 0u64;
    let mut final_candidates = Vec::new();
    let mut rounds = Vec::new();
    let mut final_outcome = None;

    for (round_index, contribution) in round_spool.into_round_contributions().into_iter().enumerate() {
        cumulative_sampled_fragments = cumulative_sampled_fragments
            .saturating_add(contribution.sampled_fragments);
        cumulative_aggregator.merge_from(contribution.aggregator)?;
        cumulative_aggregator.set_total_sampled_fragments(cumulative_sampled_fragments);

        let round_candidates = cumulative_aggregator.finalize(&qc_stats);
        let background_adjusted =
            apply_negative_control(&round_candidates, dependencies.negative_control);
        let evaluated = dependencies
            .decision_engine
            .evaluate_candidates(&background_adjusted);
        let outcome = dependencies.decision_engine.decide(
            &evaluated,
            &DecisionContext::new((round_index + 1) as u64, round_targets.len() as u64),
        );
        rounds.push(round_record(
            &evaluated,
            cumulative_sampled_fragments,
            &outcome,
        ));
        final_candidates = evaluated;
        final_outcome = Some(outcome.clone());

        if outcome.should_stop() {
            return Ok(PipelineExecution {
                qc_stats,
                final_round_sampled_fragments: cumulative_sampled_fragments,
                stop_reason: outcome
                    .stop_reason
                    .clone()
                    .unwrap_or(StopReason::MaxRoundsReached),
                final_outcome: outcome,
                candidates: final_candidates,
                rounds,
            });
        }
    }

    let final_outcome = final_outcome.ok_or_else(|| {
        RvScreenError::validation("sampling.rounds", "no representative rounds were executed")
    })?;
    Ok(PipelineExecution {
        qc_stats,
        final_round_sampled_fragments: rounds
            .last()
            .map(|round| round.sampled_fragments)
            .unwrap_or(0),
        stop_reason: final_outcome
            .stop_reason
            .clone()
            .unwrap_or(StopReason::MaxRoundsReached),
        final_outcome,
        candidates: final_candidates,
        rounds,
    })
}

fn execute_streaming_rounds(
    input: &dyn FragmentReaderFactory,
    round_targets: &[u64],
    dependencies: &RoundDependencies<'_>,
) -> Result<PipelineExecution> {
    let mut qc_baseline = None;
    let mut final_candidates = Vec::new();
    let mut rounds = Vec::new();
    let mut final_outcome = None;

    for (round_index, target) in round_targets.iter().copied().enumerate() {
        let max_fragments = usize::try_from(target).map_err(|_| {
            RvScreenError::validation(
                "sampling.rounds",
                format!("streaming round target {target} exceeds usize capacity on this platform"),
            )
        })?;
        let sampler = StreamingSampler::new(max_fragments);
        let execution = execute_round(input, &sampler, dependencies)?;
        let qc_stats = merge_or_validate_qc_baseline(&mut qc_baseline, execution.qc_stats)?;
        let background_adjusted =
            apply_negative_control(&execution.candidates, dependencies.negative_control);
        let evaluated = dependencies
            .decision_engine
            .evaluate_candidates(&background_adjusted);
        let outcome = dependencies.decision_engine.decide(
            &evaluated,
            &DecisionContext::new((round_index + 1) as u64, round_targets.len() as u64),
        );
        rounds.push(round_record(
            &evaluated,
            execution.sampled_fragments,
            &outcome,
        ));
        final_candidates = evaluated;
        final_outcome = Some(outcome.clone());

        if outcome.should_stop() {
            return Ok(PipelineExecution {
                qc_stats,
                final_round_sampled_fragments: execution.sampled_fragments,
                stop_reason: outcome
                    .stop_reason
                    .clone()
                    .unwrap_or(StopReason::MaxRoundsReached),
                final_outcome: outcome,
                candidates: final_candidates,
                rounds,
            });
        }
    }

    let final_outcome = final_outcome.ok_or_else(|| {
        RvScreenError::validation("sampling.rounds", "no streaming rounds were executed")
    })?;
    Ok(PipelineExecution {
        qc_stats: qc_baseline.unwrap_or_default(),
        final_round_sampled_fragments: rounds
            .last()
            .map(|round| round.sampled_fragments)
            .unwrap_or(0),
        stop_reason: final_outcome
            .stop_reason
            .clone()
            .unwrap_or(StopReason::MaxRoundsReached),
        final_outcome,
        candidates: final_candidates,
        rounds,
    })
}

fn execute_round(
    input: &dyn FragmentReaderFactory,
    sampler: &StreamingSampler,
    dependencies: &RoundDependencies<'_>,
) -> Result<RoundExecution> {
    let reader = input.open_reader()?;
    let mut qc_stats = QcStats::default();
    let mut qc_passing_index = 0usize;
    let mut sampled_fragments = 0u64;
    let mut submission_sequence = 0usize;

    let mut aggregator = std::thread::scope(|scope| -> Result<CandidateAggregator> {
        let executor = dependencies.parallelism.start_round_executor(
            scope,
            dependencies.aligner,
            dependencies.adjudicator,
        );

        for fragment in reader {
            let fragment = fragment?;
            let qc_fragment = bridge_qc_fragment(&fragment);
            let qc_result = dependencies
                .qc_filter
                .filter_and_record(&qc_fragment, &mut qc_stats);
            if !qc_result.is_pass() {
                continue;
            }

            let should_accept = sampler.should_accept(&fragment.fragment_key, qc_passing_index);
            qc_passing_index = qc_passing_index.saturating_add(1);
            if !should_accept {
                continue;
            }

            sampled_fragments = sampled_fragments.saturating_add(1);
            executor.submit(submission_sequence, fragment)?;
            submission_sequence = submission_sequence.saturating_add(1);

            if sampled_fragments >= sampler.max_fragments() as u64 {
                break;
            }
        }

        executor.finish()
    })?;

    aggregator.set_total_sampled_fragments(sampled_fragments);
    let candidates = aggregator.finalize(&qc_stats);

    Ok(RoundExecution {
        qc_stats,
        sampled_fragments,
        candidates,
    })
}

fn execute_representative_pass(
    input: &dyn FragmentReaderFactory,
    sampler: &RepresentativeSampler,
    dependencies: &RoundDependencies<'_>,
) -> Result<RepresentativePassExecution> {
    let reader = input.open_reader()?;
    let mut qc_stats = QcStats::default();
    let mut submission_sequence = 0usize;

    let round_spool = std::thread::scope(|scope| -> Result<RepresentativeRoundSpool> {
        let executor = dependencies.parallelism.start_representative_round_executor(
            scope,
            dependencies.aligner,
            dependencies.adjudicator,
            RepresentativeRoundSpool::new(sampler.round_count())?,
        );

        for fragment in reader {
            let fragment = fragment?;
            let qc_fragment = bridge_qc_fragment(&fragment);
            let qc_result = dependencies
                .qc_filter
                .filter_and_record(&qc_fragment, &mut qc_stats);
            if !qc_result.is_pass() {
                continue;
            }

            let bucket = sampler.bucket_for(&fragment.fragment_key);
            let Some(entry_round) = sampler.entry_round_for_bucket(bucket) else {
                continue;
            };

            executor.submit(submission_sequence, entry_round, fragment)?;
            submission_sequence = submission_sequence.saturating_add(1);
        }

        executor.finish()
    })?;

    Ok(RepresentativePassExecution {
        qc_stats,
        round_spool,
    })
}

fn execute_single_round_representative_pass(
    input: &dyn FragmentReaderFactory,
    sampler: &RepresentativeSampler,
    dependencies: &RoundDependencies<'_>,
) -> Result<RoundExecution> {
    let reader = input.open_reader()?;
    let mut qc_stats = QcStats::default();
    let mut sampled_fragments = 0u64;
    let mut submission_sequence = 0usize;

    let mut aggregator = std::thread::scope(|scope| -> Result<CandidateAggregator> {
        let executor = dependencies.parallelism.start_round_executor(
            scope,
            dependencies.aligner,
            dependencies.adjudicator,
        );

        for fragment in reader {
            let fragment = fragment?;
            let qc_fragment = bridge_qc_fragment(&fragment);
            let qc_result = dependencies
                .qc_filter
                .filter_and_record(&qc_fragment, &mut qc_stats);
            if !qc_result.is_pass() {
                continue;
            }

            if !sampler.should_accept(&fragment.fragment_key, 0) {
                continue;
            }

            sampled_fragments = sampled_fragments.saturating_add(1);
            executor.submit(submission_sequence, fragment)?;
            submission_sequence = submission_sequence.saturating_add(1);
        }

        executor.finish()
    })?;

    aggregator.set_total_sampled_fragments(sampled_fragments);
    let candidates = aggregator.finalize(&qc_stats);

    Ok(RoundExecution {
        qc_stats,
        sampled_fragments,
        candidates,
    })
}

fn validate_profile_compatibility(
    bundle: &LoadedReferenceBundle,
    profile: &crate::types::ProfileToml,
) -> Result<()> {
    if profile.reference_bundle != bundle.bundle.version {
        return Err(RvScreenError::validation(
            "calibration_profile.reference_bundle",
            format!(
                "calibration profile expects reference bundle `{}`, but loaded bundle is `{}`",
                profile.reference_bundle, bundle.bundle.version
            ),
        ));
    }

    if profile.backend != SUPPORTED_BACKEND {
        return Err(RvScreenError::validation(
            "calibration_profile.backend",
            format!(
                "unsupported backend `{}`; Task 18 only wires `{SUPPORTED_BACKEND}`",
                profile.backend
            ),
        ));
    }

    if profile.preset != SUPPORTED_PRESET {
        return Err(RvScreenError::validation(
            "calibration_profile.preset",
            format!(
                "unsupported preset `{}`; Task 18 only wires `{SUPPORTED_PRESET}`",
                profile.preset
            ),
        ));
    }

    Ok(())
}

fn validate_supported_input(supported_inputs: &[String], requested: &str) -> Result<()> {
    if supported_inputs.is_empty() || supported_inputs.iter().any(|value| value == requested) {
        return Ok(());
    }

    Err(RvScreenError::validation(
        "calibration_profile.supported_input",
        format!(
            "input kind `{requested}` is not listed in calibration profile supported_input [{}]",
            supported_inputs.join(", ")
        ),
    ))
}

fn active_round_targets(rounds: &[u64], max_rounds: u64) -> Result<Vec<u64>> {
    if rounds.is_empty() {
        return Err(RvScreenError::validation(
            "sampling.rounds",
            "screen requires at least one sampling round",
        ));
    }

    if max_rounds == 0 {
        return Err(RvScreenError::validation(
            "sampling.max_rounds",
            "sampling.max_rounds must be greater than zero",
        ));
    }

    let active_len = usize::try_from(max_rounds)
        .unwrap_or(usize::MAX)
        .min(rounds.len());
    let active = rounds[..active_len].to_vec();

    if active.is_empty() {
        return Err(RvScreenError::validation(
            "sampling.rounds",
            "sampling.max_rounds truncated all configured rounds",
        ));
    }

    if active.len() > MAX_ACTIVE_SAMPLING_ROUNDS {
        return Err(RvScreenError::validation(
            "sampling.max_rounds",
            format!(
                "screen supports at most {MAX_ACTIVE_SAMPLING_ROUNDS} active sampling rounds per run; got {}",
                active.len()
            ),
        ));
    }

    Ok(active)
}

fn bridge_qc_fragment(fragment: &IoFragmentRecord) -> QcFragmentRecord<'_> {
    QcFragmentRecord::paired(&fragment.r1_seq, &fragment.r2_seq)
}

fn merge_or_validate_qc_baseline(
    baseline: &mut Option<QcStats>,
    observed: QcStats,
) -> Result<QcStats> {
    match baseline {
        Some(existing) if *existing != observed => Err(RvScreenError::validation(
            "input",
            format!(
                "re-opened round observed inconsistent input/QC totals: first pass {:?}, later pass {:?}",
                existing, observed
            ),
        )),
        Some(existing) => Ok(*existing),
        slot @ None => {
            *slot = Some(observed);
            Ok(observed)
        }
    }
}

fn round_record(
    candidates: &[CandidateCall],
    sampled_fragments: u64,
    outcome: &DecisionOutcome,
) -> RoundRecord {
    RoundRecord {
        sampled_fragments,
        accepted_virus: candidates
            .iter()
            .map(|candidate| candidate.accepted_fragments)
            .sum(),
        decision_status: outcome.status.clone(),
    }
}

#[cfg(test)]
#[path = "../../tests/testutil/mod.rs"]
mod testutil;

#[cfg(test)]
mod tests {
    use super::testutil;
    use super::*;
    use crate::io::{FastqPairReader, FragmentStream};
    use crate::reference::{build_reference_bundle, BuildReferenceBundleRequest};
    use crate::types::{BundleManifest, ContigEntry};
    use std::fs;
    use std::io;
    use std::path::Path;
    use std::sync::{
        atomic::{AtomicUsize, Ordering},
        Arc,
    };
    use tempfile::tempdir;

    use testutil::{
        generate_fastq_pair, generate_mini_reference, FastqPairConfig, ReadComponent,
        SyntheticSource, VirusSelector,
    };

    #[test]
    fn screen_results_match_across_thread_counts() {
        let tempdir = tempdir().expect("tempdir should be created");
        let bundle =
            prepare_reference_bundle(tempdir.path()).expect("reference bundle should be prepared");
        let calibration_dir = write_calibration_profile(
            tempdir.path().join("calibration"),
            &bundle.version,
            false,
            0.0,
        )
        .expect("calibration profile should be written");
        let (r1, r2) = generate_fastq_pair(
            &FastqPairConfig::new(200, 100, 2301)
                .with_output_dir(tempdir.path().join("fastq"))
                .with_components(vec![
                    ReadComponent::new(SyntheticSource::Human, 0.95),
                    ReadComponent::new(
                        SyntheticSource::Virus(VirusSelector::Accession(
                            "NC_SYNTHV1.1".to_string(),
                        )),
                        0.05,
                    ),
                ]),
        )
        .expect("FASTQ should be generated");

        let outcome_single = run_screen(&ScreenArgs {
            input: vec![r1.clone(), r2.clone()],
            reference_bundle: bundle.bundle_dir.clone(),
            calibration_profile: calibration_dir.clone(),
            negative_control: None,
            out: tempdir.path().join("report-single-thread"),
            mode: ScreenMode::Representative,
            threads: 1,
        })
        .expect("single-thread screen should succeed");

        let outcome_parallel = run_screen(&ScreenArgs {
            input: vec![r1, r2],
            reference_bundle: bundle.bundle_dir,
            calibration_profile: calibration_dir,
            negative_control: None,
            out: tempdir.path().join("report-four-thread"),
            mode: ScreenMode::Representative,
            threads: 4,
        })
        .expect("four-thread screen should succeed");

        assert_eq!(outcome_single, outcome_parallel);
    }

    #[test]
    fn screen_parallel_results_are_stable_across_repeated_runs() {
        let tempdir = tempdir().expect("tempdir should be created");
        let bundle =
            prepare_reference_bundle(tempdir.path()).expect("reference bundle should be prepared");
        let calibration_dir = write_calibration_profile(
            tempdir.path().join("calibration"),
            &bundle.version,
            false,
            0.0,
        )
        .expect("calibration profile should be written");
        let (r1, r2) = generate_fastq_pair(
            &FastqPairConfig::new(240, 100, 9017)
                .with_output_dir(tempdir.path().join("fastq"))
                .with_components(vec![
                    ReadComponent::new(SyntheticSource::Human, 0.90),
                    ReadComponent::new(
                        SyntheticSource::Virus(VirusSelector::Accession(
                            "NC_SYNTHV1.1".to_string(),
                        )),
                        0.06,
                    ),
                    ReadComponent::new(
                        SyntheticSource::Virus(VirusSelector::Accession(
                            "NC_SYNTHV2.1".to_string(),
                        )),
                        0.04,
                    ),
                ]),
        )
        .expect("FASTQ should be generated");

        let mut baseline = None;
        for run_index in 0..5 {
            let outcome = run_screen(&ScreenArgs {
                input: vec![r1.clone(), r2.clone()],
                reference_bundle: bundle.bundle_dir.clone(),
                calibration_profile: calibration_dir.clone(),
                negative_control: None,
                out: tempdir.path().join(format!("report-repeat-{run_index}")),
                mode: ScreenMode::Representative,
                threads: 4,
            })
            .expect("repeated multi-thread screen should succeed");

            match &baseline {
                Some(expected) => assert_eq!(expected, &outcome),
                None => baseline = Some(outcome),
            }
        }
    }

    #[test]
    fn representative_mode_opens_input_only_once_per_run() {
        let tempdir = tempdir().expect("tempdir should be created");
        let bundle =
            prepare_reference_bundle(tempdir.path()).expect("reference bundle should be prepared");
        let calibration_dir = write_calibration_profile(
            tempdir.path().join("calibration"),
            &bundle.version,
            false,
            0.0,
        )
        .expect("calibration profile should be written");
        let (r1, r2) = generate_fastq_pair(
            &FastqPairConfig::new(240, 100, 4819)
                .with_output_dir(tempdir.path().join("fastq"))
                .with_components(vec![
                    ReadComponent::new(SyntheticSource::Human, 0.92),
                    ReadComponent::new(
                        SyntheticSource::Virus(VirusSelector::Accession(
                            "NC_SYNTHV1.1".to_string(),
                        )),
                        0.08,
                    ),
                ]),
        )
        .expect("FASTQ should be generated");
        let args = ScreenArgs {
            input: vec![r1.clone(), r2.clone()],
            reference_bundle: bundle.bundle_dir,
            calibration_profile: calibration_dir,
            negative_control: None,
            out: tempdir.path().join("report"),
            mode: ScreenMode::Representative,
            threads: 2,
        };
        let prepared =
            PreparedScreenRunner::prepare(&args).expect("screen runner should prepare cleanly");
        let fragments = FastqPairReader::open(&r1, &r2)
            .expect("FASTQ reader should open")
            .collect::<Result<Vec<_>>>()
            .expect("FASTQ fragments should collect");
        let input = CountingInput::new(fragments);
        let negative_control = NegativeControlDecisionInput::from_optional_path(None)
            .expect("negative control should be optional");
        let parallelism = RoundParallelism::new(args.threads)
            .expect("parallelism should respect configured threads");
        let dependencies = RoundDependencies {
            qc_filter: &prepared.qc_filter,
            aligner: &prepared.aligner,
            adjudicator: &prepared.adjudicator,
            decision_engine: &prepared.decision_engine,
            negative_control: &negative_control,
            parallelism: &parallelism,
        };

        let outcome = execute_representative_rounds(
            &input,
            &prepared.loaded_profile.profile.seed,
            &prepared.round_targets,
            &dependencies,
        )
        .expect("single-pass representative execution should succeed");

        assert_eq!(input.open_count(), 1);
        assert!(!outcome.rounds.is_empty());
    }

    #[test]
    fn prepare_prefers_prebuilt_mmi_and_normalizes_thread_budget_once() {
        let tempdir = tempdir().expect("tempdir should be created");
        let bundle =
            prepare_reference_bundle(tempdir.path()).expect("reference bundle should be prepared");
        let calibration_dir = write_calibration_profile(
            tempdir.path().join("calibration"),
            &bundle.version,
            false,
            0.0,
        )
        .expect("calibration profile should be written");
        let (r1, r2) = generate_fastq_pair(
            &FastqPairConfig::new(24, 100, 8821)
                .with_output_dir(tempdir.path().join("fastq"))
                .with_components(vec![
                    ReadComponent::new(SyntheticSource::Human, 0.92),
                    ReadComponent::new(
                        SyntheticSource::Virus(VirusSelector::Accession(
                            "NC_SYNTHV1.1".to_string(),
                        )),
                        0.08,
                    ),
                ]),
        )
        .expect("FASTQ should be generated");

        let prepared = PreparedScreenRunner::prepare(&ScreenArgs {
            input: vec![r1, r2],
            reference_bundle: bundle.bundle_dir,
            calibration_profile: calibration_dir,
            negative_control: None,
            out: tempdir.path().join("report"),
            mode: ScreenMode::Representative,
            threads: 64,
        })
        .expect("screen runner should prepare cleanly");

        assert!(prepared.aligner.uses_prebuilt_index());
        assert!(prepared
            .aligner
            .reference_input()
            .ends_with("index/minimap2/composite.mmi"));
        assert_eq!(prepared.thread_budget, 16);
        assert_eq!(prepared.aligner.thread_budget(), prepared.thread_budget);
    }

    #[test]
    fn prepared_runner_reuses_aligner_across_multiple_runs() {
        let tempdir = tempdir().expect("tempdir should be created");
        let bundle =
            prepare_reference_bundle(tempdir.path()).expect("reference bundle should be prepared");
        let calibration_dir = write_calibration_profile(
            tempdir.path().join("calibration"),
            &bundle.version,
            false,
            0.0,
        )
        .expect("calibration profile should be written");
        let (r1, r2) = generate_fastq_pair(
            &FastqPairConfig::new(128, 100, 7713)
                .with_output_dir(tempdir.path().join("fastq"))
                .with_components(vec![
                    ReadComponent::new(SyntheticSource::Human, 0.95),
                    ReadComponent::new(
                        SyntheticSource::Virus(VirusSelector::Accession(
                            "NC_SYNTHV1.1".to_string(),
                        )),
                        0.05,
                    ),
                ]),
        )
        .expect("FASTQ should be generated");
        let args = ScreenArgs {
            input: vec![r1, r2],
            reference_bundle: bundle.bundle_dir,
            calibration_profile: calibration_dir,
            negative_control: None,
            out: tempdir.path().join("report-first"),
            mode: ScreenMode::Representative,
            threads: 4,
        };
        let prepared =
            PreparedScreenRunner::prepare(&args).expect("screen runner should prepare cleanly");
        let shared_aligner_id = prepared.aligner.shared_instance_id();

        let first = prepared
            .run_without_report(&args)
            .expect("first prepared run should succeed");
        let mut second_args = args.clone();
        second_args.out = tempdir.path().join("report-second");
        let second = prepared
            .run_without_report(&second_args)
            .expect("second prepared run should reuse the same aligner");

        assert_eq!(prepared.aligner.shared_instance_id(), shared_aligner_id);
        assert_eq!(first, second);
    }

    #[test]
    fn sampling_round_count_is_bounded_before_coordinator_allocation() {
        let rounds = (1..=65).collect::<Vec<u64>>();

        let error = active_round_targets(&rounds, rounds.len() as u64)
            .expect_err("too many active rounds should be rejected");

        assert!(error
            .to_string()
            .contains("supports at most 64 active sampling rounds"));
    }

    #[derive(Clone)]
    struct CountingInput {
        fragments: Vec<IoFragmentRecord>,
        opens: Arc<AtomicUsize>,
    }

    impl CountingInput {
        fn new(fragments: Vec<IoFragmentRecord>) -> Self {
            Self {
                fragments,
                opens: Arc::new(AtomicUsize::new(0)),
            }
        }

        fn open_count(&self) -> usize {
            self.opens.load(Ordering::SeqCst)
        }
    }

    impl FragmentReaderFactory for CountingInput {
        fn open_reader(&self) -> Result<FragmentStream> {
            self.opens.fetch_add(1, Ordering::SeqCst);
            Ok(Box::new(self.fragments.clone().into_iter().map(Ok)))
        }

        fn kind_label(&self) -> &'static str {
            "fastq"
        }
    }

    struct ReferenceBundleFixture {
        bundle_dir: PathBuf,
        version: String,
    }

    fn prepare_reference_bundle(base_dir: &Path) -> io::Result<ReferenceBundleFixture> {
        let mini_reference_dir = generate_mini_reference()?;
        let manifest_path = mini_reference_dir.join("manifest.json");
        let fasta_path = mini_reference_dir.join("mini_reference.fa");
        let manifest: BundleManifest =
            serde_json::from_str(&fs::read_to_string(&manifest_path)?).map_err(io_other)?;
        let sequences = read_fasta_records(&fasta_path)?;

        let inputs_dir = base_dir.join("reference-inputs");
        fs::create_dir_all(&inputs_dir)?;
        let host_fasta = inputs_dir.join("host.fa");
        let virus_fasta = inputs_dir.join("virus.fa");
        let decoy_fasta = inputs_dir.join("decoy.fa");
        let manifest_out = inputs_dir.join("manifest.json");
        let taxonomy = inputs_dir.join("taxonomy.tsv");

        write_group_fasta(&host_fasta, &manifest.0, &sequences, "human")?;
        write_group_fasta(&virus_fasta, &manifest.0, &sequences, "virus")?;
        write_group_fasta(&decoy_fasta, &manifest.0, &sequences, "decoy")?;
        fs::copy(&manifest_path, &manifest_out)?;
        write_taxonomy(&taxonomy, &manifest.0)?;

        let bundle_dir = base_dir.join("reference-bundle");
        let outcome = build_reference_bundle(&BuildReferenceBundleRequest {
            host_fasta,
            virus_fasta,
            decoy_fasta: Some(decoy_fasta),
            manifest: manifest_out,
            taxonomy,
            out_dir: bundle_dir.clone(),
        })
        .map_err(io_other)?;

        Ok(ReferenceBundleFixture {
            bundle_dir,
            version: outcome.bundle.version,
        })
    }

    fn write_calibration_profile(
        dir: PathBuf,
        reference_bundle: &str,
        negative_control_required: bool,
        max_background_ratio: f64,
    ) -> io::Result<PathBuf> {
        fs::create_dir_all(&dir)?;
        fs::write(
            dir.join("profile.toml"),
            format!(
                "profile_id = \"rvscreen_calib_test\"\nstatus = \"release_candidate\"\nreference_bundle = \"{reference_bundle}\"\nbackend = \"minimap2\"\npreset = \"sr-conservative\"\nseed = 20260420\nsupported_input = [\"fastq\", \"fastq.gz\", \"bam\", \"ubam\", \"cram\"]\nsupported_read_type = [\"illumina_pe_shortread\"]\nnegative_control_required = {negative_control_required}\n\n[sampling]\nmode = \"representative\"\nrounds = [50, 100]\nmax_rounds = 2\n\n[fragment_rules]\nmin_mapq = 0\nmin_as_diff = 0\nmax_nm = 100\nrequire_pair_consistency = true\n\n[candidate_rules]\nmin_nonoverlap_fragments = 1\nmin_breadth = 0.0\nmax_background_ratio = {max_background_ratio}\n\n[decision_rules]\ntheta_pos = 0.01\ntheta_neg = 0.0001\nallow_indeterminate = true\n"
            ),
        )?;
        Ok(dir)
    }

    fn read_fasta_records(path: &Path) -> io::Result<Vec<(String, String)>> {
        let mut records = Vec::new();
        let mut current_header: Option<String> = None;
        let mut current_sequence = String::new();

        for line in fs::read_to_string(path)?.lines() {
            if let Some(header) = line.strip_prefix('>') {
                if let Some(previous) = current_header.replace(header.to_string()) {
                    records.push((previous, std::mem::take(&mut current_sequence)));
                }
            } else if !line.trim().is_empty() {
                current_sequence.push_str(line.trim());
            }
        }

        if let Some(header) = current_header {
            records.push((header, current_sequence));
        }

        Ok(records)
    }

    fn write_group_fasta(
        path: &Path,
        entries: &[ContigEntry],
        sequences: &[(String, String)],
        group: &str,
    ) -> io::Result<()> {
        let mut body = String::new();

        for entry in entries.iter().filter(|entry| entry.group == group) {
            let sequence = sequences
                .iter()
                .find_map(|(header, sequence)| (header == &entry.contig).then_some(sequence))
                .ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::NotFound,
                        format!("missing FASTA sequence for {}", entry.contig),
                    )
                })?;
            body.push('>');
            body.push_str(&entry.contig);
            body.push('\n');
            body.push_str(sequence);
            body.push('\n');
        }

        fs::write(path, body)
    }

    fn write_taxonomy(path: &Path, entries: &[ContigEntry]) -> io::Result<()> {
        let mut rows = String::from("taxid\tname\n");
        let mut seen = std::collections::BTreeSet::new();

        for entry in entries {
            if seen.insert((entry.taxid, entry.virus_name.clone())) {
                rows.push_str(&format!("{}\t{}\n", entry.taxid, entry.virus_name));
            }
        }

        fs::write(path, rows)
    }

    fn io_other(error: impl ToString) -> io::Error {
        io::Error::other(error.to_string())
    }
}

fn infer_sample_id(inputs: &[PathBuf]) -> String {
    inputs
        .first()
        .and_then(|path| path.file_name())
        .and_then(|name| name.to_str())
        .map(normalize_sample_name)
        .filter(|value| !value.is_empty())
        .unwrap_or_else(|| "sample".to_string())
}

fn normalize_sample_name(file_name: &str) -> String {
    let mut value = file_name.to_string();
    for suffix in [
        ".fastq.gz",
        ".fq.gz",
        ".fastq",
        ".fq",
        ".ubam",
        ".bam",
        ".cram",
    ] {
        if let Some(stripped) = value.strip_suffix(suffix) {
            value = stripped.to_string();
            break;
        }
    }

    for suffix in ["_R1", "_R2", "_1", "_2", ".R1", ".R2"] {
        if let Some(stripped) = value.strip_suffix(suffix) {
            value = stripped.to_string();
            break;
        }
    }

    value
}
