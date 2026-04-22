mod parallel;

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
    BamFragmentReader, CramFragmentReader, FastqPairReader, FragmentRecord as IoFragmentRecord,
};
use crate::pipeline::parallel::RoundParallelism;
use crate::qc::{FragmentRecord as QcFragmentRecord, QcFilter, QcStats};
use crate::report::ReportWriter;
use crate::sampling::{RepresentativeSampler, RoundManager, Sampler, StreamingSampler};
use crate::types::{CandidateCall, RoundRecord, RunManifest, SampleSummary, StopReason};
use std::path::{Path, PathBuf};

const SUPPORTED_BACKEND: &str = "minimap2";
const SUPPORTED_PRESET: &str = "sr-conservative";
const COMPOSITE_FASTA: &str = "composite.fa";

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
    input: InputSpec,
    round_targets: Vec<u64>,
    qc_filter: QcFilter,
    adjudicator: FragmentAdjudicator,
    decision_engine: DecisionEngine,
    aligner: CompetitiveAligner,
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

        let input = InputSpec::from_cli_inputs(&args.input, &loaded_bundle)?;
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
        let aligner = CompetitiveAligner::new_with_threads(
            &args.reference_bundle,
            CompetitivePreset::SrConservative,
            args.threads,
        )?;

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
        })
    }

    pub fn run(&self, args: &ScreenArgs) -> Result<ScreenRunOutcome> {
        self.execute(args, true)
    }

    pub fn run_without_report(&self, args: &ScreenArgs) -> Result<ScreenRunOutcome> {
        self.execute(args, false)
    }

    fn execute(&self, args: &ScreenArgs, write_report: bool) -> Result<ScreenRunOutcome> {
        let negative_control =
            NegativeControlDecisionInput::from_optional_path(args.negative_control.as_deref())?;
        let parallelism = RoundParallelism::new(args.threads)?;
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
            input_files: args
                .input
                .iter()
                .map(|path| path.display().to_string())
                .collect(),
        };

        if write_report {
            ReportWriter::write(
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
enum InputSpec {
    FastqPair {
        r1: PathBuf,
        r2: PathBuf,
    },
    Bam {
        path: PathBuf,
    },
    Cram {
        path: PathBuf,
        reference_fasta: PathBuf,
    },
}

impl InputSpec {
    fn from_cli_inputs(inputs: &[PathBuf], bundle: &LoadedReferenceBundle) -> Result<Self> {
        match inputs {
            [r1, r2] => {
                if !is_fastq_like(r1) || !is_fastq_like(r2) {
                    return Err(RvScreenError::validation(
                        "input",
                        format!(
                            "two-input mode currently supports only paired FASTQ files; got `{}` and `{}`",
                            r1.display(),
                            r2.display()
                        ),
                    ));
                }

                Ok(Self::FastqPair {
                    r1: r1.clone(),
                    r2: r2.clone(),
                })
            }
            [path] if is_bam_like(path) => Ok(Self::Bam { path: path.clone() }),
            [path] if is_cram_like(path) => Ok(Self::Cram {
                path: path.clone(),
                reference_fasta: bundle.bundle_dir.join(COMPOSITE_FASTA),
            }),
            [path] if is_fastq_like(path) => Err(RvScreenError::validation(
                "input",
                format!(
                    "single-file interleaved FASTQ input `{}` is not yet wired in Task 18; provide two paired FASTQ files, one BAM/uBAM, or one CRAM",
                    path.display()
                ),
            )),
            [path] => Err(RvScreenError::validation(
                "input",
                format!(
                    "unsupported single input `{}`; supported inputs are paired FASTQ, BAM/uBAM, or CRAM",
                    path.display()
                ),
            )),
            [] => Err(RvScreenError::validation(
                "input",
                "at least one input file is required",
            )),
            _ => Err(RvScreenError::validation(
                "input",
                "screen currently accepts exactly two paired FASTQs or exactly one BAM/uBAM/CRAM input",
            )),
        }
    }

    fn open_reader(&self) -> Result<Box<dyn Iterator<Item = Result<IoFragmentRecord>>>> {
        match self {
            Self::FastqPair { r1, r2 } => Ok(Box::new(FastqPairReader::open(r1, r2)?)),
            Self::Bam { path } => Ok(Box::new(BamFragmentReader::open(path)?)),
            Self::Cram {
                path,
                reference_fasta,
            } => Ok(Box::new(CramFragmentReader::open(
                path,
                Some(reference_fasta),
            )?)),
        }
    }

    fn kind_label(&self) -> &'static str {
        match self {
            Self::FastqPair { r1, .. } if is_gzip_path(r1) => "fastq.gz",
            Self::FastqPair { .. } => "fastq",
            Self::Bam { path } if is_ubam_like(path) => "ubam",
            Self::Bam { .. } => "bam",
            Self::Cram { .. } => "cram",
        }
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

enum RoundSampling<'a> {
    Representative {
        sampler: &'a RepresentativeSampler,
        round_index: usize,
    },
    Streaming {
        sampler: &'a StreamingSampler,
    },
}

impl RoundSampling<'_> {
    fn should_accept(&self, fragment_key: &str, qc_passing_index: usize) -> bool {
        match self {
            Self::Representative {
                sampler,
                round_index,
            } => sampler.should_accept(fragment_key, *round_index),
            Self::Streaming { sampler } => sampler.should_accept(fragment_key, qc_passing_index),
        }
    }

    fn sampled_limit(&self) -> Option<u64> {
        match self {
            Self::Representative { .. } => None,
            Self::Streaming { sampler } => Some(sampler.max_fragments() as u64),
        }
    }
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
    input: &InputSpec,
    seed: &u64,
    round_targets: &[u64],
    dependencies: &RoundDependencies<'_>,
) -> Result<PipelineExecution> {
    let sampler = RepresentativeSampler::new(*seed, round_targets.to_vec())?;
    let mut round_manager = RoundManager::new(round_targets.to_vec())?;
    let mut qc_baseline = None;
    let mut final_candidates = Vec::new();
    let mut rounds = Vec::new();
    let mut final_outcome = None;

    while round_manager.advance_round().is_some() {
        let round_index = round_manager
            .current_round()
            .expect("round manager should expose current round after advance");
        let execution = execute_round(
            input,
            RoundSampling::Representative {
                sampler: &sampler,
                round_index,
            },
            dependencies,
        )?;
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
        round_manager.mark_current_round_complete();

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
        RvScreenError::validation("sampling.rounds", "no representative rounds were executed")
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

fn execute_streaming_rounds(
    input: &InputSpec,
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
        let execution = execute_round(
            input,
            RoundSampling::Streaming { sampler: &sampler },
            dependencies,
        )?;
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
    input: &InputSpec,
    sampling: RoundSampling<'_>,
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

            let should_accept = sampling.should_accept(&fragment.fragment_key, qc_passing_index);
            qc_passing_index = qc_passing_index.saturating_add(1);
            if !should_accept {
                continue;
            }

            sampled_fragments = sampled_fragments.saturating_add(1);
            executor.submit(submission_sequence, fragment)?;
            submission_sequence = submission_sequence.saturating_add(1);

            if sampling
                .sampled_limit()
                .is_some_and(|limit| sampled_fragments >= limit)
            {
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
    use crate::reference::{build_reference_bundle, BuildReferenceBundleRequest};
    use crate::types::{BundleManifest, ContigEntry};
    use std::fs;
    use std::io;
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

fn is_fastq_like(path: &Path) -> bool {
    let path = path.to_string_lossy().to_ascii_lowercase();
    path.ends_with(".fastq")
        || path.ends_with(".fq")
        || path.ends_with(".fastq.gz")
        || path.ends_with(".fq.gz")
}

fn is_bam_like(path: &Path) -> bool {
    let path = path.to_string_lossy().to_ascii_lowercase();
    path.ends_with(".bam") || path.ends_with(".ubam")
}

fn is_ubam_like(path: &Path) -> bool {
    path.to_string_lossy()
        .to_ascii_lowercase()
        .ends_with(".ubam")
}

fn is_cram_like(path: &Path) -> bool {
    path.to_string_lossy()
        .to_ascii_lowercase()
        .ends_with(".cram")
}

fn is_gzip_path(path: &Path) -> bool {
    path.to_string_lossy().to_ascii_lowercase().ends_with(".gz")
}
