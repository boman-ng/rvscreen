mod parallel;
mod spool;

use crate::adjudicate::FragmentAdjudicator;
use crate::aggregate::CandidateAggregator;
use crate::align::competitive::resolve_bundle_inputs;
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
    fastq::{BorrowedFragmentRecord, FastqPairReader},
    FragmentReaderFactory, FragmentRecord as IoFragmentRecord, ScreenInput,
};
use crate::pipeline::parallel::RoundParallelism;
use crate::pipeline::spool::RepresentativeRoundSpool;
use crate::qc::{
    CheapQcPass, CheapQcResult, FragmentRecord as QcFragmentRecord, QcFilter, QcFilterReason,
    QcStats,
};
use crate::report::write_report_bundle;
use crate::sampling::{
    resolve_proportional_round_targets, ProportionalRoundResolution, RankedRepresentative,
    RepresentativeRetention, RepresentativeSampler, RepresentativeSamplingWarning, Sampler,
    StreamingSampler, PROPORTIONAL_MAX_EFFECTIVE_COUNT,
};
use crate::types::{
    CandidateCall, RoundRecord, RunManifest, SampleSummary, SamplingConfig, SamplingRoundMode,
    SamplingRoundPlanReport, SamplingWarningReport, StopReason,
};
use crossbeam_channel::{bounded, unbounded, Receiver};
use std::collections::BTreeMap;
use std::path::PathBuf;
use std::time::{Duration, Instant};

const SUPPORTED_BACKEND: &str = "minimap2";
const SUPPORTED_PRESET: &str = "sr-conservative";
const MAX_ACTIVE_SAMPLING_ROUNDS: usize = 64;
const FASTQ_RETENTION_WORK_QUEUE_MULTIPLIER: usize = 2;
pub const PREPARE_TOTAL_STAGE_LABEL: &str = "prepare_total";
pub const MINIMAP2_INDEX_LOAD_STAGE_LABEL: &str = "minimap2_index_load";
pub const RETENTION_SCAN_STAGE_LABEL: &str = "retention_scan";
pub const RETAINED_ALIGNMENT_STAGE_LABEL: &str = "retained_alignment";
pub const REPORTING_STAGE_LABEL: &str = "reporting";
pub const REPRESENTATIVE_STARTUP_STAGE_LABELS: [&str; 5] = [
    PREPARE_TOTAL_STAGE_LABEL,
    MINIMAP2_INDEX_LOAD_STAGE_LABEL,
    RETENTION_SCAN_STAGE_LABEL,
    RETAINED_ALIGNMENT_STAGE_LABEL,
    REPORTING_STAGE_LABEL,
];

#[derive(Debug, Clone, PartialEq)]
pub struct ScreenRunOutcome {
    pub summary: SampleSummary,
    pub candidates: Vec<CandidateCall>,
    pub rounds: Vec<RoundRecord>,
    pub manifest: RunManifest,
    pub perf_metrics: ScreenPerfMetrics,
    pub sampling_warnings: Vec<RepresentativeSamplingWarning>,
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct ScreenPerfMetrics {
    pub input_fragments: u64,
    pub qc_passing_fragments: u64,
    pub representative_sampled_final: u64,
    pub representative_round_sampled: u64,
    pub minimap2_calls: u64,
    pub io_wall_duration: Duration,
    pub qc_hash_wall_duration: Duration,
    pub align_wall_duration: Duration,
    pub aggregate_wall_duration: Duration,
    pub executor_queue_wait_duration: Duration,
    pub peak_heap_bytes: u64,
    pub aggregate_candidate_count: u64,
    pub stage_timings: Vec<ScreenStageTiming>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ScreenStageTiming {
    pub label: String,
    pub wall_duration: Duration,
}

impl ScreenStageTiming {
    pub fn new(label: impl Into<String>, wall_duration: Duration) -> Self {
        Self {
            label: label.into(),
            wall_duration,
        }
    }

    pub fn wall_millis(&self) -> f64 {
        self.wall_duration.as_secs_f64() * 1_000.0
    }
}

impl ScreenPerfMetrics {
    fn push_stage_timing(&mut self, label: &'static str, wall_duration: Duration) {
        self.stage_timings
            .push(ScreenStageTiming::new(label, wall_duration));
    }
}

#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct SamplingRoundOverrides<'a> {
    pub rounds: Option<&'a [u64]>,
    pub round_proportions: Option<&'a [f64]>,
}

impl<'a> SamplingRoundOverrides<'a> {
    pub fn from_screen_args(args: &'a ScreenArgs) -> Self {
        Self {
            rounds: args.rounds.as_deref(),
            round_proportions: args.round_proportions.as_deref(),
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum ResolvedSamplingRoundPlan {
    Absolute(Vec<u64>),
    Proportional(Vec<f64>),
}

impl ResolvedSamplingRoundPlan {
    pub fn mode(&self) -> SamplingRoundMode {
        match self {
            Self::Absolute(_) => SamplingRoundMode::Absolute,
            Self::Proportional(_) => SamplingRoundMode::Proportional,
        }
    }
}

pub struct PreparedScreenRunner {
    loaded_bundle: LoadedReferenceBundle,
    loaded_profile: LoadedProfile,
    benchmark_gates: Option<CalibrationReleaseGate>,
    input: ScreenInput,
    round_plan: ResolvedSamplingRoundPlan,
    qc_filter: QcFilter,
    adjudicator: FragmentAdjudicator,
    decision_engine: DecisionEngine,
    aligner: CompetitiveAligner,
    thread_budget: usize,
    stage_timings: Vec<ScreenStageTiming>,
}

pub fn run_screen(args: &ScreenArgs) -> Result<ScreenRunOutcome> {
    PreparedScreenRunner::prepare(args)?.run(args)
}

impl PreparedScreenRunner {
    pub fn prepare(args: &ScreenArgs) -> Result<Self> {
        let prepare_started = Instant::now();
        args.validate_sampling_overrides()
            .map_err(|error| RvScreenError::validation("sampling overrides", error))?;
        let loaded_bundle = load_reference_bundle(&args.reference_bundle)?;
        let loaded_profile = load_profile(&args.calibration_profile)?;
        let benchmark_gates = load_benchmark_gates(&loaded_profile.profile_dir)?;
        validate_profile_compatibility(&loaded_bundle, &loaded_profile.profile)?;

        let input = ScreenInput::from_cli_inputs(&args.input, &loaded_bundle.bundle_dir)?;
        validate_supported_input(&loaded_profile.profile.supported_input, input.kind_label())?;

        let round_plan = resolve_sampling_round_plan(
            &loaded_profile.profile.sampling,
            SamplingRoundOverrides::from_screen_args(args),
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
        let minimap2_index_load_duration = aligner.minimap2_index_load_duration();
        let prepare_total_duration = prepare_started.elapsed();
        let stage_timings = vec![
            ScreenStageTiming::new(PREPARE_TOTAL_STAGE_LABEL, prepare_total_duration),
            ScreenStageTiming::new(
                MINIMAP2_INDEX_LOAD_STAGE_LABEL,
                minimap2_index_load_duration,
            ),
        ];

        Ok(Self {
            loaded_bundle,
            loaded_profile,
            benchmark_gates,
            input,
            round_plan,
            qc_filter,
            adjudicator,
            decision_engine,
            aligner,
            thread_budget,
            stage_timings,
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

        let execution = match (args.mode, &self.round_plan) {
            (ScreenMode::Representative, ResolvedSamplingRoundPlan::Absolute(round_targets)) => {
                execute_representative_rounds(
                    &self.input,
                    &self.loaded_profile.profile.seed,
                    round_targets,
                    &round_dependencies,
                )?
            }
            (
                ScreenMode::Representative,
                ResolvedSamplingRoundPlan::Proportional(round_proportions),
            ) => execute_proportional_representative_rounds(
                &self.input,
                &self.loaded_profile.profile.seed,
                round_proportions,
                args.allow_sampling_threshold_override,
                &round_dependencies,
            )?,
            (ScreenMode::Streaming, ResolvedSamplingRoundPlan::Absolute(round_targets)) => {
                execute_streaming_rounds(&self.input, round_targets, &round_dependencies)?
            }
            (ScreenMode::Streaming, ResolvedSamplingRoundPlan::Proportional(_)) => {
                return Err(RvScreenError::validation(
                    "sampling.round_proportions",
                    "proportional sampling rounds are only supported for representative mode",
                ));
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
            sampling_round_plan: execution.sampling_round_plan.clone(),
            negative_control: negative_control
                .manifest(self.loaded_profile.profile.negative_control_required),
            input_files: args
                .input
                .iter()
                .map(|path| path.display().to_string())
                .collect(),
        };

        let mut reporting_duration = Duration::default();
        if write_report {
            let reporting_started = Instant::now();
            write_report_bundle(
                &args.out,
                &summary,
                &execution.candidates,
                &execution.rounds,
                &manifest,
            )?;
            reporting_duration = reporting_started.elapsed();
        }

        let mut perf_metrics = execution.perf_metrics;
        let mut stage_timings = self.stage_timings.clone();
        stage_timings.append(&mut perf_metrics.stage_timings);
        perf_metrics.stage_timings = stage_timings;
        perf_metrics.push_stage_timing(REPORTING_STAGE_LABEL, reporting_duration);

        Ok(ScreenRunOutcome {
            summary,
            candidates: execution.candidates,
            rounds: execution.rounds,
            manifest,
            perf_metrics: with_runtime_memory_metrics(perf_metrics),
            sampling_warnings: execution.sampling_warnings,
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
    perf_metrics: ScreenPerfMetrics,
    sampling_warnings: Vec<RepresentativeSamplingWarning>,
    sampling_round_plan: Option<SamplingRoundPlanReport>,
}

#[derive(Debug, Clone)]
struct RoundExecution {
    qc_stats: QcStats,
    sampled_fragments: u64,
    candidates: Vec<CandidateCall>,
    perf_metrics: ScreenPerfMetrics,
}

struct RepresentativePassExecution {
    qc_stats: QcStats,
    round_spool: RepresentativeRoundSpool,
    perf_metrics: ScreenPerfMetrics,
}

struct ProportionalRepresentativePassExecution {
    qc_stats: QcStats,
    round_spool: RepresentativeRoundSpool,
    perf_metrics: ScreenPerfMetrics,
    effective_round_targets: Vec<u64>,
    sampling_warnings: Vec<RepresentativeSamplingWarning>,
    sampling_round_plan: SamplingRoundPlanReport,
}

struct FastqRepresentativeRetentionMetrics {
    qc_stats: QcStats,
    retained_fragments: Vec<RankedRepresentative<IoFragmentRecord>>,
    io_wall_duration: Duration,
    qc_hash_wall_duration: Duration,
}

#[derive(Debug, Clone, Copy)]
enum FastqRetentionMode {
    Absolute,
    Proportional,
}

struct FastqRetentionBatch {
    sequence: usize,
    fragments: Vec<IoFragmentRecord>,
}

struct FastqRetentionBatchResult {
    sequence: usize,
    entries: Vec<FastqRetentionBatchEntry>,
    qc_hash_wall_duration: Duration,
}

enum FastqRetentionBatchEntry {
    CheapFiltered(QcFilterReason),
    CheapPassed {
        fragment: IoFragmentRecord,
        cheap_pass: CheapQcPass,
        rank_hash: u64,
    },
    Filtered(QcFilterReason),
    Passed {
        fragment: IoFragmentRecord,
        rank_hash: u64,
    },
}

struct RepresentativeBatchExecution {
    round_spool: RepresentativeRoundSpool,
    executor_queue_wait_duration: Duration,
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
    if let Some(reader) = input.open_fastq_pair_reader() {
        return execute_representative_fastq_rounds(reader?, seed, round_targets, dependencies);
    }

    let sampling_round_plan = absolute_sampling_round_plan(round_targets);
    let sampler = RepresentativeSampler::new(*seed, round_targets.to_vec())?;

    if sampler.round_count() == 1 {
        let execution = execute_representative_single_round_pass(input, &sampler, dependencies)?;
        let background_adjusted =
            apply_negative_control(&execution.candidates, dependencies.negative_control);
        let evaluated = dependencies
            .decision_engine
            .evaluate_candidates(&background_adjusted);
        let outcome = dependencies
            .decision_engine
            .decide(&evaluated, &DecisionContext::new(1, 1));
        let rounds = vec![round_record(
            &evaluated,
            execution.sampled_fragments,
            &outcome,
        )];

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
            perf_metrics: execution.perf_metrics,
            sampling_warnings: Vec::new(),
            sampling_round_plan: Some(sampling_round_plan),
        });
    }

    let RepresentativePassExecution {
        qc_stats,
        round_spool,
        perf_metrics: representative_perf_metrics,
    } = execute_representative_pass(input, &sampler, dependencies)?;
    let mut cumulative_aggregator = CandidateAggregator::new();
    let mut cumulative_sampled_fragments = 0u64;
    let mut final_candidates = Vec::new();
    let mut rounds = Vec::new();
    let mut final_outcome = None;

    let aggregate_started = Instant::now();
    for (round_index, contribution) in round_spool
        .into_round_contributions()
        .into_iter()
        .enumerate()
    {
        cumulative_sampled_fragments =
            cumulative_sampled_fragments.saturating_add(contribution.sampled_fragments);
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
            let candidate_count = final_candidates.len();
            let mut perf_metrics = finalize_aggregate_perf_metrics(
                representative_perf_metrics.clone(),
                candidate_count,
            );
            perf_metrics.aggregate_wall_duration = aggregate_started.elapsed();
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
                perf_metrics,
                sampling_warnings: Vec::new(),
                sampling_round_plan: Some(sampling_round_plan.clone()),
            });
        }
    }

    let final_outcome = final_outcome.ok_or_else(|| {
        RvScreenError::validation("sampling.rounds", "no representative rounds were executed")
    })?;
    let final_candidate_count = final_candidates.len();
    let mut perf_metrics =
        finalize_aggregate_perf_metrics(representative_perf_metrics, final_candidate_count);
    perf_metrics.aggregate_wall_duration = aggregate_started.elapsed();
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
        perf_metrics,
        sampling_warnings: Vec::new(),
        sampling_round_plan: Some(sampling_round_plan),
    })
}

fn execute_representative_fastq_rounds(
    mut reader: FastqPairReader,
    seed: &u64,
    round_targets: &[u64],
    dependencies: &RoundDependencies<'_>,
) -> Result<PipelineExecution> {
    let sampling_round_plan = absolute_sampling_round_plan(round_targets);
    let sampler = RepresentativeSampler::new(*seed, round_targets.to_vec())?;

    if sampler.round_count() == 1 {
        let execution =
            execute_representative_single_round_fastq_pass(&mut reader, &sampler, dependencies)?;
        let background_adjusted =
            apply_negative_control(&execution.candidates, dependencies.negative_control);
        let evaluated = dependencies
            .decision_engine
            .evaluate_candidates(&background_adjusted);
        let outcome = dependencies
            .decision_engine
            .decide(&evaluated, &DecisionContext::new(1, 1));
        let rounds = vec![round_record(
            &evaluated,
            execution.sampled_fragments,
            &outcome,
        )];

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
            perf_metrics: execution.perf_metrics,
            sampling_warnings: Vec::new(),
            sampling_round_plan: Some(sampling_round_plan),
        });
    }

    let RepresentativePassExecution {
        qc_stats,
        round_spool,
        perf_metrics: representative_perf_metrics,
    } = execute_fastq_representative_pass(&mut reader, &sampler, dependencies)?;
    let mut cumulative_aggregator = CandidateAggregator::new();
    let mut cumulative_sampled_fragments = 0u64;
    let mut final_candidates = Vec::new();
    let mut rounds = Vec::new();
    let mut final_outcome = None;

    let aggregate_started = Instant::now();
    for (round_index, contribution) in round_spool
        .into_round_contributions()
        .into_iter()
        .enumerate()
    {
        cumulative_sampled_fragments =
            cumulative_sampled_fragments.saturating_add(contribution.sampled_fragments);
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
            let candidate_count = final_candidates.len();
            let mut perf_metrics = finalize_aggregate_perf_metrics(
                representative_perf_metrics.clone(),
                candidate_count,
            );
            perf_metrics.aggregate_wall_duration = aggregate_started.elapsed();
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
                perf_metrics,
                sampling_warnings: Vec::new(),
                sampling_round_plan: Some(sampling_round_plan.clone()),
            });
        }
    }

    let final_outcome = final_outcome.ok_or_else(|| {
        RvScreenError::validation("sampling.rounds", "no representative rounds were executed")
    })?;
    let final_candidate_count = final_candidates.len();
    let mut perf_metrics =
        finalize_aggregate_perf_metrics(representative_perf_metrics, final_candidate_count);
    perf_metrics.aggregate_wall_duration = aggregate_started.elapsed();
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
        perf_metrics,
        sampling_warnings: Vec::new(),
        sampling_round_plan: Some(sampling_round_plan),
    })
}

fn execute_proportional_representative_rounds(
    input: &dyn FragmentReaderFactory,
    seed: &u64,
    round_proportions: &[f64],
    allow_sampling_threshold_override: bool,
    dependencies: &RoundDependencies<'_>,
) -> Result<PipelineExecution> {
    if let Some(reader) = input.open_fastq_pair_reader() {
        return execute_proportional_representative_fastq_rounds(
            reader?,
            seed,
            round_proportions,
            allow_sampling_threshold_override,
            dependencies,
        );
    }

    let sampler = RepresentativeSampler::for_retention_capacity(
        *seed,
        proportional_retention_capacity(allow_sampling_threshold_override)?,
    );
    let execution = execute_proportional_representative_pass(
        input,
        &sampler,
        round_proportions,
        allow_sampling_threshold_override,
        dependencies,
    )?;

    finalize_representative_spool_execution(execution, dependencies)
}

fn execute_proportional_representative_fastq_rounds(
    mut reader: FastqPairReader,
    seed: &u64,
    round_proportions: &[f64],
    allow_sampling_threshold_override: bool,
    dependencies: &RoundDependencies<'_>,
) -> Result<PipelineExecution> {
    let sampler = RepresentativeSampler::for_retention_capacity(
        *seed,
        proportional_retention_capacity(allow_sampling_threshold_override)?,
    );
    let execution = execute_fastq_proportional_representative_pass(
        &mut reader,
        &sampler,
        round_proportions,
        allow_sampling_threshold_override,
        dependencies,
    )?;

    finalize_representative_spool_execution(execution, dependencies)
}

fn finalize_representative_spool_execution(
    execution: ProportionalRepresentativePassExecution,
    dependencies: &RoundDependencies<'_>,
) -> Result<PipelineExecution> {
    let ProportionalRepresentativePassExecution {
        qc_stats,
        round_spool,
        perf_metrics: representative_perf_metrics,
        effective_round_targets,
        sampling_warnings,
        sampling_round_plan,
    } = execution;
    let mut cumulative_aggregator = CandidateAggregator::new();
    let mut cumulative_sampled_fragments = 0u64;
    let mut final_candidates = Vec::new();
    let mut rounds = Vec::new();
    let mut final_outcome = None;

    let aggregate_started = Instant::now();
    for (round_index, contribution) in round_spool
        .into_round_contributions()
        .into_iter()
        .enumerate()
    {
        cumulative_sampled_fragments =
            cumulative_sampled_fragments.saturating_add(contribution.sampled_fragments);
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
            &DecisionContext::new(
                (round_index + 1) as u64,
                effective_round_targets.len() as u64,
            ),
        );
        rounds.push(round_record(
            &evaluated,
            cumulative_sampled_fragments,
            &outcome,
        ));
        final_candidates = evaluated;
        final_outcome = Some(outcome.clone());

        if outcome.should_stop() {
            let candidate_count = final_candidates.len();
            let mut perf_metrics = finalize_aggregate_perf_metrics(
                representative_perf_metrics.clone(),
                candidate_count,
            );
            perf_metrics.aggregate_wall_duration = aggregate_started.elapsed();
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
                perf_metrics,
                sampling_warnings,
                sampling_round_plan: Some(sampling_round_plan),
            });
        }
    }

    let final_outcome = final_outcome.ok_or_else(|| {
        RvScreenError::validation(
            "sampling.round_proportions",
            "no proportional representative rounds were executed",
        )
    })?;
    let final_candidate_count = final_candidates.len();
    let mut perf_metrics =
        finalize_aggregate_perf_metrics(representative_perf_metrics, final_candidate_count);
    perf_metrics.aggregate_wall_duration = aggregate_started.elapsed();
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
        perf_metrics,
        sampling_warnings,
        sampling_round_plan: Some(sampling_round_plan),
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
                perf_metrics: execution.perf_metrics,
                sampling_warnings: Vec::new(),
                sampling_round_plan: None,
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
        perf_metrics: ScreenPerfMetrics::default(),
        sampling_warnings: Vec::new(),
        sampling_round_plan: None,
    })
}

fn execute_round(
    input: &dyn FragmentReaderFactory,
    sampler: &StreamingSampler,
    dependencies: &RoundDependencies<'_>,
) -> Result<RoundExecution> {
    let io_started = Instant::now();
    let reader = input.open_reader()?;
    let mut io_wall_duration = io_started.elapsed();
    let mut qc_stats = QcStats::default();
    let mut qc_passing_index = 0usize;
    let mut sampled_fragments = 0u64;
    let mut submission_sequence = 0usize;
    let mut qc_hash_wall_duration = Duration::default();
    let mut executor_queue_wait_duration = Duration::default();

    let mut aggregator = std::thread::scope(|scope| -> Result<CandidateAggregator> {
        let executor = dependencies.parallelism.start_round_executor(
            scope,
            dependencies.aligner,
            dependencies.adjudicator,
        );

        for fragment in reader {
            let io_step_started = Instant::now();
            let fragment = fragment?;
            io_wall_duration = io_wall_duration.saturating_add(io_step_started.elapsed());

            let qc_hash_started = Instant::now();
            let qc_fragment = bridge_qc_fragment(&fragment);
            let qc_result = dependencies
                .qc_filter
                .filter_and_record(&qc_fragment, &mut qc_stats);
            qc_hash_wall_duration = qc_hash_wall_duration.saturating_add(qc_hash_started.elapsed());
            if !qc_result.is_pass() {
                continue;
            }

            let hash_started = Instant::now();
            let should_accept = sampler.should_accept(&fragment.fragment_key, qc_passing_index);
            qc_passing_index = qc_passing_index.saturating_add(1);
            qc_hash_wall_duration = qc_hash_wall_duration.saturating_add(hash_started.elapsed());
            if !should_accept {
                continue;
            }

            sampled_fragments = sampled_fragments.saturating_add(1);
            let queue_wait_started = Instant::now();
            executor.submit(submission_sequence, fragment)?;
            executor_queue_wait_duration =
                executor_queue_wait_duration.saturating_add(queue_wait_started.elapsed());
            submission_sequence = submission_sequence.saturating_add(1);

            if sampled_fragments >= sampler.max_fragments() as u64 {
                break;
            }
        }

        executor.finish()
    })?;

    aggregator.set_total_sampled_fragments(sampled_fragments);
    let aggregate_started = Instant::now();
    let candidates = aggregator.finalize(&qc_stats);
    let aggregate_duration = aggregate_started.elapsed();
    let candidate_count = candidates.len();

    Ok(RoundExecution {
        qc_stats,
        sampled_fragments,
        candidates,
        perf_metrics: finalize_aggregate_perf_metrics(
            ScreenPerfMetrics {
                input_fragments: qc_stats.total_fragments,
                qc_passing_fragments: qc_stats.passed_fragments,
                representative_sampled_final: sampled_fragments,
                representative_round_sampled: sampled_fragments,
                minimap2_calls: sampled_fragments,
                io_wall_duration,
                qc_hash_wall_duration,
                align_wall_duration: Duration::default(),
                aggregate_wall_duration: aggregate_duration,
                executor_queue_wait_duration,
                peak_heap_bytes: 0,
                aggregate_candidate_count: 0,
                stage_timings: Vec::new(),
            },
            candidate_count,
        ),
    })
}

fn execute_representative_pass(
    input: &dyn FragmentReaderFactory,
    sampler: &RepresentativeSampler,
    dependencies: &RoundDependencies<'_>,
) -> Result<RepresentativePassExecution> {
    let retention_scan_started = Instant::now();
    let io_started = Instant::now();
    let reader = input.open_reader()?;
    let mut io_wall_duration = io_started.elapsed();
    let mut qc_stats = QcStats::default();
    let mut retention = RepresentativeRetention::new(sampler.max_fragments());
    let mut qc_hash_wall_duration = Duration::default();

    for fragment in reader {
        let io_step_started = Instant::now();
        let fragment = fragment?;
        io_wall_duration = io_wall_duration.saturating_add(io_step_started.elapsed());
        let qc_hash_started = Instant::now();
        if !retain_owned_representative_fragment(
            &fragment,
            sampler,
            dependencies.qc_filter,
            &mut qc_stats,
            &mut retention,
        ) {
            qc_hash_wall_duration = qc_hash_wall_duration.saturating_add(qc_hash_started.elapsed());
            continue;
        }
        qc_hash_wall_duration = qc_hash_wall_duration.saturating_add(qc_hash_started.elapsed());
    }

    let retained_fragments = retention.into_sorted_vec();
    let retention_scan_duration = retention_scan_started.elapsed();
    let prefix_sample_sizes = sampler.prefix_sample_sizes(retained_fragments.len());
    let representative_round_sampled = prefix_last_sampled(&prefix_sample_sizes);
    let align_started = Instant::now();
    let RepresentativeBatchExecution {
        round_spool,
        executor_queue_wait_duration,
    } = execute_representative_batches(
        retained_fragments,
        prefix_sample_sizes.clone(),
        dependencies,
    )?;
    let align_duration = align_started.elapsed();

    Ok(RepresentativePassExecution {
        qc_stats,
        round_spool,
        perf_metrics: representative_perf_metrics(
            qc_stats,
            prefix_sample_sizes.last().copied().unwrap_or(0),
            representative_round_sampled,
            io_wall_duration,
            qc_hash_wall_duration,
            align_duration,
            executor_queue_wait_duration,
            retention_scan_duration,
        ),
    })
}

fn execute_fastq_representative_pass(
    reader: &mut FastqPairReader,
    sampler: &RepresentativeSampler,
    dependencies: &RoundDependencies<'_>,
) -> Result<RepresentativePassExecution> {
    let retention_scan_started = Instant::now();
    let FastqRepresentativeRetentionMetrics {
        qc_stats,
        retained_fragments,
        io_wall_duration,
        qc_hash_wall_duration,
    } = retain_fastq_representative_fragments_with_parallelism(
        reader,
        sampler,
        dependencies.qc_filter,
        dependencies.parallelism,
    )?;
    let retention_scan_duration = retention_scan_started.elapsed();
    let prefix_sample_sizes = sampler.prefix_sample_sizes(retained_fragments.len());
    let representative_round_sampled = prefix_last_sampled(&prefix_sample_sizes);
    let align_started = Instant::now();
    let RepresentativeBatchExecution {
        round_spool,
        executor_queue_wait_duration,
    } = execute_representative_batches(
        retained_fragments,
        prefix_sample_sizes.clone(),
        dependencies,
    )?;
    let align_duration = align_started.elapsed();

    Ok(RepresentativePassExecution {
        qc_stats,
        round_spool,
        perf_metrics: representative_perf_metrics(
            qc_stats,
            prefix_sample_sizes.last().copied().unwrap_or(0),
            representative_round_sampled,
            io_wall_duration,
            qc_hash_wall_duration,
            align_duration,
            executor_queue_wait_duration,
            retention_scan_duration,
        ),
    })
}

fn execute_proportional_representative_pass(
    input: &dyn FragmentReaderFactory,
    sampler: &RepresentativeSampler,
    round_proportions: &[f64],
    allow_sampling_threshold_override: bool,
    dependencies: &RoundDependencies<'_>,
) -> Result<ProportionalRepresentativePassExecution> {
    let retention_scan_started = Instant::now();
    let io_started = Instant::now();
    let reader = input.open_reader()?;
    let mut io_wall_duration = io_started.elapsed();
    let mut qc_stats = QcStats::default();
    let mut retention = RepresentativeRetention::new(sampler.max_fragments());
    let mut qc_hash_wall_duration = Duration::default();

    for fragment in reader {
        let io_step_started = Instant::now();
        let fragment = fragment?;
        io_wall_duration = io_wall_duration.saturating_add(io_step_started.elapsed());
        let qc_hash_started = Instant::now();
        retain_owned_proportional_representative_fragment(
            &fragment,
            sampler,
            dependencies.qc_filter,
            &mut qc_stats,
            &mut retention,
        );
        qc_hash_wall_duration = qc_hash_wall_duration.saturating_add(qc_hash_started.elapsed());
    }

    finish_proportional_representative_pass(
        qc_stats,
        retention.into_sorted_vec(),
        round_proportions,
        allow_sampling_threshold_override,
        io_wall_duration,
        qc_hash_wall_duration,
        retention_scan_started.elapsed(),
        dependencies,
    )
}

fn execute_fastq_proportional_representative_pass(
    reader: &mut FastqPairReader,
    sampler: &RepresentativeSampler,
    round_proportions: &[f64],
    allow_sampling_threshold_override: bool,
    dependencies: &RoundDependencies<'_>,
) -> Result<ProportionalRepresentativePassExecution> {
    let retention_scan_started = Instant::now();
    let FastqRepresentativeRetentionMetrics {
        qc_stats,
        retained_fragments,
        io_wall_duration,
        qc_hash_wall_duration,
    } = retain_fastq_proportional_representative_fragments_with_parallelism(
        reader,
        sampler,
        dependencies.qc_filter,
        dependencies.parallelism,
    )?;

    finish_proportional_representative_pass(
        qc_stats,
        retained_fragments,
        round_proportions,
        allow_sampling_threshold_override,
        io_wall_duration,
        qc_hash_wall_duration,
        retention_scan_started.elapsed(),
        dependencies,
    )
}

fn finish_proportional_representative_pass(
    qc_stats: QcStats,
    mut retained_fragments: Vec<RankedRepresentative<IoFragmentRecord>>,
    round_proportions: &[f64],
    allow_sampling_threshold_override: bool,
    io_wall_duration: Duration,
    qc_hash_wall_duration: Duration,
    retention_scan_duration: Duration,
    dependencies: &RoundDependencies<'_>,
) -> Result<ProportionalRepresentativePassExecution> {
    let resolution = resolve_proportional_round_targets(
        round_proportions,
        qc_stats.passed_fragments,
        allow_sampling_threshold_override,
    );
    let sampling_round_plan =
        proportional_sampling_round_plan(round_proportions, qc_stats.passed_fragments, &resolution);
    let prefix_sample_sizes = resolution.effective_counts.clone();
    let final_prefix_len = usize::try_from(prefix_sample_sizes.last().copied().unwrap_or(0))
        .map_err(|_| {
            RvScreenError::validation(
                "sampling.round_proportions",
                "proportional representative effective count exceeds platform usize capacity",
            )
        })?;
    retained_fragments.truncate(final_prefix_len);
    let representative_round_sampled = prefix_last_sampled(&prefix_sample_sizes);
    let align_started = Instant::now();
    let RepresentativeBatchExecution {
        round_spool,
        executor_queue_wait_duration,
    } = execute_representative_batches(
        retained_fragments,
        prefix_sample_sizes.clone(),
        dependencies,
    )?;
    let align_duration = align_started.elapsed();

    Ok(ProportionalRepresentativePassExecution {
        qc_stats,
        round_spool,
        perf_metrics: representative_perf_metrics(
            qc_stats,
            prefix_sample_sizes.last().copied().unwrap_or(0),
            representative_round_sampled,
            io_wall_duration,
            qc_hash_wall_duration,
            align_duration,
            executor_queue_wait_duration,
            retention_scan_duration,
        ),
        effective_round_targets: prefix_sample_sizes,
        sampling_warnings: resolution.warnings,
        sampling_round_plan,
    })
}

fn execute_representative_single_round_pass(
    input: &dyn FragmentReaderFactory,
    sampler: &RepresentativeSampler,
    dependencies: &RoundDependencies<'_>,
) -> Result<RoundExecution> {
    let RepresentativePassExecution {
        qc_stats,
        round_spool,
        perf_metrics,
    } = execute_representative_pass(input, sampler, dependencies)?;
    finalize_single_round_representative_execution(qc_stats, round_spool, perf_metrics)
}

fn execute_representative_single_round_fastq_pass(
    reader: &mut FastqPairReader,
    sampler: &RepresentativeSampler,
    dependencies: &RoundDependencies<'_>,
) -> Result<RoundExecution> {
    let RepresentativePassExecution {
        qc_stats,
        round_spool,
        perf_metrics,
    } = execute_fastq_representative_pass(reader, sampler, dependencies)?;
    finalize_single_round_representative_execution(qc_stats, round_spool, perf_metrics)
}

fn execute_representative_batches(
    retained_fragments: Vec<RankedRepresentative<IoFragmentRecord>>,
    prefix_sample_sizes: Vec<u64>,
    dependencies: &RoundDependencies<'_>,
) -> Result<RepresentativeBatchExecution> {
    std::thread::scope(|scope| -> Result<RepresentativeBatchExecution> {
        let executor = dependencies
            .parallelism
            .start_representative_round_executor(
                scope,
                dependencies.aligner,
                dependencies.adjudicator,
                RepresentativeRoundSpool::new(prefix_sample_sizes.clone())?,
            );
        let batch_size = dependencies.parallelism.representative_batch_size();
        let mut batch_sequence = 0usize;
        let mut batch = Vec::with_capacity(batch_size);
        let mut executor_queue_wait_duration = Duration::default();

        for (prefix_index, retained_fragment) in retained_fragments.into_iter().enumerate() {
            let prefix_index = u64::try_from(prefix_index).map_err(|_| {
                RvScreenError::validation(
                    "sampling.rounds",
                    "representative retained prefix index exceeds u64 capacity",
                )
            })?;
            let entry_round = RepresentativeRoundSpool::entry_round_for_prefix_index(
                &prefix_sample_sizes,
                prefix_index,
            )
            .ok_or_else(|| {
                RvScreenError::validation(
                    "sampling.rounds",
                    format!(
                        "representative retained prefix index {} does not map to any configured round",
                        prefix_index
                    ),
                )
            })?;
            batch.push((entry_round, retained_fragment.into_payload()));

            if batch.len() == batch_size {
                let queue_wait_started = Instant::now();
                executor.submit_batch(batch_sequence, std::mem::take(&mut batch))?;
                executor_queue_wait_duration =
                    executor_queue_wait_duration.saturating_add(queue_wait_started.elapsed());
                batch_sequence = batch_sequence.saturating_add(1);
            }
        }

        if !batch.is_empty() {
            let queue_wait_started = Instant::now();
            executor.submit_batch(batch_sequence, batch)?;
            executor_queue_wait_duration =
                executor_queue_wait_duration.saturating_add(queue_wait_started.elapsed());
        }

        Ok(RepresentativeBatchExecution {
            round_spool: executor.finish()?,
            executor_queue_wait_duration,
        })
    })
}

fn finalize_single_round_representative_execution(
    qc_stats: QcStats,
    round_spool: RepresentativeRoundSpool,
    mut perf_metrics: ScreenPerfMetrics,
) -> Result<RoundExecution> {
    let mut contributions = round_spool.into_round_contributions();
    if contributions.len() != 1 {
        return Err(RvScreenError::validation(
            "sampling.rounds",
            format!(
                "single-round representative execution expected exactly one contribution, but observed {}",
                contributions.len()
            ),
        ));
    }

    let contribution = contributions
        .pop()
        .expect("single-round representative execution should expose one contribution");
    let sampled_fragments = contribution.sampled_fragments;
    let mut aggregator = contribution.aggregator;
    aggregator.set_total_sampled_fragments(sampled_fragments);
    let aggregate_started = Instant::now();
    let candidates = aggregator.finalize(&qc_stats);
    perf_metrics.aggregate_wall_duration = aggregate_started.elapsed();
    perf_metrics.aggregate_candidate_count = candidates.len() as u64;
    perf_metrics.representative_sampled_final = sampled_fragments;
    perf_metrics.representative_round_sampled = sampled_fragments;
    perf_metrics.minimap2_calls = sampled_fragments;

    Ok(RoundExecution {
        qc_stats,
        sampled_fragments,
        candidates,
        perf_metrics,
    })
}

fn finalize_aggregate_perf_metrics(
    mut perf_metrics: ScreenPerfMetrics,
    candidate_count: usize,
) -> ScreenPerfMetrics {
    perf_metrics.aggregate_candidate_count = candidate_count as u64;
    perf_metrics
}

fn representative_perf_metrics(
    qc_stats: QcStats,
    sampled_final: u64,
    representative_round_sampled: u64,
    io_wall_duration: Duration,
    qc_hash_wall_duration: Duration,
    align_wall_duration: Duration,
    executor_queue_wait_duration: Duration,
    retention_scan_duration: Duration,
) -> ScreenPerfMetrics {
    let mut perf_metrics = ScreenPerfMetrics {
        input_fragments: qc_stats.total_fragments,
        qc_passing_fragments: qc_stats.passed_fragments,
        representative_sampled_final: sampled_final,
        representative_round_sampled,
        minimap2_calls: sampled_final,
        io_wall_duration,
        qc_hash_wall_duration,
        align_wall_duration,
        aggregate_wall_duration: Duration::default(),
        executor_queue_wait_duration,
        peak_heap_bytes: 0,
        aggregate_candidate_count: 0,
        stage_timings: Vec::new(),
    };
    perf_metrics.push_stage_timing(RETENTION_SCAN_STAGE_LABEL, retention_scan_duration);
    perf_metrics.push_stage_timing(RETAINED_ALIGNMENT_STAGE_LABEL, align_wall_duration);
    perf_metrics
}

fn with_runtime_memory_metrics(mut perf_metrics: ScreenPerfMetrics) -> ScreenPerfMetrics {
    perf_metrics.peak_heap_bytes = current_rss_bytes().unwrap_or(0);
    perf_metrics
}

fn current_rss_bytes() -> std::io::Result<u64> {
    let status = std::fs::read_to_string("/proc/self/status")?;
    status
        .lines()
        .find_map(|line| {
            line.strip_prefix("VmHWM:")
                .or_else(|| line.strip_prefix("VmRSS:"))
        })
        .and_then(|value| value.split_whitespace().next())
        .and_then(|value| value.parse::<u64>().ok())
        .map(|value| value.saturating_mul(1024))
        .ok_or_else(|| std::io::Error::other("VmRSS/VmHWM not found in /proc/self/status"))
}

fn prefix_last_sampled(prefix_sample_sizes: &[u64]) -> u64 {
    prefix_sample_sizes.last().copied().unwrap_or(0)
}

fn proportional_retention_capacity(allow_sampling_threshold_override: bool) -> Result<usize> {
    if allow_sampling_threshold_override {
        return Ok(usize::MAX);
    }

    usize::try_from(PROPORTIONAL_MAX_EFFECTIVE_COUNT).map_err(|_| {
        RvScreenError::validation(
            "sampling.round_proportions",
            "proportional representative maximum threshold exceeds platform usize capacity",
        )
    })
}

fn retain_fastq_representative_fragments(
    reader: &mut FastqPairReader,
    sampler: &RepresentativeSampler,
    qc_filter: &QcFilter,
) -> Result<FastqRepresentativeRetentionMetrics> {
    let mut qc_stats = QcStats::default();
    let mut retention = RepresentativeRetention::new(sampler.max_fragments());
    let mut qc_hash_wall_duration = Duration::default();

    let io_wall_duration = reader.for_each_borrowed_pair_with_timing(|fragment| {
        let qc_hash_started = Instant::now();
        retain_borrowed_fastq_representative_fragment(
            fragment,
            sampler,
            qc_filter,
            &mut qc_stats,
            &mut retention,
        );
        qc_hash_wall_duration = qc_hash_wall_duration.saturating_add(qc_hash_started.elapsed());
        Ok(())
    })?;

    Ok(FastqRepresentativeRetentionMetrics {
        qc_stats,
        retained_fragments: retention.into_sorted_vec(),
        io_wall_duration,
        qc_hash_wall_duration,
    })
}

fn retain_fastq_representative_fragments_with_parallelism(
    reader: &mut FastqPairReader,
    sampler: &RepresentativeSampler,
    qc_filter: &QcFilter,
    parallelism: &RoundParallelism,
) -> Result<FastqRepresentativeRetentionMetrics> {
    if parallelism.threads() <= 1 {
        return retain_fastq_representative_fragments(reader, sampler, qc_filter);
    }

    retain_fastq_representative_fragments_parallel(
        reader,
        sampler,
        qc_filter,
        parallelism.threads(),
        parallelism.representative_batch_size(),
        FastqRetentionMode::Absolute,
    )
}

fn retain_fastq_proportional_representative_fragments(
    reader: &mut FastqPairReader,
    sampler: &RepresentativeSampler,
    qc_filter: &QcFilter,
) -> Result<FastqRepresentativeRetentionMetrics> {
    let mut qc_stats = QcStats::default();
    let mut retention = RepresentativeRetention::new(sampler.max_fragments());
    let mut qc_hash_wall_duration = Duration::default();

    let io_wall_duration = reader.for_each_borrowed_pair_with_timing(|fragment| {
        let qc_hash_started = Instant::now();
        retain_borrowed_fastq_proportional_representative_fragment(
            fragment,
            sampler,
            qc_filter,
            &mut qc_stats,
            &mut retention,
        );
        qc_hash_wall_duration = qc_hash_wall_duration.saturating_add(qc_hash_started.elapsed());
        Ok(())
    })?;

    Ok(FastqRepresentativeRetentionMetrics {
        qc_stats,
        retained_fragments: retention.into_sorted_vec(),
        io_wall_duration,
        qc_hash_wall_duration,
    })
}

fn retain_fastq_proportional_representative_fragments_with_parallelism(
    reader: &mut FastqPairReader,
    sampler: &RepresentativeSampler,
    qc_filter: &QcFilter,
    parallelism: &RoundParallelism,
) -> Result<FastqRepresentativeRetentionMetrics> {
    if parallelism.threads() <= 1 {
        return retain_fastq_proportional_representative_fragments(reader, sampler, qc_filter);
    }

    retain_fastq_representative_fragments_parallel(
        reader,
        sampler,
        qc_filter,
        parallelism.threads(),
        parallelism.representative_batch_size(),
        FastqRetentionMode::Proportional,
    )
}

fn retain_fastq_representative_fragments_parallel(
    reader: &mut FastqPairReader,
    sampler: &RepresentativeSampler,
    qc_filter: &QcFilter,
    worker_count: usize,
    batch_size: usize,
    mode: FastqRetentionMode,
) -> Result<FastqRepresentativeRetentionMetrics> {
    let worker_count = worker_count.max(1);
    let batch_size = batch_size.max(1);
    let mut qc_stats = QcStats::default();
    let mut retention = RepresentativeRetention::new(sampler.max_fragments());
    let mut qc_hash_wall_duration = Duration::default();
    let mut pending_results = BTreeMap::new();
    let mut next_result_sequence = 0usize;
    let mut submitted_batches = 0usize;
    let mut received_batches = 0usize;

    let io_wall_duration = std::thread::scope(|scope| -> Result<Duration> {
        let (work_tx, work_rx) = bounded(worker_count * FASTQ_RETENTION_WORK_QUEUE_MULTIPLIER);
        let (result_tx, result_rx) = unbounded();
        let mut worker_handles = Vec::with_capacity(worker_count);

        for _ in 0..worker_count {
            let work_rx = work_rx.clone();
            let result_tx = result_tx.clone();
            worker_handles.push(scope.spawn(move || {
                while let Ok(batch) = work_rx.recv() {
                    let result = score_fastq_retention_batch(batch, sampler, qc_filter, mode);
                    if result_tx.send(result).is_err() {
                        break;
                    }
                }
            }));
        }
        drop(result_tx);

        let mut batch_sequence = 0usize;
        let mut batch = Vec::with_capacity(batch_size);
        let read_result = reader.for_each_borrowed_pair_with_timing(|fragment| {
            batch.push(fragment.into_owned_without_qualities());
            if batch.len() == batch_size {
                submit_fastq_retention_batch(
                    &work_tx,
                    &mut batch_sequence,
                    &mut submitted_batches,
                    std::mem::take(&mut batch),
                )?;
                batch = Vec::with_capacity(batch_size);
                receive_available_fastq_retention_results(
                    &result_rx,
                    &mut pending_results,
                    &mut received_batches,
                )?;
                merge_ready_fastq_retention_results(
                    &mut pending_results,
                    &mut next_result_sequence,
                    &mut qc_stats,
                    &mut retention,
                    sampler,
                    qc_filter,
                    &mut qc_hash_wall_duration,
                )?;
            }
            Ok(())
        });

        let io_wall_duration = match read_result {
            Ok(duration) => duration,
            Err(error) => {
                drop(work_tx);
                join_fastq_retention_workers(worker_handles)?;
                return Err(error);
            }
        };

        if !batch.is_empty() {
            submit_fastq_retention_batch(
                &work_tx,
                &mut batch_sequence,
                &mut submitted_batches,
                batch,
            )?;
        }
        drop(work_tx);

        while received_batches < submitted_batches {
            let result = result_rx.recv().map_err(|_| {
                RvScreenError::validation(
                    "pipeline.parallel",
                    "FASTQ retention worker result channel closed before all batches were merged",
                )
            })?;
            insert_fastq_retention_result(&mut pending_results, &mut received_batches, result)?;
            merge_ready_fastq_retention_results(
                &mut pending_results,
                &mut next_result_sequence,
                &mut qc_stats,
                &mut retention,
                sampler,
                qc_filter,
                &mut qc_hash_wall_duration,
            )?;
        }

        join_fastq_retention_workers(worker_handles)?;
        Ok(io_wall_duration)
    })?;

    Ok(FastqRepresentativeRetentionMetrics {
        qc_stats,
        retained_fragments: retention.into_sorted_vec(),
        io_wall_duration,
        qc_hash_wall_duration,
    })
}

fn submit_fastq_retention_batch(
    work_tx: &crossbeam_channel::Sender<FastqRetentionBatch>,
    batch_sequence: &mut usize,
    submitted_batches: &mut usize,
    fragments: Vec<IoFragmentRecord>,
) -> Result<()> {
    work_tx
        .send(FastqRetentionBatch {
            sequence: *batch_sequence,
            fragments,
        })
        .map_err(|_| {
            RvScreenError::validation(
                "pipeline.parallel",
                "FASTQ retention worker queue closed before batch submission completed",
            )
        })?;
    *batch_sequence = batch_sequence.saturating_add(1);
    *submitted_batches = submitted_batches.saturating_add(1);
    Ok(())
}

fn score_fastq_retention_batch(
    batch: FastqRetentionBatch,
    sampler: &RepresentativeSampler,
    qc_filter: &QcFilter,
    mode: FastqRetentionMode,
) -> FastqRetentionBatchResult {
    let qc_hash_started = Instant::now();
    let entries = batch
        .fragments
        .into_iter()
        .map(|fragment| score_fastq_retention_fragment(fragment, sampler, qc_filter, mode))
        .collect();

    FastqRetentionBatchResult {
        sequence: batch.sequence,
        entries,
        qc_hash_wall_duration: qc_hash_started.elapsed(),
    }
}

fn score_fastq_retention_fragment(
    fragment: IoFragmentRecord,
    sampler: &RepresentativeSampler,
    qc_filter: &QcFilter,
    mode: FastqRetentionMode,
) -> FastqRetentionBatchEntry {
    let qc_fragment = bridge_qc_fragment(&fragment);
    let cheap_pass = match qc_filter.filter_cheap(&qc_fragment) {
        CheapQcResult::Pass(pass) => pass,
        CheapQcResult::Filtered(failure) => {
            return FastqRetentionBatchEntry::CheapFiltered(failure.reason);
        }
    };

    match mode {
        FastqRetentionMode::Absolute => {
            let rank_hash = sampler.hash_fragment_key(&fragment.fragment_key);
            FastqRetentionBatchEntry::CheapPassed {
                fragment,
                cheap_pass,
                rank_hash,
            }
        }
        FastqRetentionMode::Proportional => {
            match qc_filter.complete_entropy_validation(&qc_fragment, cheap_pass) {
                crate::qc::QcResult::Pass(_) => {
                    let rank_hash = sampler.hash_fragment_key(&fragment.fragment_key);
                    FastqRetentionBatchEntry::Passed {
                        fragment,
                        rank_hash,
                    }
                }
                crate::qc::QcResult::Filtered(failure) => {
                    FastqRetentionBatchEntry::Filtered(failure.reason)
                }
            }
        }
    }
}

fn receive_available_fastq_retention_results(
    result_rx: &Receiver<FastqRetentionBatchResult>,
    pending_results: &mut BTreeMap<usize, FastqRetentionBatchResult>,
    received_batches: &mut usize,
) -> Result<()> {
    while let Ok(result) = result_rx.try_recv() {
        insert_fastq_retention_result(pending_results, received_batches, result)?;
    }

    Ok(())
}

fn insert_fastq_retention_result(
    pending_results: &mut BTreeMap<usize, FastqRetentionBatchResult>,
    received_batches: &mut usize,
    result: FastqRetentionBatchResult,
) -> Result<()> {
    if pending_results.insert(result.sequence, result).is_some() {
        return Err(RvScreenError::validation(
            "pipeline.parallel",
            "received duplicate FASTQ retention batch result sequence",
        ));
    }
    *received_batches = received_batches.saturating_add(1);
    Ok(())
}

fn merge_ready_fastq_retention_results(
    pending_results: &mut BTreeMap<usize, FastqRetentionBatchResult>,
    next_result_sequence: &mut usize,
    qc_stats: &mut QcStats,
    retention: &mut RepresentativeRetention<IoFragmentRecord>,
    sampler: &RepresentativeSampler,
    qc_filter: &QcFilter,
    qc_hash_wall_duration: &mut Duration,
) -> Result<()> {
    while let Some(result) = pending_results.remove(next_result_sequence) {
        merge_fastq_retention_result(
            result,
            qc_stats,
            retention,
            sampler,
            qc_filter,
            qc_hash_wall_duration,
        );
        *next_result_sequence = next_result_sequence.saturating_add(1);
    }

    Ok(())
}

fn merge_fastq_retention_result(
    result: FastqRetentionBatchResult,
    qc_stats: &mut QcStats,
    retention: &mut RepresentativeRetention<IoFragmentRecord>,
    sampler: &RepresentativeSampler,
    qc_filter: &QcFilter,
    qc_hash_wall_duration: &mut Duration,
) {
    *qc_hash_wall_duration = qc_hash_wall_duration.saturating_add(result.qc_hash_wall_duration);
    for entry in result.entries {
        qc_stats.record_total_fragment();
        match entry {
            FastqRetentionBatchEntry::CheapFiltered(reason)
            | FastqRetentionBatchEntry::Filtered(reason) => {
                qc_stats.record_failure_reason(&reason);
            }
            FastqRetentionBatchEntry::CheapPassed {
                fragment,
                cheap_pass,
                rank_hash,
            } => {
                merge_absolute_fastq_retention_entry(
                    fragment,
                    cheap_pass,
                    rank_hash,
                    qc_stats,
                    retention,
                    sampler,
                    qc_filter,
                    qc_hash_wall_duration,
                );
            }
            FastqRetentionBatchEntry::Passed {
                fragment,
                rank_hash,
            } => {
                qc_stats.record_passed_fragment();
                let fragment_key = fragment.fragment_key.as_str();
                if retention.would_retain_hash(rank_hash, fragment_key) {
                    let rank = sampler.rank_fragment(fragment_key);
                    retention.consider(RankedRepresentative::new(rank, fragment));
                }
            }
        }
    }
}

fn merge_absolute_fastq_retention_entry(
    fragment: IoFragmentRecord,
    cheap_pass: CheapQcPass,
    rank_hash: u64,
    qc_stats: &mut QcStats,
    retention: &mut RepresentativeRetention<IoFragmentRecord>,
    sampler: &RepresentativeSampler,
    qc_filter: &QcFilter,
    qc_hash_wall_duration: &mut Duration,
) {
    let fragment_key = fragment.fragment_key.as_str();
    if !retention.would_retain_hash(rank_hash, fragment_key) {
        qc_stats.record_passed_fragment();
        return;
    }

    let entropy_started = Instant::now();
    let qc_fragment = bridge_qc_fragment(&fragment);
    let qc_result = qc_filter.complete_entropy_validation(&qc_fragment, cheap_pass);
    *qc_hash_wall_duration = qc_hash_wall_duration.saturating_add(entropy_started.elapsed());

    match qc_result {
        crate::qc::QcResult::Pass(_) => {
            qc_stats.record_passed_fragment();
            let rank = sampler.rank_fragment(fragment_key);
            retention.consider(RankedRepresentative::new(rank, fragment));
        }
        crate::qc::QcResult::Filtered(failure) => {
            qc_stats.record_failure_reason(&failure.reason);
        }
    }
}

fn join_fastq_retention_workers(
    worker_handles: Vec<std::thread::ScopedJoinHandle<'_, ()>>,
) -> Result<()> {
    for handle in worker_handles {
        if handle.join().is_err() {
            return Err(RvScreenError::validation(
                "pipeline.parallel",
                "FASTQ retention worker panicked",
            ));
        }
    }
    Ok(())
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

pub fn resolve_sampling_round_plan(
    profile_sampling: &SamplingConfig,
    overrides: SamplingRoundOverrides<'_>,
) -> Result<ResolvedSamplingRoundPlan> {
    if overrides.rounds.is_some() && overrides.round_proportions.is_some() {
        return Err(RvScreenError::validation(
            "sampling overrides",
            "--rounds and --round-proportions are mutually exclusive sampling overrides",
        ));
    }

    if let Some(rounds) = overrides.rounds {
        return active_round_targets(rounds, rounds.len() as u64)
            .map(ResolvedSamplingRoundPlan::Absolute);
    }

    if let Some(proportions) = overrides.round_proportions {
        return active_round_proportions(proportions, proportions.len() as u64)
            .map(ResolvedSamplingRoundPlan::Proportional);
    }

    match profile_sampling.round_mode {
        Some(SamplingRoundMode::Absolute) => {
            active_round_targets(&profile_sampling.rounds, profile_sampling.max_rounds)
                .map(ResolvedSamplingRoundPlan::Absolute)
        }
        Some(SamplingRoundMode::Proportional) => {
            let proportions = profile_sampling
                .round_proportions
                .as_deref()
                .ok_or_else(|| {
                    RvScreenError::validation(
                        "sampling.round_proportions",
                        "sampling.round_mode = \"proportional\" requires sampling.round_proportions",
                    )
                })?;
            active_round_proportions(proportions, profile_sampling.max_rounds)
                .map(ResolvedSamplingRoundPlan::Proportional)
        }
        None => {
            if let Some(proportions) = profile_sampling
                .round_proportions
                .as_deref()
                .filter(|proportions| !proportions.is_empty())
            {
                active_round_proportions(proportions, profile_sampling.max_rounds)
                    .map(ResolvedSamplingRoundPlan::Proportional)
            } else {
                active_round_targets(&profile_sampling.rounds, profile_sampling.max_rounds)
                    .map(ResolvedSamplingRoundPlan::Absolute)
            }
        }
    }
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

fn active_round_proportions(proportions: &[f64], max_rounds: u64) -> Result<Vec<f64>> {
    if proportions.is_empty() {
        return Err(RvScreenError::validation(
            "sampling.round_proportions",
            "screen requires at least one proportional sampling round",
        ));
    }

    if max_rounds == 0 {
        return Err(RvScreenError::validation(
            "sampling.max_rounds",
            "sampling.max_rounds must be greater than zero",
        ));
    }

    crate::types::validate_round_proportion_values(proportions)
        .map_err(|error| RvScreenError::validation("sampling.round_proportions", error))?;

    let active_len = usize::try_from(max_rounds)
        .unwrap_or(usize::MAX)
        .min(proportions.len());
    let active = proportions[..active_len].to_vec();

    if active.is_empty() {
        return Err(RvScreenError::validation(
            "sampling.round_proportions",
            "sampling.max_rounds truncated all configured proportional rounds",
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

fn bridge_fastq_qc_fragment<'a>(fragment: &'a BorrowedFragmentRecord<'a>) -> QcFragmentRecord<'a> {
    QcFragmentRecord::paired(fragment.r1_seq, fragment.r2_seq)
}

fn retain_owned_representative_fragment(
    fragment: &IoFragmentRecord,
    sampler: &RepresentativeSampler,
    qc_filter: &QcFilter,
    qc_stats: &mut QcStats,
    retention: &mut RepresentativeRetention<IoFragmentRecord>,
) -> bool {
    qc_stats.record_total_fragment();

    let qc_fragment = bridge_qc_fragment(fragment);
    let cheap_pass = match qc_filter.filter_cheap(&qc_fragment) {
        CheapQcResult::Pass(pass) => pass,
        CheapQcResult::Filtered(failure) => {
            qc_stats.record_failure_reason(&failure.reason);
            return false;
        }
    };

    let fragment_key = fragment.fragment_key.as_str();
    let rank_hash = sampler.hash_fragment_key(fragment_key);
    if !retention.would_retain_hash(rank_hash, fragment_key) {
        qc_stats.record_passed_fragment();
        return false;
    }

    match qc_filter.complete_entropy_validation(&qc_fragment, cheap_pass) {
        crate::qc::QcResult::Pass(_) => {
            qc_stats.record_passed_fragment();
            let rank = sampler.rank_fragment(fragment_key);
            retention.consider(RankedRepresentative::new(rank, fragment.clone()));
            true
        }
        crate::qc::QcResult::Filtered(failure) => {
            qc_stats.record_failure_reason(&failure.reason);
            false
        }
    }
}

fn retain_owned_proportional_representative_fragment(
    fragment: &IoFragmentRecord,
    sampler: &RepresentativeSampler,
    qc_filter: &QcFilter,
    qc_stats: &mut QcStats,
    retention: &mut RepresentativeRetention<IoFragmentRecord>,
) {
    qc_stats.record_total_fragment();

    let qc_fragment = bridge_qc_fragment(fragment);
    let cheap_pass = match qc_filter.filter_cheap(&qc_fragment) {
        CheapQcResult::Pass(pass) => pass,
        CheapQcResult::Filtered(failure) => {
            qc_stats.record_failure_reason(&failure.reason);
            return;
        }
    };

    match qc_filter.complete_entropy_validation(&qc_fragment, cheap_pass) {
        crate::qc::QcResult::Pass(_) => {
            qc_stats.record_passed_fragment();
            let fragment_key = fragment.fragment_key.as_str();
            let rank_hash = sampler.hash_fragment_key(fragment_key);
            if retention.would_retain_hash(rank_hash, fragment_key) {
                let rank = sampler.rank_fragment(fragment_key);
                retention.consider(RankedRepresentative::new(rank, fragment.clone()));
            }
        }
        crate::qc::QcResult::Filtered(failure) => {
            qc_stats.record_failure_reason(&failure.reason);
        }
    }
}

fn retain_borrowed_fastq_representative_fragment(
    fragment: BorrowedFragmentRecord<'_>,
    sampler: &RepresentativeSampler,
    qc_filter: &QcFilter,
    qc_stats: &mut QcStats,
    retention: &mut RepresentativeRetention<IoFragmentRecord>,
) {
    qc_stats.record_total_fragment();

    let qc_fragment = bridge_fastq_qc_fragment(&fragment);
    let cheap_pass = match qc_filter.filter_cheap(&qc_fragment) {
        CheapQcResult::Pass(pass) => pass,
        CheapQcResult::Filtered(failure) => {
            qc_stats.record_failure_reason(&failure.reason);
            return;
        }
    };

    let fragment_key = fragment.fragment_key.as_ref();
    let rank_hash = sampler.hash_fragment_key(fragment_key);
    if !retention.would_retain_hash(rank_hash, fragment_key) {
        qc_stats.record_passed_fragment();
        return;
    }

    match qc_filter.complete_entropy_validation(&qc_fragment, cheap_pass) {
        crate::qc::QcResult::Pass(_) => {
            qc_stats.record_passed_fragment();
            let rank = sampler.rank_fragment(fragment_key);
            retention.consider(RankedRepresentative::new(
                rank,
                fragment.into_owned_without_qualities(),
            ));
        }
        crate::qc::QcResult::Filtered(failure) => {
            qc_stats.record_failure_reason(&failure.reason);
        }
    }
}

fn retain_borrowed_fastq_proportional_representative_fragment(
    fragment: BorrowedFragmentRecord<'_>,
    sampler: &RepresentativeSampler,
    qc_filter: &QcFilter,
    qc_stats: &mut QcStats,
    retention: &mut RepresentativeRetention<IoFragmentRecord>,
) {
    qc_stats.record_total_fragment();

    let qc_fragment = bridge_fastq_qc_fragment(&fragment);
    let cheap_pass = match qc_filter.filter_cheap(&qc_fragment) {
        CheapQcResult::Pass(pass) => pass,
        CheapQcResult::Filtered(failure) => {
            qc_stats.record_failure_reason(&failure.reason);
            return;
        }
    };

    match qc_filter.complete_entropy_validation(&qc_fragment, cheap_pass) {
        crate::qc::QcResult::Pass(_) => {
            qc_stats.record_passed_fragment();
            let fragment_key = fragment.fragment_key.as_ref();
            let rank_hash = sampler.hash_fragment_key(fragment_key);
            if retention.would_retain_hash(rank_hash, fragment_key) {
                let rank = sampler.rank_fragment(fragment_key);
                retention.consider(RankedRepresentative::new(
                    rank,
                    fragment.into_owned_without_qualities(),
                ));
            }
        }
        crate::qc::QcResult::Filtered(failure) => {
            qc_stats.record_failure_reason(&failure.reason);
        }
    }
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

fn absolute_sampling_round_plan(round_targets: &[u64]) -> SamplingRoundPlanReport {
    SamplingRoundPlanReport {
        round_mode: SamplingRoundMode::Absolute,
        requested_rounds: Some(round_targets.to_vec()),
        requested_proportions: None,
        raw_proportional_counts: None,
        effective_counts: round_targets.to_vec(),
        qc_passing_denominator: None,
        warnings: Vec::new(),
    }
}

fn proportional_sampling_round_plan(
    round_proportions: &[f64],
    qc_passing_denominator: u64,
    resolution: &ProportionalRoundResolution,
) -> SamplingRoundPlanReport {
    SamplingRoundPlanReport {
        round_mode: SamplingRoundMode::Proportional,
        requested_rounds: None,
        requested_proportions: Some(round_proportions.to_vec()),
        raw_proportional_counts: Some(resolution.raw_counts.clone()),
        effective_counts: resolution.effective_counts.clone(),
        qc_passing_denominator: Some(qc_passing_denominator),
        warnings: resolution
            .warnings
            .iter()
            .map(sampling_warning_report)
            .collect(),
    }
}

fn sampling_warning_report(warning: &RepresentativeSamplingWarning) -> SamplingWarningReport {
    SamplingWarningReport {
        kind: warning.kind.as_str().to_string(),
        round: warning.round,
        message: warning.message.clone(),
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
            rounds: None,
            round_proportions: None,
            allow_sampling_threshold_override: false,
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
            rounds: None,
            round_proportions: None,
            allow_sampling_threshold_override: false,
            threads: 4,
        })
        .expect("four-thread screen should succeed");

        assert_screen_outcomes_equal(&outcome_single, &outcome_parallel);
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
                rounds: None,
                round_proportions: None,
                allow_sampling_threshold_override: false,
                threads: 4,
            })
            .expect("repeated multi-thread screen should succeed");

            match &baseline {
                Some(expected) => assert_screen_outcomes_equal(expected, &outcome),
                None => baseline = Some(outcome),
            }
        }
    }

    #[test]
    fn representative_results_match_across_thread_counts_when_all_rounds_run() {
        let tempdir = tempdir().expect("tempdir should be created");
        let bundle =
            prepare_reference_bundle(tempdir.path()).expect("reference bundle should be prepared");
        let calibration_dir = write_calibration_profile(
            tempdir.path().join("calibration-two-round"),
            &bundle.version,
            false,
            1.0,
        )
        .expect("calibration profile should be written");
        let (r1, r2) = generate_fastq_pair(
            &FastqPairConfig::new(240, 100, 6127)
                .with_output_dir(tempdir.path().join("fastq-two-round"))
                .with_components(vec![
                    ReadComponent::new(SyntheticSource::Human, 0.98),
                    ReadComponent::new(
                        SyntheticSource::Virus(VirusSelector::Accession(
                            "NC_SYNTHV1.1".to_string(),
                        )),
                        0.02,
                    ),
                ]),
        )
        .expect("FASTQ should be generated");

        let outcome_single = run_screen(&ScreenArgs {
            input: vec![r1.clone(), r2.clone()],
            reference_bundle: bundle.bundle_dir.clone(),
            calibration_profile: calibration_dir.clone(),
            negative_control: None,
            out: tempdir.path().join("report-two-round-single-thread"),
            mode: ScreenMode::Representative,
            rounds: None,
            round_proportions: None,
            allow_sampling_threshold_override: false,
            threads: 1,
        })
        .expect("single-thread two-round representative screen should succeed");

        let outcome_parallel = run_screen(&ScreenArgs {
            input: vec![r1, r2],
            reference_bundle: bundle.bundle_dir,
            calibration_profile: calibration_dir,
            negative_control: None,
            out: tempdir.path().join("report-two-round-four-thread"),
            mode: ScreenMode::Representative,
            rounds: None,
            round_proportions: None,
            allow_sampling_threshold_override: false,
            threads: 4,
        })
        .expect("four-thread two-round representative screen should succeed");

        assert_eq!(outcome_single.rounds.len(), 2);
        assert_eq!(outcome_parallel.rounds.len(), 2);
        assert_eq!(outcome_single.rounds[0].sampled_fragments, 50);
        assert_eq!(outcome_single.rounds[1].sampled_fragments, 100);
        assert_screen_outcomes_equal(&outcome_single, &outcome_parallel);
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
            rounds: None,
            round_proportions: None,
            allow_sampling_threshold_override: false,
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

        let round_targets = match &prepared.round_plan {
            ResolvedSamplingRoundPlan::Absolute(round_targets) => round_targets.as_slice(),
            ResolvedSamplingRoundPlan::Proportional(_) => {
                panic!("single-pass representative test expects absolute rounds")
            }
        };
        let outcome = execute_representative_rounds(
            &input,
            &prepared.loaded_profile.profile.seed,
            round_targets,
            &dependencies,
        )
        .expect("single-pass representative execution should succeed");

        assert_eq!(input.open_count(), 1);
        assert!(!outcome.rounds.is_empty());
    }

    #[test]
    fn representative_fastq_retention_clones_sequences_without_decoding_qualities() {
        let tempdir = tempdir().expect("tempdir should be created");
        let r1 = tempdir.path().join("reads_R1.fastq");
        let r2 = tempdir.path().join("reads_R2.fastq");
        let seq1 = "ACGT".repeat(25);
        let seq2 = "TGCA".repeat(25);
        let invalid_qual = String::from_utf8(vec![31; 100]).expect("quality bytes should build");
        fs::write(&r1, format!("@frag/1\n{seq1}\n+\n{invalid_qual}\n"))
            .expect("R1 FASTQ should be written");
        fs::write(&r2, format!("@frag/2\n{seq2}\n+\n{invalid_qual}\n"))
            .expect("R2 FASTQ should be written");

        let mut reader = FastqPairReader::open(&r1, &r2).expect("FASTQ reader should open");
        let sampler =
            RepresentativeSampler::new(7, vec![1]).expect("representative sampler should build");
        let qc_filter = QcFilter::default();

        let FastqRepresentativeRetentionMetrics {
            qc_stats,
            retained_fragments: retained,
            ..
        } = retain_fastq_representative_fragments(&mut reader, &sampler, &qc_filter)
            .expect("representative FASTQ retention should not decode qualities");

        assert_eq!(qc_stats.total_fragments, 1);
        assert_eq!(qc_stats.passed_fragments, 1);
        assert_eq!(retained.len(), 1);

        let fragment = retained
            .into_iter()
            .next()
            .expect("one fragment should be retained")
            .into_payload();
        assert_eq!(fragment.fragment_key, "frag");
        assert_eq!(fragment.r1_seq, seq1.as_bytes());
        assert_eq!(fragment.r2_seq, seq2.as_bytes());
        assert!(fragment.r1_qual.is_empty());
        assert!(fragment.r2_qual.is_empty());
    }

    #[test]
    fn representative_fastq_retention_filters_low_complexity_after_cheap_admission() {
        let tempdir = tempdir().expect("tempdir should be created");
        let r1 = tempdir.path().join("reads_R1.fastq");
        let r2 = tempdir.path().join("reads_R2.fastq");
        let read1 = format!("{}C", "A".repeat(99));
        let read2 = "ACGT".repeat(25);
        let qual = "I".repeat(100);
        fs::write(&r1, format!("@frag/1\n{read1}\n+\n{qual}\n"))
            .expect("R1 FASTQ should be written");
        fs::write(&r2, format!("@frag/2\n{read2}\n+\n{qual}\n"))
            .expect("R2 FASTQ should be written");

        let mut reader = FastqPairReader::open(&r1, &r2).expect("FASTQ reader should open");
        let sampler =
            RepresentativeSampler::new(7, vec![1]).expect("representative sampler should build");
        let qc_filter = QcFilter::default();

        let FastqRepresentativeRetentionMetrics {
            qc_stats,
            retained_fragments: retained,
            ..
        } = retain_fastq_representative_fragments(&mut reader, &sampler, &qc_filter)
            .expect("representative FASTQ retention should apply full entropy before retaining");

        assert_eq!(qc_stats.total_fragments, 1);
        assert_eq!(qc_stats.passed_fragments, 0);
        assert_eq!(qc_stats.low_complexity_filtered, 1);
        assert!(retained.is_empty());
    }

    #[test]
    fn representative_fastq_retention_reports_nonzero_io_and_qc_hash_timing() {
        let tempdir = tempdir().expect("tempdir should be created");
        let r1 = tempdir.path().join("reads_R1.fastq");
        let r2 = tempdir.path().join("reads_R2.fastq");
        let read1 = "ACGT".repeat(25);
        let read2 = "TGCA".repeat(25);
        let qual = "I".repeat(100);
        fs::write(&r1, format!("@frag/1\n{read1}\n+\n{qual}\n"))
            .expect("R1 FASTQ should be written");
        fs::write(&r2, format!("@frag/2\n{read2}\n+\n{qual}\n"))
            .expect("R2 FASTQ should be written");

        let mut reader = FastqPairReader::open(&r1, &r2).expect("FASTQ reader should open");
        let sampler =
            RepresentativeSampler::new(7, vec![1]).expect("representative sampler should build");
        let qc_filter = QcFilter::default();

        let metrics = retain_fastq_representative_fragments(&mut reader, &sampler, &qc_filter)
            .expect("representative FASTQ retention should collect timing metrics");

        assert!(
            metrics.io_wall_duration > Duration::default(),
            "borrowed FASTQ representative path should report nonzero I/O timing"
        );
        assert!(
            metrics.qc_hash_wall_duration > Duration::default(),
            "borrowed FASTQ representative path should report nonzero QC/hash timing"
        );
    }

    #[test]
    fn representative_fastq_parallel_retention_preserves_ids_and_round_membership() {
        let tempdir = tempdir().expect("tempdir should be created");
        let (r1, r2) = generate_fastq_pair(
            &FastqPairConfig::new(180, 100, 9721)
                .with_output_dir(tempdir.path().join("fastq-retention-threadcount"))
                .with_components(vec![ReadComponent::new(SyntheticSource::Human, 1.0)]),
        )
        .expect("FASTQ should be generated");
        let sampler = RepresentativeSampler::new(7, vec![25, 80])
            .expect("representative sampler should build");
        let qc_filter = QcFilter::default();

        let mut sequential_reader = FastqPairReader::open(&r1, &r2).expect("FASTQ reader opens");
        let sequential =
            retain_fastq_representative_fragments(&mut sequential_reader, &sampler, &qc_filter)
                .expect("sequential retention should succeed");
        let sequential_qc_stats = sequential.qc_stats;
        let sequential_memberships = retained_round_memberships(sequential, &sampler);

        let parallelism = RoundParallelism::new(16).expect("parallelism should build");
        let mut parallel_reader = FastqPairReader::open(&r1, &r2).expect("FASTQ reader opens");
        let parallel = retain_fastq_representative_fragments_with_parallelism(
            &mut parallel_reader,
            &sampler,
            &qc_filter,
            &parallelism,
        )
        .expect("parallel retention should succeed");
        let parallel_qc_stats = parallel.qc_stats;
        let parallel_memberships = retained_round_memberships(parallel, &sampler);

        assert_eq!(sequential_qc_stats, parallel_qc_stats);
        assert_eq!(sequential_memberships, parallel_memberships);
        let mut membership_counts = vec![0usize; sampler.round_count()];
        for (entry_round, _) in &sequential_memberships {
            membership_counts[*entry_round] += 1;
        }
        println!(
            "retained_id_threadcount_equality threads_1_vs_16=true retained={} round_membership_counts={:?} first_memberships={:?}",
            sequential_memberships.len(),
            membership_counts,
            sequential_memberships.iter().take(8).collect::<Vec<_>>()
        );
    }

    #[test]
    fn run_screen_populates_runtime_peak_heap_metric() {
        let tempdir = tempdir().expect("tempdir should be created");
        let bundle =
            prepare_reference_bundle(tempdir.path()).expect("reference bundle should be prepared");
        let calibration_dir = write_calibration_profile(
            tempdir.path().join("calibration-runtime-metrics"),
            &bundle.version,
            false,
            0.0,
        )
        .expect("calibration profile should be written");
        let (r1, r2) = generate_fastq_pair(
            &FastqPairConfig::new(32, 100, 5511)
                .with_output_dir(tempdir.path().join("fastq-runtime-metrics"))
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

        let outcome = run_screen(&ScreenArgs {
            input: vec![r1, r2],
            reference_bundle: bundle.bundle_dir,
            calibration_profile: calibration_dir,
            negative_control: None,
            out: tempdir.path().join("report-runtime-metrics"),
            mode: ScreenMode::Representative,
            rounds: None,
            round_proportions: None,
            allow_sampling_threshold_override: false,
            threads: 1,
        })
        .expect("representative screen should succeed");

        assert!(
            outcome.perf_metrics.peak_heap_bytes > 0,
            "runtime perf metrics should populate peak_heap_bytes before bench serialization"
        );
        assert!(
            outcome.perf_metrics.io_wall_duration > Duration::default(),
            "runtime perf metrics should carry representative FASTQ I/O timing"
        );
        let labels = outcome
            .perf_metrics
            .stage_timings
            .iter()
            .map(|timing| timing.label.as_str())
            .collect::<Vec<_>>();
        for label in REPRESENTATIVE_STARTUP_STAGE_LABELS {
            assert!(
                labels.contains(&label),
                "runtime perf metrics should include stage timing label {label}"
            );
        }
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
            rounds: None,
            round_proportions: None,
            allow_sampling_threshold_override: false,
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
            rounds: None,
            round_proportions: None,
            allow_sampling_threshold_override: false,
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
        assert_screen_outcomes_equal(&first, &second);
    }

    fn retained_round_memberships(
        metrics: FastqRepresentativeRetentionMetrics,
        sampler: &RepresentativeSampler,
    ) -> Vec<(usize, String)> {
        let prefix_sample_sizes = sampler.prefix_sample_sizes(metrics.retained_fragments.len());
        metrics
            .retained_fragments
            .into_iter()
            .enumerate()
            .map(|(prefix_index, retained)| {
                let entry_round = RepresentativeRoundSpool::entry_round_for_prefix_index(
                    &prefix_sample_sizes,
                    prefix_index as u64,
                )
                .expect("retained prefix should map to a representative round");
                (entry_round, retained.into_payload().fragment_key)
            })
            .collect()
    }

    fn assert_screen_outcomes_equal(left: &ScreenRunOutcome, right: &ScreenRunOutcome) {
        let mut left = left.clone();
        let mut right = right.clone();
        left.perf_metrics = ScreenPerfMetrics::default();
        right.perf_metrics = ScreenPerfMetrics::default();

        assert_eq!(left, right);
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

    #[test]
    fn sampling_resolution_preserves_legacy_absolute_profile_rounds() {
        let plan = resolve_sampling_round_plan(
            &sampling_config(vec![50, 100], None, None, 2),
            SamplingRoundOverrides::default(),
        )
        .expect("legacy absolute sampling should resolve");

        assert_eq!(plan, ResolvedSamplingRoundPlan::Absolute(vec![50, 100]));
        assert_eq!(plan.mode(), SamplingRoundMode::Absolute);
    }

    #[test]
    fn sampling_resolution_prefers_non_empty_proportions_without_round_mode() {
        let plan = resolve_sampling_round_plan(
            &sampling_config(vec![50, 100], None, Some(vec![0.25, 0.5, 1.0]), 2),
            SamplingRoundOverrides::default(),
        )
        .expect("implicit proportional sampling should resolve");

        assert_eq!(
            plan,
            ResolvedSamplingRoundPlan::Proportional(vec![0.25, 0.5])
        );
        assert_eq!(plan.mode(), SamplingRoundMode::Proportional);
    }

    #[test]
    fn sampling_resolution_honors_explicit_round_modes() {
        let profile = sampling_config(
            vec![50, 100],
            Some(SamplingRoundMode::Absolute),
            Some(vec![0.25, 0.5]),
            2,
        );
        let absolute = resolve_sampling_round_plan(&profile, SamplingRoundOverrides::default())
            .expect("explicit absolute mode should use rounds");
        assert_eq!(absolute, ResolvedSamplingRoundPlan::Absolute(vec![50, 100]));

        let profile = sampling_config(
            vec![50, 100],
            Some(SamplingRoundMode::Proportional),
            Some(vec![0.25, 0.5]),
            2,
        );
        let proportional = resolve_sampling_round_plan(&profile, SamplingRoundOverrides::default())
            .expect("explicit proportional mode should use proportions");
        assert_eq!(
            proportional,
            ResolvedSamplingRoundPlan::Proportional(vec![0.25, 0.5])
        );
    }

    #[test]
    fn sampling_resolution_prefers_cli_over_profile() {
        let profile = sampling_config(
            vec![50, 100],
            Some(SamplingRoundMode::Proportional),
            Some(vec![0.25, 0.5]),
            2,
        );
        let cli_rounds = [10, 20, 40];
        let plan = resolve_sampling_round_plan(
            &profile,
            SamplingRoundOverrides {
                rounds: Some(&cli_rounds),
                round_proportions: None,
            },
        )
        .expect("CLI absolute rounds should override profile proportional rounds");

        assert_eq!(plan, ResolvedSamplingRoundPlan::Absolute(vec![10, 20, 40]));
    }

    #[test]
    fn task8_sampling_resolution_covers_cli_absolute_and_proportional_overrides() {
        let profile = sampling_config(
            vec![50, 100],
            Some(SamplingRoundMode::Absolute),
            Some(vec![0.25, 0.5]),
            2,
        );
        let cli_rounds = [10, 20, 40];
        let absolute = resolve_sampling_round_plan(
            &profile,
            SamplingRoundOverrides {
                rounds: Some(&cli_rounds),
                round_proportions: None,
            },
        )
        .expect("CLI absolute rounds should override profile sampling");
        assert_eq!(
            absolute,
            ResolvedSamplingRoundPlan::Absolute(vec![10, 20, 40])
        );
        assert_eq!(absolute.mode(), SamplingRoundMode::Absolute);

        let cli_proportions = [0.25, 0.5, 1.0];
        let proportional = resolve_sampling_round_plan(
            &profile,
            SamplingRoundOverrides {
                rounds: None,
                round_proportions: Some(&cli_proportions),
            },
        )
        .expect("CLI proportional rounds should override profile sampling");
        assert_eq!(
            proportional,
            ResolvedSamplingRoundPlan::Proportional(vec![0.25, 0.5, 1.0])
        );
        assert_eq!(proportional.mode(), SamplingRoundMode::Proportional);
    }

    #[test]
    fn sampling_resolution_rejects_mixed_cli_overrides() {
        let profile = sampling_config(vec![50, 100], None, None, 2);
        let cli_rounds = [10, 20];
        let cli_proportions = [0.5, 1.0];
        let error = resolve_sampling_round_plan(
            &profile,
            SamplingRoundOverrides {
                rounds: Some(&cli_rounds),
                round_proportions: Some(&cli_proportions),
            },
        )
        .expect_err("mixed CLI overrides should be rejected");

        let error = error.to_string();
        assert!(error.contains("--rounds"), "{error}");
        assert!(error.contains("--round-proportions"), "{error}");
        assert!(error.contains("mutually exclusive"), "{error}");
    }

    fn sampling_config(
        rounds: Vec<u64>,
        round_mode: Option<SamplingRoundMode>,
        round_proportions: Option<Vec<f64>>,
        max_rounds: u64,
    ) -> SamplingConfig {
        SamplingConfig {
            mode: "representative".to_string(),
            rounds,
            round_mode,
            round_proportions,
            max_rounds,
        }
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
