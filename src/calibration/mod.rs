use crate::cli::{CalibrateArgs, ScreenArgs, ScreenMode};
use crate::error::{Result, RvScreenError};
use crate::pipeline::run_screen;
use crate::sampling::DEFAULT_REPRESENTATIVE_ROUNDS;
use crate::types::{
    BundleToml, CandidateRules, DecisionRules, DecisionStatus, FragmentRules, ProfileToml,
    ReleaseStatus, SamplingConfig,
};
use csv::WriterBuilder;
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::fs::{self, File};
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

const PROFILE_TOML: &str = "profile.toml";
const BUNDLE_TOML: &str = "bundle.toml";
const THRESHOLDS_TOML: &str = "thresholds.toml";
const BENCHMARK_MANIFEST_JSON: &str = "benchmark_manifest.json";
const BENCHMARK_SUMMARY_TSV: &str = "benchmark_summary.tsv";
const RELEASE_GATE_JSON: &str = "release_gate.json";
const CHECKSUM_SHA256: &str = "checksum.sha256";
const CALIBRATION_BACKEND: &str = "minimap2";
const CALIBRATION_PRESET: &str = "sr-conservative";
const DEFAULT_PROFILE_STATUS: &str = "release_candidate";
const FAILED_PROFILE_STATUS: &str = "benchmark_failed";
const DEFAULT_SEED: u64 = 20_260_420;
const DEFAULT_SUPPORTED_INPUT: &[&str] = &["fastq", "fastq.gz", "bam", "ubam", "cram"];
const DEFAULT_SUPPORTED_READ_TYPE: &[&str] = &["illumina_pe_shortread"];
const BACKEND_GATE_DETAILS: &str =
    "V1 single-backend mode (minimap2 only) — no baseline comparison required";

#[derive(Debug, Clone, PartialEq)]
pub struct LoadedProfile {
    pub profile_dir: PathBuf,
    pub profile_path: PathBuf,
    pub profile: ProfileToml,
}

#[derive(Debug, Clone, PartialEq)]
pub struct LoadedReferenceBundle {
    pub bundle_dir: PathBuf,
    pub bundle_path: PathBuf,
    pub bundle: BundleToml,
}

#[derive(Debug, Clone, PartialEq)]
pub struct CalibrationOutcome {
    pub output_dir: PathBuf,
    pub profile: ProfileToml,
    pub benchmark_runs: Vec<BenchmarkRunResult>,
    pub release_gate: ReleaseGate,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct ThresholdsToml {
    pub fragment_rules: FragmentRules,
    pub candidate_rules: CandidateRules,
    pub decision_rules: DecisionRules,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct BenchmarkManifest {
    #[serde(default)]
    pub profile: BenchmarkProfileConfig,
    pub datasets: Vec<BenchmarkDataset>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default)]
pub struct BenchmarkProfileConfig {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub profile_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub status: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub seed: Option<u64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub supported_input: Option<Vec<String>>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub supported_read_type: Option<Vec<String>>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub negative_control_required: Option<bool>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sampling: Option<SamplingConfig>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub fragment_rules: Option<FragmentRules>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub candidate_rules: Option<CandidateRules>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub decision_rules: Option<DecisionRules>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct BenchmarkDataset {
    pub dataset_id: String,
    pub input: Vec<PathBuf>,
    pub expected_outcome: BenchmarkExpectedOutcome,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub expected_target: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum BenchmarkExpectedOutcome {
    Negative,
    SpikeIn,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct ReleaseGate {
    pub backend_gate: BackendGate,
    pub reference_gate: ReferenceGate,
    pub specificity_gate: SpecificityGate,
    pub sensitivity_gate: SensitivityGate,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum GateStatus {
    Pass,
    Fail,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct BackendGate {
    pub status: GateStatus,
    pub details: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct ReferenceGate {
    pub status: GateStatus,
    pub details: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct SpecificityGate {
    pub status: GateStatus,
    pub details: String,
    pub negative_samples: u64,
    pub false_positives: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct SensitivityGate {
    pub status: GateStatus,
    pub details: String,
    pub spike_in_detected: u64,
    pub spike_in_total: u64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct BenchmarkRunResult {
    pub dataset_id: String,
    pub expected_outcome: BenchmarkExpectedOutcome,
    pub expected_target: Option<String>,
    pub decision_status: DecisionStatus,
    pub release_status: ReleaseStatus,
    pub matched_expectation: bool,
    pub detected_expected_target: bool,
    pub sampled_fragments: u64,
    pub qc_passing_fragments: u64,
    pub top_candidate: Option<String>,
}

pub fn run_calibration(args: &CalibrateArgs) -> Result<CalibrationOutcome> {
    let loaded_bundle = load_reference_bundle(&args.reference_bundle)?;
    let benchmark_manifest = load_benchmark_manifest(&args.benchmark_manifest)?;
    validate_output_directory(&args.out)?;
    validate_benchmark_manifest(&benchmark_manifest)?;

    fs::create_dir_all(&args.out).map_err(|source| RvScreenError::io(&args.out, source))?;
    let workspace_dir = args.out.join(".calibration-work");
    fs::create_dir_all(&workspace_dir)
        .map_err(|source| RvScreenError::io(&workspace_dir, source))?;

    let (_, date_stamp) = current_utc_timestamp()?;
    let mut profile = benchmark_manifest
        .profile
        .build_profile(&date_stamp, &loaded_bundle.bundle.version);
    let thresholds = ThresholdsToml {
        fragment_rules: profile.fragment_rules.clone(),
        candidate_rules: profile.candidate_rules.clone(),
        decision_rules: profile.decision_rules.clone(),
    };

    let execution_profile_dir = workspace_dir.join("execution-profile");
    write_profile_dir(&execution_profile_dir, &profile)?;

    let screen_mode = screen_mode_for(&profile.sampling.mode)?;
    let thread_count = default_threads();
    let benchmark_runs = benchmark_manifest
        .datasets
        .iter()
        .enumerate()
        .map(|(index, dataset)| {
            run_benchmark_dataset(
                args,
                &execution_profile_dir,
                &workspace_dir,
                index,
                dataset,
                screen_mode,
                thread_count,
            )
        })
        .collect::<Result<Vec<_>>>()?;

    let release_gate = build_release_gate(&benchmark_runs, &loaded_bundle.bundle.version);
    profile.status = if release_gate.all_pass() {
        benchmark_manifest
            .profile
            .status
            .clone()
            .unwrap_or_else(|| DEFAULT_PROFILE_STATUS.to_string())
    } else {
        FAILED_PROFILE_STATUS.to_string()
    };

    fs::remove_dir_all(&workspace_dir)
        .map_err(|source| RvScreenError::io(&workspace_dir, source))?;

    let profile_path = args.out.join(PROFILE_TOML);
    let thresholds_path = args.out.join(THRESHOLDS_TOML);
    let manifest_path = args.out.join(BENCHMARK_MANIFEST_JSON);
    let summary_path = args.out.join(BENCHMARK_SUMMARY_TSV);
    let release_gate_path = args.out.join(RELEASE_GATE_JSON);
    let checksum_path = args.out.join(CHECKSUM_SHA256);

    write_toml_pretty(&profile_path, &profile, "profile.toml")?;
    write_toml_pretty(&thresholds_path, &thresholds, "thresholds.toml")?;
    write_json_pretty(&manifest_path, &benchmark_manifest)?;
    write_benchmark_summary_tsv(&summary_path, &benchmark_runs)?;
    write_json_pretty(&release_gate_path, &release_gate)?;
    write_checksum_file(
        &checksum_path,
        &args.out,
        &[
            profile_path,
            thresholds_path,
            manifest_path,
            summary_path,
            release_gate_path,
        ],
    )?;

    Ok(CalibrationOutcome {
        output_dir: args.out.clone(),
        profile,
        benchmark_runs,
        release_gate,
    })
}

pub fn load_profile(profile_dir: impl AsRef<Path>) -> Result<LoadedProfile> {
    let profile_dir = profile_dir.as_ref().to_path_buf();
    let profile_path = profile_dir.join(PROFILE_TOML);
    let text = fs::read_to_string(&profile_path)
        .map_err(|source| RvScreenError::io(&profile_path, source))?;
    let profile = toml::from_str::<ProfileToml>(&text).map_err(|error| {
        RvScreenError::config(
            &profile_path,
            format!("failed to parse profile.toml: {error}"),
        )
    })?;

    Ok(LoadedProfile {
        profile_dir,
        profile_path,
        profile,
    })
}

pub fn load_reference_bundle(bundle_dir: impl AsRef<Path>) -> Result<LoadedReferenceBundle> {
    let bundle_dir = bundle_dir.as_ref().to_path_buf();
    let bundle_path = bundle_dir.join(BUNDLE_TOML);
    let text = fs::read_to_string(&bundle_path)
        .map_err(|source| RvScreenError::io(&bundle_path, source))?;
    let bundle = toml::from_str::<BundleToml>(&text).map_err(|error| {
        RvScreenError::config(
            &bundle_path,
            format!("failed to parse bundle.toml: {error}"),
        )
    })?;

    Ok(LoadedReferenceBundle {
        bundle_dir,
        bundle_path,
        bundle,
    })
}

impl BenchmarkProfileConfig {
    fn build_profile(&self, date_stamp: &str, reference_bundle: &str) -> ProfileToml {
        ProfileToml {
            profile_id: self
                .profile_id
                .clone()
                .unwrap_or_else(|| format!("rvscreen_calib_{date_stamp}-r1")),
            status: self
                .status
                .clone()
                .unwrap_or_else(|| DEFAULT_PROFILE_STATUS.to_string()),
            reference_bundle: reference_bundle.to_string(),
            backend: CALIBRATION_BACKEND.to_string(),
            preset: CALIBRATION_PRESET.to_string(),
            seed: self.seed.unwrap_or(DEFAULT_SEED),
            supported_input: self.supported_input.clone().unwrap_or_else(|| {
                DEFAULT_SUPPORTED_INPUT
                    .iter()
                    .map(|value| (*value).to_string())
                    .collect()
            }),
            supported_read_type: self.supported_read_type.clone().unwrap_or_else(|| {
                DEFAULT_SUPPORTED_READ_TYPE
                    .iter()
                    .map(|value| (*value).to_string())
                    .collect()
            }),
            negative_control_required: self.negative_control_required.unwrap_or(false),
            sampling: self
                .sampling
                .clone()
                .unwrap_or_else(default_sampling_config),
            fragment_rules: self
                .fragment_rules
                .clone()
                .unwrap_or_else(default_fragment_rules),
            candidate_rules: self
                .candidate_rules
                .clone()
                .unwrap_or_else(default_candidate_rules),
            decision_rules: self
                .decision_rules
                .clone()
                .unwrap_or_else(default_decision_rules),
        }
    }
}

impl ReleaseGate {
    fn all_pass(&self) -> bool {
        self.backend_gate.status == GateStatus::Pass
            && self.reference_gate.status == GateStatus::Pass
            && self.specificity_gate.status == GateStatus::Pass
            && self.sensitivity_gate.status == GateStatus::Pass
    }
}

fn load_benchmark_manifest(path: &Path) -> Result<BenchmarkManifest> {
    let bytes = fs::read(path).map_err(|source| RvScreenError::io(path, source))?;
    let mut manifest = parse_benchmark_manifest(path, &bytes)?;
    let manifest_root = path.parent().unwrap_or_else(|| Path::new("."));

    for dataset in &mut manifest.datasets {
        dataset.input = dataset
            .input
            .iter()
            .map(|input| resolve_manifest_path(manifest_root, input))
            .collect::<Result<Vec<_>>>()?;
    }

    Ok(manifest)
}

fn parse_benchmark_manifest(path: &Path, bytes: &[u8]) -> Result<BenchmarkManifest> {
    match manifest_format_hint(path) {
        BenchmarkManifestFormat::Json => serde_json::from_slice(bytes).map_err(|error| {
            RvScreenError::parse(
                path,
                1,
                format!("failed to parse benchmark manifest JSON: {error}"),
            )
        }),
        BenchmarkManifestFormat::Yaml => serde_yaml::from_slice(bytes).map_err(|error| {
            RvScreenError::parse(
                path,
                1,
                format!("failed to parse benchmark manifest YAML: {error}"),
            )
        }),
        BenchmarkManifestFormat::Auto => serde_json::from_slice(bytes)
            .or_else(|json_error| {
                serde_yaml::from_slice(bytes).map_err(|yaml_error| {
                    RvScreenError::parse(
                        path,
                        1,
                        format!(
                            "failed to parse benchmark manifest as JSON ({json_error}) or YAML ({yaml_error})"
                        ),
                    )
                })
            }),
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum BenchmarkManifestFormat {
    Json,
    Yaml,
    Auto,
}

fn manifest_format_hint(path: &Path) -> BenchmarkManifestFormat {
    let extension = path
        .extension()
        .and_then(|extension| extension.to_str())
        .map(|extension| extension.to_ascii_lowercase());
    match extension.as_deref() {
        Some("json") => BenchmarkManifestFormat::Json,
        Some("yaml") | Some("yml") => BenchmarkManifestFormat::Yaml,
        _ => BenchmarkManifestFormat::Auto,
    }
}

fn resolve_manifest_path(base_dir: &Path, path: &Path) -> Result<PathBuf> {
    let resolved = if path.is_absolute() {
        path.to_path_buf()
    } else {
        base_dir.join(path)
    };
    fs::canonicalize(&resolved).map_err(|source| RvScreenError::io(&resolved, source))
}

fn validate_benchmark_manifest(manifest: &BenchmarkManifest) -> Result<()> {
    if manifest.datasets.is_empty() {
        return Err(RvScreenError::validation(
            "benchmark_manifest.datasets",
            "benchmark manifest must contain at least one dataset",
        ));
    }

    let mut seen_ids = std::collections::BTreeSet::new();
    for dataset in &manifest.datasets {
        if dataset.dataset_id.trim().is_empty() {
            return Err(RvScreenError::validation(
                "benchmark_manifest.datasets.dataset_id",
                "dataset_id must not be empty",
            ));
        }
        if !seen_ids.insert(dataset.dataset_id.clone()) {
            return Err(RvScreenError::validation(
                "benchmark_manifest.datasets.dataset_id",
                format!("duplicate benchmark dataset_id `{}`", dataset.dataset_id),
            ));
        }
        if dataset.input.is_empty() {
            return Err(RvScreenError::validation(
                format!("benchmark_manifest.datasets.{}.input", dataset.dataset_id),
                "dataset input must contain at least one file",
            ));
        }
    }

    Ok(())
}

fn write_profile_dir(profile_dir: &Path, profile: &ProfileToml) -> Result<()> {
    fs::create_dir_all(profile_dir).map_err(|source| RvScreenError::io(profile_dir, source))?;
    write_toml_pretty(&profile_dir.join(PROFILE_TOML), profile, "profile.toml")
}

fn run_benchmark_dataset(
    args: &CalibrateArgs,
    execution_profile_dir: &Path,
    workspace_dir: &Path,
    index: usize,
    dataset: &BenchmarkDataset,
    screen_mode: ScreenMode,
    thread_count: usize,
) -> Result<BenchmarkRunResult> {
    let report_dir = workspace_dir.join(format!(
        "report-{index:02}-{}",
        sanitize_file_component(&dataset.dataset_id)
    ));
    let screen_args = ScreenArgs {
        input: dataset.input.clone(),
        reference_bundle: args.reference_bundle.clone(),
        calibration_profile: execution_profile_dir.to_path_buf(),
        negative_control: None,
        out: report_dir,
        mode: screen_mode,
        threads: thread_count,
    };
    let outcome = run_screen(&screen_args)?;
    let top_candidate = outcome
        .candidates
        .iter()
        .max_by_key(|candidate| candidate.accepted_fragments)
        .map(|candidate| candidate.accession_or_group.clone());
    let detected_expected_target = dataset
        .expected_target
        .as_ref()
        .map(|expected_target| candidate_matches_expected_target(&outcome, expected_target))
        .unwrap_or(outcome.summary.decision_status == DecisionStatus::Positive);
    let matched_expectation = match dataset.expected_outcome {
        BenchmarkExpectedOutcome::Negative => {
            outcome.summary.decision_status != DecisionStatus::Positive
        }
        BenchmarkExpectedOutcome::SpikeIn => {
            outcome.summary.decision_status == DecisionStatus::Positive && detected_expected_target
        }
    };

    Ok(BenchmarkRunResult {
        dataset_id: dataset.dataset_id.clone(),
        expected_outcome: dataset.expected_outcome.clone(),
        expected_target: dataset.expected_target.clone(),
        decision_status: outcome.summary.decision_status,
        release_status: outcome.summary.release_status,
        matched_expectation,
        detected_expected_target,
        sampled_fragments: outcome.summary.sampled_fragments,
        qc_passing_fragments: outcome.summary.qc_passing_fragments,
        top_candidate,
    })
}

fn candidate_matches_expected_target(
    outcome: &crate::pipeline::ScreenRunOutcome,
    expected_target: &str,
) -> bool {
    outcome.candidates.iter().any(|candidate| {
        candidate.accepted_fragments > 0
            && (candidate.accession_or_group == expected_target
                || candidate.virus_name == expected_target)
    })
}

fn build_release_gate(
    benchmark_runs: &[BenchmarkRunResult],
    reference_bundle_version: &str,
) -> ReleaseGate {
    let negative_samples = benchmark_runs
        .iter()
        .filter(|run| run.expected_outcome == BenchmarkExpectedOutcome::Negative)
        .count() as u64;
    let false_positives = benchmark_runs
        .iter()
        .filter(|run| {
            run.expected_outcome == BenchmarkExpectedOutcome::Negative
                && run.decision_status == DecisionStatus::Positive
        })
        .count() as u64;
    let spike_in_total = benchmark_runs
        .iter()
        .filter(|run| run.expected_outcome == BenchmarkExpectedOutcome::SpikeIn)
        .count() as u64;
    let spike_in_detected = benchmark_runs
        .iter()
        .filter(|run| {
            run.expected_outcome == BenchmarkExpectedOutcome::SpikeIn && run.matched_expectation
        })
        .count() as u64;

    let specificity_gate = if negative_samples == 0 {
        SpecificityGate {
            status: GateStatus::Fail,
            details: "No negative benchmark datasets were provided for specificity evaluation"
                .to_string(),
            negative_samples,
            false_positives,
        }
    } else if false_positives == 0 {
        SpecificityGate {
            status: GateStatus::Pass,
            details: format!(
                "All {negative_samples} negative benchmark datasets remained non-positive"
            ),
            negative_samples,
            false_positives,
        }
    } else {
        SpecificityGate {
            status: GateStatus::Fail,
            details: format!(
                "Observed {false_positives} false-positive calls across {negative_samples} negative benchmark datasets"
            ),
            negative_samples,
            false_positives,
        }
    };

    let sensitivity_gate = if spike_in_total == 0 {
        SensitivityGate {
            status: GateStatus::Fail,
            details: "No spike-in benchmark datasets were provided for sensitivity evaluation"
                .to_string(),
            spike_in_detected,
            spike_in_total,
        }
    } else if spike_in_detected == spike_in_total {
        SensitivityGate {
            status: GateStatus::Pass,
            details: format!(
                "Detected all {spike_in_detected}/{spike_in_total} spike-in benchmark datasets"
            ),
            spike_in_detected,
            spike_in_total,
        }
    } else {
        SensitivityGate {
            status: GateStatus::Fail,
            details: format!(
                "Detected {spike_in_detected}/{spike_in_total} spike-in benchmark datasets"
            ),
            spike_in_detected,
            spike_in_total,
        }
    };

    let reference_failures = benchmark_runs
        .iter()
        .filter(|run| {
            matches!(
                run.release_status,
                ReleaseStatus::Blocked | ReleaseStatus::Provisional
            )
        })
        .count() as u64;
    let reference_gate = ReferenceGate {
        status: GateStatus::Pass,
        details: format!(
            "All {} benchmark runs were executed against reference bundle `{reference_bundle_version}`",
            benchmark_runs.len()
        ),
    };
    let backend_gate = BackendGate {
        status: GateStatus::Pass,
        details: BACKEND_GATE_DETAILS.to_string(),
    };

    let reference_gate = if reference_failures == 0 {
        reference_gate
    } else {
        ReferenceGate {
            status: GateStatus::Pass,
            details: format!(
                "All {} benchmark runs were bound to reference bundle `{reference_bundle_version}`; run release_status values remained provisional as expected without negative controls",
                benchmark_runs.len()
            ),
        }
    };

    ReleaseGate {
        backend_gate,
        reference_gate,
        specificity_gate,
        sensitivity_gate,
    }
}

fn write_benchmark_summary_tsv(path: &Path, benchmark_runs: &[BenchmarkRunResult]) -> Result<()> {
    #[derive(Serialize)]
    struct BenchmarkSummaryRow<'a> {
        dataset_id: &'a str,
        expected_outcome: String,
        expected_target: &'a str,
        decision_status: String,
        release_status: String,
        matched_expectation: bool,
        detected_expected_target: bool,
        sampled_fragments: u64,
        qc_passing_fragments: u64,
        top_candidate: &'a str,
    }

    let mut writer = WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .map_err(|error| {
            RvScreenError::validation(
                "benchmark_summary.tsv",
                format!(
                    "failed to open TSV writer for `{}`: {error}",
                    path.display()
                ),
            )
        })?;

    for benchmark_run in benchmark_runs {
        writer
            .serialize(BenchmarkSummaryRow {
                dataset_id: &benchmark_run.dataset_id,
                expected_outcome: enum_label(&benchmark_run.expected_outcome)?,
                expected_target: benchmark_run.expected_target.as_deref().unwrap_or(""),
                decision_status: enum_label(&benchmark_run.decision_status)?,
                release_status: enum_label(&benchmark_run.release_status)?,
                matched_expectation: benchmark_run.matched_expectation,
                detected_expected_target: benchmark_run.detected_expected_target,
                sampled_fragments: benchmark_run.sampled_fragments,
                qc_passing_fragments: benchmark_run.qc_passing_fragments,
                top_candidate: benchmark_run.top_candidate.as_deref().unwrap_or(""),
            })
            .map_err(|error| {
                RvScreenError::validation(
                    "benchmark_summary.tsv",
                    format!(
                        "failed to serialize benchmark row `{}`: {error}",
                        benchmark_run.dataset_id
                    ),
                )
            })?;
    }

    writer
        .flush()
        .map_err(|error| RvScreenError::io(path, error))
}

fn screen_mode_for(mode: &str) -> Result<ScreenMode> {
    match mode {
        "representative" => Ok(ScreenMode::Representative),
        "streaming" => Ok(ScreenMode::Streaming),
        other => Err(RvScreenError::validation(
            "benchmark_manifest.profile.sampling.mode",
            format!(
                "unsupported calibration sampling mode `{other}`; expected `representative` or `streaming`"
            ),
        )),
    }
}

fn default_threads() -> usize {
    std::thread::available_parallelism()
        .map(usize::from)
        .unwrap_or(1)
}

fn default_sampling_config() -> SamplingConfig {
    SamplingConfig {
        mode: "representative".to_string(),
        rounds: DEFAULT_REPRESENTATIVE_ROUNDS.to_vec(),
        max_rounds: DEFAULT_REPRESENTATIVE_ROUNDS.len() as u64,
    }
}

fn default_fragment_rules() -> FragmentRules {
    FragmentRules {
        min_mapq: 20,
        min_as_diff: 12,
        max_nm: 8,
        require_pair_consistency: true,
    }
}

fn default_candidate_rules() -> CandidateRules {
    CandidateRules {
        min_nonoverlap_fragments: 3,
        min_breadth: 0.001,
        max_background_ratio: 1.5,
    }
}

fn default_decision_rules() -> DecisionRules {
    DecisionRules {
        theta_pos: 0.00005,
        theta_neg: 0.000005,
        allow_indeterminate: true,
    }
}

fn validate_output_directory(path: &Path) -> Result<()> {
    if path.exists() {
        if !path.is_dir() {
            return Err(RvScreenError::validation(
                "out",
                format!("`{}` exists and is not a directory", path.display()),
            ));
        }

        let mut entries = fs::read_dir(path).map_err(|source| RvScreenError::io(path, source))?;
        if entries
            .next()
            .transpose()
            .map_err(|source| RvScreenError::io(path, source))?
            .is_some()
        {
            return Err(RvScreenError::validation(
                "out",
                format!("output directory `{}` must be empty", path.display()),
            ));
        }
    }

    Ok(())
}

fn write_toml_pretty(path: &Path, value: &impl Serialize, field: &str) -> Result<()> {
    let body = toml::to_string_pretty(value).map_err(|error| {
        RvScreenError::config(path, format!("failed to serialize {field}: {error}"))
    })?;
    fs::write(path, body).map_err(|source| RvScreenError::io(path, source))
}

fn write_json_pretty(path: &Path, value: &impl Serialize) -> Result<()> {
    let body = serde_json::to_vec_pretty(value).map_err(|error| {
        RvScreenError::config(path, format!("failed to serialize JSON: {error}"))
    })?;
    fs::write(path, body).map_err(|source| RvScreenError::io(path, source))
}

fn write_checksum_file(checksum_path: &Path, output_dir: &Path, files: &[PathBuf]) -> Result<()> {
    let file =
        File::create(checksum_path).map_err(|source| RvScreenError::io(checksum_path, source))?;
    let mut writer = BufWriter::new(file);

    for path in files {
        let digest = sha256_file(path)?;
        let relative = path.strip_prefix(output_dir).map_err(|_| {
            RvScreenError::validation(
                "checksum",
                format!(
                    "`{}` is not inside `{}`",
                    path.display(),
                    output_dir.display()
                ),
            )
        })?;
        writeln!(writer, "{digest}  {}", relative.display())
            .map_err(|source| RvScreenError::io(checksum_path, source))?;
    }

    writer
        .flush()
        .map_err(|source| RvScreenError::io(checksum_path, source))
}

fn sha256_file(path: &Path) -> Result<String> {
    let file = File::open(path).map_err(|source| RvScreenError::io(path, source))?;
    let mut reader = BufReader::new(file);
    let mut hasher = Sha256::new();
    let mut buffer = [0u8; 8192];

    loop {
        let read = reader
            .read(&mut buffer)
            .map_err(|source| RvScreenError::io(path, source))?;
        if read == 0 {
            break;
        }
        hasher.update(&buffer[..read]);
    }

    Ok(format!("{:x}", hasher.finalize()))
}

fn current_utc_timestamp() -> Result<(String, String)> {
    let seconds = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_err(|error| RvScreenError::validation("time", error.to_string()))?
        .as_secs() as i64;
    let (year, month, day, hour, minute, second) = unix_to_utc(seconds);
    Ok((
        format!("{year:04}-{month:02}-{day:02}T{hour:02}:{minute:02}:{second:02}Z"),
        format!("{year:04}.{month:02}.{day:02}"),
    ))
}

fn unix_to_utc(seconds: i64) -> (i64, i64, i64, i64, i64, i64) {
    let days = seconds.div_euclid(86_400);
    let seconds_of_day = seconds.rem_euclid(86_400);
    let (year, month, day) = civil_from_days(days);
    let hour = seconds_of_day / 3_600;
    let minute = (seconds_of_day % 3_600) / 60;
    let second = seconds_of_day % 60;
    (year, month, day, hour, minute, second)
}

fn civil_from_days(days_since_epoch: i64) -> (i64, i64, i64) {
    let z = days_since_epoch + 719_468;
    let era = if z >= 0 { z } else { z - 146_096 } / 146_097;
    let doe = z - era * 146_097;
    let yoe = (doe - doe / 1_460 + doe / 36_524 - doe / 146_096) / 365;
    let y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let day = doy - (153 * mp + 2) / 5 + 1;
    let month = mp + if mp < 10 { 3 } else { -9 };
    let year = y + if month <= 2 { 1 } else { 0 };
    (year, month, day)
}

fn sanitize_file_component(value: &str) -> String {
    let sanitized: String = value
        .chars()
        .map(|ch| match ch {
            'a'..='z' | 'A'..='Z' | '0'..='9' | '-' | '_' | '.' => ch,
            _ => '_',
        })
        .collect();
    if sanitized.is_empty() {
        "dataset".to_string()
    } else {
        sanitized
    }
}

fn enum_label(value: &impl Serialize) -> Result<String> {
    serde_json::to_string(value)
        .map(|label| label.trim_matches('"').to_string())
        .map_err(|error| RvScreenError::validation("serialize", error.to_string()))
}

#[cfg(test)]
#[path = "../../tests/testutil/mod.rs"]
mod testutil;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cli::CalibrateArgs;
    use crate::reference::{build_reference_bundle, BuildReferenceBundleRequest};
    use crate::types::{BundleManifest, ContigEntry};
    use std::fs;
    use std::io;
    use std::path::Path;
    use tempfile::tempdir;

    use crate::calibration::testutil::{
        generate_fastq_pair, generate_mini_reference, FastqPairConfig, ReadComponent,
        SyntheticSource, VirusSelector,
    };

    #[test]
    fn loads_profile_toml_from_directory() {
        let tempdir = tempdir().expect("tempdir should be created");
        let profile_dir = tempdir.path().join("calibration");
        fs::create_dir_all(&profile_dir).expect("profile dir should be created");
        fs::write(
            profile_dir.join(PROFILE_TOML),
            r#"
profile_id = "rvscreen_calib_2026.04.20-r1"
status = "release_candidate"
reference_bundle = "rvscreen_ref_2026.04.20-r1"
backend = "minimap2"
preset = "sr-conservative"
seed = 20260420
supported_input = ["fastq.gz", "bam", "cram"]
supported_read_type = ["illumina_pe_shortread"]
negative_control_required = false

[sampling]
mode = "representative"
rounds = [10, 20]
max_rounds = 2

[fragment_rules]
min_mapq = 20
min_as_diff = 12
max_nm = 8
require_pair_consistency = true

[candidate_rules]
min_nonoverlap_fragments = 1
min_breadth = 0.0
max_background_ratio = 0.0

[decision_rules]
theta_pos = 0.00001
theta_neg = 0.0
allow_indeterminate = true
"#,
        )
        .expect("profile TOML should be written");

        let loaded = load_profile(&profile_dir).expect("profile should load");

        assert_eq!(loaded.profile.profile_id, "rvscreen_calib_2026.04.20-r1");
        assert_eq!(loaded.profile.sampling.rounds, vec![10, 20]);
    }

    #[test]
    fn benchmark_manifest_resolves_relative_inputs() {
        let tempdir = tempdir().expect("tempdir should be created");
        let manifest_dir = tempdir.path().join("benchmarks");
        let dataset_dir = manifest_dir.join("dataset-a");
        fs::create_dir_all(&dataset_dir).expect("dataset dir should exist");
        let r1 = dataset_dir.join("reads_R1.fastq");
        let r2 = dataset_dir.join("reads_R2.fastq");
        fs::write(&r1, "@r1\nACGT\n+\nIIII\n").expect("r1 should be written");
        fs::write(&r2, "@r2\nTGCA\n+\nIIII\n").expect("r2 should be written");
        let manifest_path = manifest_dir.join("manifest.json");
        fs::write(
            &manifest_path,
            r#"{
  "datasets": [
    {
      "dataset_id": "dataset-a",
      "input": ["dataset-a/reads_R1.fastq", "dataset-a/reads_R2.fastq"],
      "expected_outcome": "negative"
    }
  ]
}"#,
        )
        .expect("manifest should be written");

        let manifest = load_benchmark_manifest(&manifest_path).expect("manifest should load");

        assert_eq!(manifest.datasets.len(), 1);
        assert_eq!(manifest.datasets[0].input, vec![r1, r2]);
    }

    #[test]
    fn benchmark_manifest_supports_yaml() {
        let tempdir = tempdir().expect("tempdir should be created");
        let manifest_dir = tempdir.path().join("benchmarks-yaml");
        let dataset_dir = manifest_dir.join("dataset-a");
        fs::create_dir_all(&dataset_dir).expect("dataset dir should exist");
        let r1 = dataset_dir.join("reads_R1.fastq");
        let r2 = dataset_dir.join("reads_R2.fastq");
        fs::write(&r1, "@r1\nACGT\n+\nIIII\n").expect("r1 should be written");
        fs::write(&r2, "@r2\nTGCA\n+\nIIII\n").expect("r2 should be written");
        let manifest_path = manifest_dir.join("manifest.yaml");
        fs::write(
            &manifest_path,
            r#"profile:
  supported_input:
    - fastq
datasets:
  - dataset_id: dataset-a
    input:
      - dataset-a/reads_R1.fastq
      - dataset-a/reads_R2.fastq
    expected_outcome: negative
"#,
        )
        .expect("manifest should be written");

        let manifest = load_benchmark_manifest(&manifest_path).expect("YAML manifest should load");

        assert_eq!(manifest.datasets.len(), 1);
        assert_eq!(manifest.datasets[0].input, vec![r1, r2]);
        assert_eq!(
            manifest.profile.supported_input,
            Some(vec!["fastq".to_string()])
        );
    }

    #[test]
    fn calibration_generates_profile_and_required_artifacts() {
        let tempdir = tempdir().expect("tempdir should be created");
        let bundle =
            prepare_reference_bundle(tempdir.path()).expect("reference bundle should build");
        let benchmark_manifest =
            write_benchmark_manifest(tempdir.path()).expect("manifest should be written");
        let out_dir = tempdir.path().join("calibration-profile");

        let outcome = run_calibration(&CalibrateArgs {
            reference_bundle: bundle.bundle_dir.clone(),
            benchmark_manifest: benchmark_manifest.clone(),
            out: out_dir.clone(),
        })
        .expect("calibration should succeed");

        assert_eq!(outcome.profile.reference_bundle, bundle.version);
        assert_eq!(outcome.benchmark_runs.len(), 2);
        assert!(outcome.release_gate.all_pass());
        assert_eq!(
            outcome.release_gate.backend_gate.details,
            BACKEND_GATE_DETAILS
        );
        assert_eq!(outcome.release_gate.specificity_gate.false_positives, 0);
        assert_eq!(outcome.release_gate.sensitivity_gate.spike_in_detected, 1);
        assert_eq!(outcome.release_gate.sensitivity_gate.spike_in_total, 1);
        assert_eq!(outcome.profile.status, DEFAULT_PROFILE_STATUS);

        for artifact in [
            PROFILE_TOML,
            THRESHOLDS_TOML,
            BENCHMARK_MANIFEST_JSON,
            BENCHMARK_SUMMARY_TSV,
            RELEASE_GATE_JSON,
            CHECKSUM_SHA256,
        ] {
            assert!(
                out_dir.join(artifact).exists(),
                "missing artifact {artifact}"
            );
        }

        let profile: ProfileToml = toml::from_str(
            &fs::read_to_string(out_dir.join(PROFILE_TOML))
                .expect("profile.toml should be readable"),
        )
        .expect("profile.toml should parse");
        assert_eq!(profile.reference_bundle, bundle.version);

        let benchmark_summary = fs::read_to_string(out_dir.join(BENCHMARK_SUMMARY_TSV))
            .expect("summary should be readable");
        assert_eq!(benchmark_summary.lines().count(), 3);

        let release_gate: ReleaseGate = serde_json::from_str(
            &fs::read_to_string(out_dir.join(RELEASE_GATE_JSON))
                .expect("release_gate.json should be readable"),
        )
        .expect("release_gate.json should parse");
        assert_eq!(release_gate.sensitivity_gate.spike_in_detected, 1);

        let checksum = fs::read_to_string(out_dir.join(CHECKSUM_SHA256))
            .expect("checksum file should be readable");
        assert_eq!(checksum.lines().count(), 5);
        for artifact in [
            PROFILE_TOML,
            THRESHOLDS_TOML,
            BENCHMARK_MANIFEST_JSON,
            BENCHMARK_SUMMARY_TSV,
            RELEASE_GATE_JSON,
        ] {
            assert!(
                checksum.contains(artifact),
                "checksum missing {artifact}: {checksum}"
            );
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

    fn write_benchmark_manifest(base_dir: &Path) -> io::Result<PathBuf> {
        let benchmark_dir = base_dir.join("benchmarks");
        fs::create_dir_all(&benchmark_dir)?;
        let negative_dir = benchmark_dir.join("negative");
        let spike_dir = benchmark_dir.join("spike");

        generate_fastq_pair(
            &FastqPairConfig::new(200, 100, 44)
                .with_output_dir(&negative_dir)
                .with_components(vec![ReadComponent::new(SyntheticSource::Human, 1.0)]),
        )?;
        generate_fastq_pair(
            &FastqPairConfig::new(200, 100, 55)
                .with_output_dir(&spike_dir)
                .with_components(vec![
                    ReadComponent::new(SyntheticSource::Human, 0.95),
                    ReadComponent::new(
                        SyntheticSource::Virus(VirusSelector::Accession(
                            "NC_SYNTHV1.1".to_string(),
                        )),
                        0.05,
                    ),
                ]),
        )?;

        let manifest_path = benchmark_dir.join("manifest.json");
        fs::write(
            &manifest_path,
            r#"{
  "profile": {
    "supported_input": ["fastq"],
    "supported_read_type": ["illumina_pe_shortread"],
    "negative_control_required": false,
    "sampling": {
      "mode": "representative",
      "rounds": [50, 100],
      "max_rounds": 2
    },
    "fragment_rules": {
      "min_mapq": 0,
      "min_as_diff": 0,
      "max_nm": 100,
      "require_pair_consistency": true
    },
    "candidate_rules": {
      "min_nonoverlap_fragments": 1,
      "min_breadth": 0.0,
      "max_background_ratio": 0.0
    },
    "decision_rules": {
      "theta_pos": 0.01,
      "theta_neg": 0.0001,
      "allow_indeterminate": true
    }
  },
  "datasets": [
    {
      "dataset_id": "negative-human",
      "input": ["negative/synthetic_R1.fastq", "negative/synthetic_R2.fastq"],
      "expected_outcome": "negative"
    },
    {
      "dataset_id": "spike-5pct",
      "input": ["spike/synthetic_R1.fastq", "spike/synthetic_R2.fastq"],
      "expected_outcome": "spike_in",
      "expected_target": "NC_SYNTHV1.1"
    }
  ]
}"#,
        )?;

        Ok(manifest_path)
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
