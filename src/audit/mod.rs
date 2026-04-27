use crate::calibration::{load_profile, load_reference_bundle};
use crate::cli::AuditVerifyArgs;
use crate::decision::{
    expected_release_status, NegativeControlDecisionInput, NegativeControlResult,
};
use crate::error::{Result, RvScreenError};
use crate::report::contract::{
    sampling_only_ci_label_violation, ReportArtifactKind, CANDIDATE_CALLS_TSV,
    CANONICAL_CANDIDATE_CALLS_ACCESSION_TSV, CANONICAL_CHECKSUM_SHA256,
    CANONICAL_RESULT_OVERVIEW_JSON, CANONICAL_RUN_MANIFEST_JSON, CANONICAL_SAMPLE_RUN_SUMMARY_JSON,
    CANONICAL_SAMPLING_ROUNDS_TSV, CANONICAL_SCHEMA_VERSIONS_JSON, CANONICAL_VIRUS_SUMMARY_TSV,
    CHECKSUM_SHA256, REPORT_ARTIFACT_CONTRACTS, REQUIRED_REPORT_BUNDLE_FILES, ROUNDS_TSV,
    RUN_MANIFEST_JSON, SAMPLE_SUMMARY_JSON,
};
use crate::report::summary::VIRUS_SUMMARY_HEADER;
use crate::types::{
    DecisionStatus, NegativeControlStatus, ProfileToml, ReleaseStatus, RoundRecord, RunManifest,
    SampleSummary, StopReason,
};
use csv::ReaderBuilder;
use serde::Deserialize;
use sha2::{Digest, Sha256};
use std::collections::BTreeMap;
use std::fs::{self, File};
use std::io::{BufReader, Read};
use std::path::{Component, Path, PathBuf};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AuditVerifyReport {
    bundle_dir: PathBuf,
    checks: Vec<AuditCheck>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct AuditCheck {
    label: &'static str,
    status: AuditCheckStatus,
    detail: String,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum AuditCheckStatus {
    Pass,
    Fail,
    Skip,
}

#[derive(Debug, Deserialize)]
struct CandidateCallTsvRow {
    accession_or_group: String,
    fraction_ci_95: String,
    decision_reasons: String,
}

#[derive(Debug, Clone, PartialEq)]
struct CandidateCiRecord {
    accession_or_group: String,
    fraction_ci_95: [f64; 2],
    decision_reasons: Vec<String>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct ReportBundleLayout {
    has_v2_artifacts: bool,
}

#[derive(Debug, Deserialize)]
struct SchemaVersions {
    report_bundle_schema: String,
    canonical_outputs: BTreeMap<String, String>,
    legacy_outputs: BTreeMap<String, String>,
}

pub fn run_audit_verify(args: &AuditVerifyArgs) -> Result<AuditVerifyReport> {
    let mut report = AuditVerifyReport::new(args.report_bundle.clone());
    let layout = detect_report_bundle_layout(&args.report_bundle);

    verify_legacy_completeness(&args.report_bundle, &mut report);
    verify_canonical_completeness(&args.report_bundle, layout, &mut report);
    verify_schema_versions(&args.report_bundle, layout, &mut report);
    verify_duplicate_layout_artifacts(&args.report_bundle, layout, &mut report);
    verify_checksum_files(&args.report_bundle, layout, &mut report);

    let candidate_calls_path = if layout.has_v2_artifacts
        && args
            .report_bundle
            .join(CANONICAL_CANDIDATE_CALLS_ACCESSION_TSV)
            .is_file()
    {
        args.report_bundle
            .join(CANONICAL_CANDIDATE_CALLS_ACCESSION_TSV)
    } else {
        args.report_bundle.join(CANDIDATE_CALLS_TSV)
    };
    let rounds_path = if layout.has_v2_artifacts
        && args
            .report_bundle
            .join(CANONICAL_SAMPLING_ROUNDS_TSV)
            .is_file()
    {
        args.report_bundle.join(CANONICAL_SAMPLING_ROUNDS_TSV)
    } else {
        args.report_bundle.join(ROUNDS_TSV)
    };
    let run_manifest_path = if layout.has_v2_artifacts
        && args
            .report_bundle
            .join(CANONICAL_RUN_MANIFEST_JSON)
            .is_file()
    {
        args.report_bundle.join(CANONICAL_RUN_MANIFEST_JSON)
    } else {
        args.report_bundle.join(RUN_MANIFEST_JSON)
    };

    let summary = parse_json_file::<SampleSummary>(
        &args.report_bundle.join(SAMPLE_SUMMARY_JSON),
        "sample_summary.json",
        &mut report,
        |summary| {
            format!(
                "parsed sample `{}` with release_status `{}`",
                summary.sample_id,
                release_status_label(&summary.release_status)
            )
        },
    );
    let manifest = parse_json_file::<RunManifest>(
        &run_manifest_path,
        "run_manifest.json",
        &mut report,
        |manifest| {
            format!(
                "parsed backend `{}` bound to reference bundle `{}`",
                manifest.backend, manifest.reference_bundle_version
            )
        },
    );
    let rounds = parse_rounds_file(&rounds_path, &mut report);
    let candidates = parse_candidate_calls_file(&candidate_calls_path, &mut report);

    verify_canonical_summary_artifacts(
        &args.report_bundle,
        layout,
        summary.as_ref(),
        manifest.as_ref(),
        rounds.as_deref(),
        &mut report,
    );

    if let (Some(summary), Some(manifest)) = (summary.as_ref(), manifest.as_ref()) {
        verify_internal_bindings(summary, manifest, &mut report);
    } else {
        report.skip(
            "report bundle bindings",
            "skipped because sample_summary.json or run_manifest.json did not parse",
        );
    }

    if let (Some(summary), Some(rounds)) = (summary.as_ref(), rounds.as_ref()) {
        if summary.rounds_run == rounds.len() as u64 {
            report.pass(
                "rounds.tsv",
                format!(
                    "parsed {} round record(s) matching sample_summary.rounds_run",
                    rounds.len()
                ),
            );
        } else {
            report.fail(
                "rounds.tsv",
                format!(
                    "sample_summary.rounds_run `{}` does not match rounds.tsv record count `{}`",
                    summary.rounds_run,
                    rounds.len()
                ),
            );
        }
    } else {
        report.skip(
            "rounds.tsv",
            "skipped round-count validation because sample_summary.json or rounds.tsv did not parse",
        );
    }

    if let Some(candidates) = candidates.as_ref() {
        verify_candidate_ci_labels(candidates, &mut report);
    } else {
        report.skip(
            "candidate CI labeling",
            "skipped because candidate_calls.tsv did not parse",
        );
    }

    if let Some(reference_bundle_dir) = args.reference_bundle.as_ref() {
        verify_reference_bundle_binding(
            reference_bundle_dir,
            summary.as_ref(),
            manifest.as_ref(),
            &mut report,
        );
    } else {
        report.skip(
            "reference bundle binding",
            "skipped because --reference-bundle was not provided",
        );
    }

    let loaded_profile = if let Some(calibration_profile_dir) = args.calibration_profile.as_ref() {
        verify_calibration_profile_binding(
            calibration_profile_dir,
            summary.as_ref(),
            manifest.as_ref(),
            &mut report,
        )
    } else {
        report.skip(
            "calibration profile binding",
            "skipped because --calibration-profile was not provided",
        );
        None
    };

    verify_release_status(
        summary.as_ref(),
        manifest.as_ref(),
        loaded_profile.as_ref(),
        args.calibration_profile.as_ref(),
        &mut report,
    );

    Ok(report)
}

impl AuditVerifyReport {
    fn new(bundle_dir: PathBuf) -> Self {
        Self {
            bundle_dir,
            checks: Vec::new(),
        }
    }

    pub fn passed(&self) -> bool {
        self.checks
            .iter()
            .all(|check| check.status != AuditCheckStatus::Fail)
    }

    pub fn render(&self) -> String {
        let overall = if self.passed() { "PASS" } else { "FAIL" };
        let mut body = vec![format!(
            "{overall} rvscreen audit verify ({})",
            self.bundle_dir.display()
        )];
        body.extend(self.checks.iter().map(|check| {
            format!(
                "[{}] {}: {}",
                check.status.label(),
                check.label,
                check.detail
            )
        }));
        body.join("\n")
    }

    fn push(&mut self, label: &'static str, status: AuditCheckStatus, detail: impl Into<String>) {
        self.checks.push(AuditCheck {
            label,
            status,
            detail: detail.into(),
        });
    }

    fn pass(&mut self, label: &'static str, detail: impl Into<String>) {
        self.push(label, AuditCheckStatus::Pass, detail);
    }

    fn fail(&mut self, label: &'static str, detail: impl Into<String>) {
        self.push(label, AuditCheckStatus::Fail, detail);
    }

    fn skip(&mut self, label: &'static str, detail: impl Into<String>) {
        self.push(label, AuditCheckStatus::Skip, detail);
    }
}

impl AuditCheckStatus {
    fn label(self) -> &'static str {
        match self {
            Self::Pass => "PASS",
            Self::Fail => "FAIL",
            Self::Skip => "SKIP",
        }
    }
}

fn detect_report_bundle_layout(report_bundle: &Path) -> ReportBundleLayout {
    let has_v2_artifacts = REPORT_ARTIFACT_CONTRACTS
        .iter()
        .any(|contract| report_bundle.join(contract.canonical_path).exists());

    ReportBundleLayout { has_v2_artifacts }
}

fn verify_legacy_completeness(report_bundle: &Path, report: &mut AuditVerifyReport) {
    let missing_files = REQUIRED_REPORT_BUNDLE_FILES
        .iter()
        .filter(|name| !report_bundle.join(name).is_file())
        .copied()
        .collect::<Vec<_>>();
    if missing_files.is_empty() {
        report.pass(
            "report bundle completeness",
            format!(
                "found required files: {}",
                REQUIRED_REPORT_BUNDLE_FILES.join(", ")
            ),
        );
    } else {
        report.fail(
            "report bundle completeness",
            format!("missing required file(s): {}", missing_files.join(", ")),
        );
    }
}

fn verify_canonical_completeness(
    report_bundle: &Path,
    layout: ReportBundleLayout,
    report: &mut AuditVerifyReport,
) {
    if !layout.has_v2_artifacts {
        report.skip(
            "canonical report bundle completeness",
            "legacy-only bundle; canonical v2 artifacts are not required",
        );
        return;
    }

    let missing_artifacts = REPORT_ARTIFACT_CONTRACTS
        .iter()
        .filter(|contract| {
            let path = report_bundle.join(contract.canonical_path);
            match contract.kind {
                ReportArtifactKind::File => !path.is_file(),
                ReportArtifactKind::Directory => !path.is_dir(),
            }
        })
        .map(|contract| contract.canonical_path)
        .collect::<Vec<_>>();

    if missing_artifacts.is_empty() {
        report.pass(
            "canonical report bundle completeness",
            format!(
                "found mandatory canonical v2 artifacts: {}",
                REPORT_ARTIFACT_CONTRACTS
                    .iter()
                    .map(|contract| contract.canonical_path)
                    .collect::<Vec<_>>()
                    .join(", ")
            ),
        );
    } else {
        report.fail(
            "canonical report bundle completeness",
            format!(
                "missing canonical mandatory artifact(s): {}",
                missing_artifacts.join(", ")
            ),
        );
    }
}

fn verify_schema_versions(
    report_bundle: &Path,
    layout: ReportBundleLayout,
    report: &mut AuditVerifyReport,
) {
    if !layout.has_v2_artifacts {
        report.skip(
            "schema metadata",
            "legacy-only bundle; audit/schema_versions.json is not required",
        );
        return;
    }

    let schema_path = report_bundle.join(CANONICAL_SCHEMA_VERSIONS_JSON);
    let Some(schema_versions) = parse_json_file::<SchemaVersions>(
        &schema_path,
        "schema metadata",
        report,
        |schema_versions| {
            format!(
                "parsed report bundle schema `{}`",
                schema_versions.report_bundle_schema
            )
        },
    ) else {
        return;
    };

    let expected_canonical_outputs = REPORT_ARTIFACT_CONTRACTS
        .iter()
        .map(|contract| {
            (
                contract.key.to_string(),
                contract.canonical_path.to_string(),
            )
        })
        .collect::<BTreeMap<_, _>>();
    let expected_legacy_outputs = REPORT_ARTIFACT_CONTRACTS
        .iter()
        .filter_map(|contract| {
            contract
                .legacy_alias_path
                .map(|legacy_alias_path| (contract.key.to_string(), legacy_alias_path.to_string()))
        })
        .collect::<BTreeMap<_, _>>();

    let mut mismatches = Vec::new();
    if schema_versions.report_bundle_schema != "rvscreen.report_bundle.v2" {
        mismatches.push(format!(
            "report_bundle_schema expected `rvscreen.report_bundle.v2`, observed `{}`",
            schema_versions.report_bundle_schema
        ));
    }
    if schema_versions.canonical_outputs != expected_canonical_outputs {
        mismatches.push("canonical_outputs does not match REPORT_ARTIFACT_CONTRACTS".to_string());
    }
    if schema_versions.legacy_outputs != expected_legacy_outputs {
        mismatches.push("legacy_outputs does not match REPORT_ARTIFACT_CONTRACTS".to_string());
    }

    if mismatches.is_empty() {
        report.pass(
            "schema metadata contract",
            "audit/schema_versions.json matches REPORT_ARTIFACT_CONTRACTS",
        );
    } else {
        report.fail("schema metadata contract", mismatches.join("; "));
    }
}

fn verify_duplicate_layout_artifacts(
    report_bundle: &Path,
    layout: ReportBundleLayout,
    report: &mut AuditVerifyReport,
) {
    if !layout.has_v2_artifacts {
        report.skip(
            "dual-layout duplicate artifacts",
            "legacy-only bundle; no canonical aliases to cross-check",
        );
        return;
    }

    compare_duplicate_file_bytes(
        report_bundle,
        CANDIDATE_CALLS_TSV,
        CANONICAL_CANDIDATE_CALLS_ACCESSION_TSV,
        "dual-layout candidate TSV",
        "candidate TSV mismatch",
        report,
    );
    compare_duplicate_file_bytes(
        report_bundle,
        ROUNDS_TSV,
        CANONICAL_SAMPLING_ROUNDS_TSV,
        "dual-layout rounds TSV",
        "rounds TSV mismatch",
        report,
    );
    compare_duplicate_json_values(
        report_bundle,
        RUN_MANIFEST_JSON,
        CANONICAL_RUN_MANIFEST_JSON,
        "dual-layout run manifest",
        "run manifest mismatch",
        report,
    );
    compare_checksum_inventory(report_bundle, report);
}

fn compare_duplicate_file_bytes(
    report_bundle: &Path,
    legacy_path: &str,
    canonical_path: &str,
    label: &'static str,
    mismatch_message: &'static str,
    report: &mut AuditVerifyReport,
) {
    let legacy = report_bundle.join(legacy_path);
    let canonical = report_bundle.join(canonical_path);
    if !legacy.is_file() || !canonical.is_file() {
        report.skip(
            label,
            format!(
                "skipped duplicate comparison because `{}` or `{}` is missing",
                legacy_path, canonical_path
            ),
        );
        return;
    }

    match (fs::read(&legacy), fs::read(&canonical)) {
        (Ok(legacy_bytes), Ok(canonical_bytes)) if legacy_bytes == canonical_bytes => {
            report.pass(label, format!("`{legacy_path}` matches `{canonical_path}`"))
        }
        (Ok(_), Ok(_)) => report.fail(
            label,
            format!("{mismatch_message} between `{legacy_path}` and `{canonical_path}`"),
        ),
        (Err(error), _) => report.fail(
            label,
            format!("failed to read `{}`: {error}", legacy.display()),
        ),
        (_, Err(error)) => report.fail(
            label,
            format!("failed to read `{}`: {error}", canonical.display()),
        ),
    }
}

fn compare_duplicate_json_values(
    report_bundle: &Path,
    legacy_path: &str,
    canonical_path: &str,
    label: &'static str,
    mismatch_message: &'static str,
    report: &mut AuditVerifyReport,
) {
    let legacy = report_bundle.join(legacy_path);
    let canonical = report_bundle.join(canonical_path);
    if !legacy.is_file() || !canonical.is_file() {
        report.skip(
            label,
            format!(
                "skipped duplicate comparison because `{}` or `{}` is missing",
                legacy_path, canonical_path
            ),
        );
        return;
    }

    let legacy_value = read_json_value(&legacy);
    let canonical_value = read_json_value(&canonical);
    match (legacy_value, canonical_value) {
        (Ok(legacy_value), Ok(canonical_value)) if legacy_value == canonical_value => {
            report.pass(label, format!("`{legacy_path}` matches `{canonical_path}`"))
        }
        (Ok(_), Ok(_)) => report.fail(
            label,
            format!("{mismatch_message} between `{legacy_path}` and `{canonical_path}`"),
        ),
        (Err(error), _) => report.fail(label, error),
        (_, Err(error)) => report.fail(label, error),
    }
}

fn compare_checksum_inventory(report_bundle: &Path, report: &mut AuditVerifyReport) {
    let legacy = report_bundle.join(CHECKSUM_SHA256);
    let canonical = report_bundle.join(CANONICAL_CHECKSUM_SHA256);
    if !legacy.is_file() || !canonical.is_file() {
        report.skip(
            "dual-layout checksum inventory",
            format!(
                "skipped checksum inventory comparison because `{}` or `{}` is missing",
                CHECKSUM_SHA256, CANONICAL_CHECKSUM_SHA256
            ),
        );
        return;
    }

    match (
        checksum_inventory(&legacy),
        checksum_inventory(&canonical),
    ) {
        (Ok(legacy_inventory), Ok(canonical_inventory)) if legacy_inventory == canonical_inventory => {
            report.pass(
                "dual-layout checksum inventory",
                "checksum.sha256 and audit/checksum.sha256 list identical covered paths and digests",
            )
        }
        (Ok(_), Ok(_)) => report.fail(
            "dual-layout checksum inventory",
            "checksum-covered inventory mismatch between `checksum.sha256` and `audit/checksum.sha256`",
        ),
        (Err(error), _) => report.fail("dual-layout checksum inventory", error.to_string()),
        (_, Err(error)) => report.fail("dual-layout checksum inventory", error.to_string()),
    }
}

fn verify_checksum_files(
    report_bundle: &Path,
    layout: ReportBundleLayout,
    report: &mut AuditVerifyReport,
) {
    verify_checksum_if_present(report_bundle, CHECKSUM_SHA256, "checksum.sha256", report);
    if layout.has_v2_artifacts {
        verify_checksum_if_present(
            report_bundle,
            CANONICAL_CHECKSUM_SHA256,
            "audit/checksum.sha256",
            report,
        );
    }
}

fn verify_checksum_if_present(
    report_bundle: &Path,
    relative_checksum_path: &str,
    label: &'static str,
    report: &mut AuditVerifyReport,
) {
    let checksum_path = report_bundle.join(relative_checksum_path);
    if !checksum_path.exists() {
        report.skip(
            label,
            format!("{relative_checksum_path} not present; checksum validation skipped"),
        );
        return;
    }

    let result = verify_checksum_file(&checksum_path, report_bundle);
    match result {
        Ok(verified_files) => report.pass(
            label,
            format!(
                "verified {verified_files} checksum entr{}",
                if verified_files == 1 { "y" } else { "ies" }
            ),
        ),
        Err(error) => report.fail(label, error.to_string()),
    }
}

fn verify_checksum_file(checksum_path: &Path, report_bundle: &Path) -> Result<usize> {
    let content = fs::read_to_string(checksum_path)
        .map_err(|source| RvScreenError::io(checksum_path, source))?;
    let mut verified_files = 0usize;

    for (index, raw_line) in content.lines().enumerate() {
        let line_number = index + 1;
        let line = raw_line.trim();
        if line.is_empty() {
            continue;
        }

        let (expected_digest, relative_path) = line.split_once("  ").ok_or_else(|| {
            RvScreenError::parse(
                checksum_path,
                line_number as u64,
                "expected `<sha256>  <relative-path>` entry",
            )
        })?;

        if expected_digest.len() != 64 || !expected_digest.chars().all(|ch| ch.is_ascii_hexdigit())
        {
            return Err(RvScreenError::parse(
                checksum_path,
                line_number as u64,
                format!("invalid SHA-256 digest `{expected_digest}`"),
            ));
        }

        let relative_path = Path::new(relative_path);
        if relative_path.is_absolute()
            || relative_path
                .components()
                .any(|component| matches!(component, Component::ParentDir | Component::RootDir))
        {
            return Err(RvScreenError::parse(
                checksum_path,
                line_number as u64,
                format!(
                    "checksum entry `{}` must stay within the report bundle",
                    relative_path.display()
                ),
            ));
        }

        let file_path = report_bundle.join(relative_path);
        if !file_path.is_file() {
            return Err(RvScreenError::validation(
                format!("checksum.{}", relative_path.display()),
                format!("listed file `{}` is missing", file_path.display()),
            ));
        }

        let observed_digest = sha256_file(&file_path)?;
        if observed_digest != expected_digest.to_ascii_lowercase() {
            return Err(RvScreenError::validation(
                format!("checksum.{}", relative_path.display()),
                format!(
                    "checksum mismatch for `{}`: expected {expected_digest}, observed {observed_digest}",
                    relative_path.display()
                ),
            ));
        }

        verified_files += 1;
    }

    Ok(verified_files)
}

fn checksum_inventory(checksum_path: &Path) -> Result<Vec<(String, String)>> {
    let content = fs::read_to_string(checksum_path)
        .map_err(|source| RvScreenError::io(checksum_path, source))?;
    let mut entries = Vec::new();

    for (index, raw_line) in content.lines().enumerate() {
        let line_number = index + 1;
        let line = raw_line.trim();
        if line.is_empty() {
            continue;
        }

        let (expected_digest, relative_path) = line.split_once("  ").ok_or_else(|| {
            RvScreenError::parse(
                checksum_path,
                line_number as u64,
                "expected `<sha256>  <relative-path>` entry",
            )
        })?;
        if expected_digest.len() != 64 || !expected_digest.chars().all(|ch| ch.is_ascii_hexdigit())
        {
            return Err(RvScreenError::parse(
                checksum_path,
                line_number as u64,
                format!("invalid SHA-256 digest `{expected_digest}`"),
            ));
        }
        entries.push((
            relative_path.replace('\\', "/"),
            expected_digest.to_ascii_lowercase(),
        ));
    }

    entries.sort();
    Ok(entries)
}

fn verify_canonical_summary_artifacts(
    report_bundle: &Path,
    layout: ReportBundleLayout,
    summary: Option<&SampleSummary>,
    manifest: Option<&RunManifest>,
    rounds: Option<&[RoundRecord]>,
    report: &mut AuditVerifyReport,
) {
    if !layout.has_v2_artifacts {
        report.skip(
            "canonical summary artifacts",
            "legacy-only bundle; canonical summary artifacts are not required",
        );
        return;
    }

    let sample_run_summary = parse_json_value_file(
        &report_bundle.join(CANONICAL_SAMPLE_RUN_SUMMARY_JSON),
        "summary/sample_run_summary.json",
        report,
    );
    let result_overview = parse_json_value_file(
        &report_bundle.join(CANONICAL_RESULT_OVERVIEW_JSON),
        "summary/result_overview.json",
        report,
    );
    verify_virus_summary_tsv(
        &report_bundle.join(CANONICAL_VIRUS_SUMMARY_TSV),
        "summary/virus_summary.tsv",
        report,
    );

    if let (Some(summary), Some(sample_run_summary)) = (summary, sample_run_summary.as_ref()) {
        verify_sample_run_summary_facts(sample_run_summary, summary, rounds, report);
    } else {
        report.skip(
            "canonical sample summary facts",
            "skipped because sample_summary.json or summary/sample_run_summary.json did not parse",
        );
    }

    if let (Some(summary), Some(manifest), Some(result_overview)) =
        (summary, manifest, result_overview.as_ref())
    {
        verify_result_overview_facts(result_overview, summary, manifest, report);
    } else {
        report.skip(
            "canonical result overview facts",
            "skipped because sample_summary.json, run_manifest.json, or summary/result_overview.json did not parse",
        );
    }
}

fn parse_json_value_file(
    path: &Path,
    label: &'static str,
    report: &mut AuditVerifyReport,
) -> Option<serde_json::Value> {
    match read_json_value(path) {
        Ok(value) => {
            report.pass(label, format!("parsed `{}`", path.display()));
            Some(value)
        }
        Err(error) => {
            report.fail(label, error);
            None
        }
    }
}

fn read_json_value(path: &Path) -> std::result::Result<serde_json::Value, String> {
    let content = fs::read_to_string(path)
        .map_err(|error| format!("failed to read `{}`: {error}", path.display()))?;
    serde_json::from_str::<serde_json::Value>(&content)
        .map_err(|error| format!("failed to parse `{}`: {error}", path.display()))
}

fn verify_virus_summary_tsv(path: &Path, label: &'static str, report: &mut AuditVerifyReport) {
    let mut reader = match ReaderBuilder::new().delimiter(b'\t').from_path(path) {
        Ok(reader) => reader,
        Err(error) => {
            report.fail(
                label,
                format!("failed to open `{}`: {error}", path.display()),
            );
            return;
        }
    };

    let headers = match reader.headers() {
        Ok(headers) => headers.clone(),
        Err(error) => {
            report.fail(
                label,
                format!(
                    "failed to parse TSV header in `{}`: {error}",
                    path.display()
                ),
            );
            return;
        }
    };
    let expected_headers = VIRUS_SUMMARY_HEADER.into_iter().collect::<Vec<_>>();
    if headers.iter().collect::<Vec<_>>() != expected_headers {
        report.fail(
            label,
            format!("unexpected TSV header in `{}`", path.display()),
        );
        return;
    }

    match reader.records().collect::<std::result::Result<Vec<_>, _>>() {
        Ok(rows) => report.pass(label, format!("parsed {} virus summary row(s)", rows.len())),
        Err(error) => report.fail(
            label,
            format!("failed to parse TSV rows in `{}`: {error}", path.display()),
        ),
    }
}

fn verify_sample_run_summary_facts(
    sample_run_summary: &serde_json::Value,
    summary: &SampleSummary,
    rounds: Option<&[RoundRecord]>,
    report: &mut AuditVerifyReport,
) {
    let mut mismatches = Vec::new();
    push_json_string_mismatch(
        &mut mismatches,
        sample_run_summary,
        "sample_id",
        &summary.sample_id,
    );
    push_json_string_mismatch(
        &mut mismatches,
        sample_run_summary,
        "decision_status",
        decision_status_label(&summary.decision_status),
    );
    push_json_string_mismatch(
        &mut mismatches,
        sample_run_summary,
        "release_status",
        release_status_label(&summary.release_status),
    );
    push_json_string_mismatch(
        &mut mismatches,
        sample_run_summary,
        "stop_reason",
        stop_reason_label(&summary.stop_reason),
    );
    push_json_u64_mismatch(
        &mut mismatches,
        sample_run_summary,
        "input_fragments",
        summary.input_fragments,
    );
    push_json_u64_mismatch(
        &mut mismatches,
        sample_run_summary,
        "qc_passing_fragments",
        summary.qc_passing_fragments,
    );
    push_json_u64_mismatch(
        &mut mismatches,
        sample_run_summary,
        "sampled_fragments",
        summary.sampled_fragments,
    );
    push_json_u64_mismatch(
        &mut mismatches,
        sample_run_summary,
        "rounds_run",
        rounds.map_or(summary.rounds_run, |rounds| rounds.len() as u64),
    );

    if mismatches.is_empty() {
        report.pass(
            "canonical sample summary facts",
            "summary/sample_run_summary.json agrees with sample_summary.json",
        );
    } else {
        report.fail(
            "canonical sample summary facts",
            format!("summary facts mismatch: {}", mismatches.join("; ")),
        );
    }
}

fn verify_result_overview_facts(
    result_overview: &serde_json::Value,
    summary: &SampleSummary,
    manifest: &RunManifest,
    report: &mut AuditVerifyReport,
) {
    let mut mismatches = Vec::new();
    push_nested_json_string_mismatch(
        &mut mismatches,
        result_overview,
        &["sample", "sample_id"],
        &summary.sample_id,
    );
    push_nested_json_string_mismatch(
        &mut mismatches,
        result_overview,
        &["sample", "decision_status"],
        decision_status_label(&summary.decision_status),
    );
    push_nested_json_string_mismatch(
        &mut mismatches,
        result_overview,
        &["sample", "release_status"],
        release_status_label(&summary.release_status),
    );
    push_nested_json_string_mismatch(
        &mut mismatches,
        result_overview,
        &["run", "reference_bundle_version"],
        &manifest.reference_bundle_version,
    );
    push_nested_json_string_mismatch(
        &mut mismatches,
        result_overview,
        &["run", "calibration_profile_version"],
        &manifest.calibration_profile_version,
    );
    push_nested_json_string_mismatch(
        &mut mismatches,
        result_overview,
        &["run", "backend"],
        &manifest.backend,
    );
    push_nested_json_string_mismatch(
        &mut mismatches,
        result_overview,
        &["run", "sampling_mode"],
        &manifest.sampling_mode,
    );
    push_nested_json_string_mismatch(
        &mut mismatches,
        result_overview,
        &["run", "negative_control_status"],
        negative_control_status_label(&manifest.negative_control.status),
    );

    if mismatches.is_empty() {
        report.pass(
            "canonical result overview facts",
            "summary/result_overview.json agrees with sample_summary.json and run_manifest.json",
        );
    } else {
        report.fail(
            "canonical result overview facts",
            format!("summary facts mismatch: {}", mismatches.join("; ")),
        );
    }
}

fn push_json_string_mismatch(
    mismatches: &mut Vec<String>,
    value: &serde_json::Value,
    key: &str,
    expected: &str,
) {
    let observed = value.get(key).and_then(serde_json::Value::as_str);
    if observed != Some(expected) {
        mismatches.push(format!(
            "{key} expected `{expected}`, observed `{}`",
            observed.unwrap_or("<missing or non-string>")
        ));
    }
}

fn push_json_u64_mismatch(
    mismatches: &mut Vec<String>,
    value: &serde_json::Value,
    key: &str,
    expected: u64,
) {
    let observed = value.get(key).and_then(serde_json::Value::as_u64);
    if observed != Some(expected) {
        mismatches.push(format!(
            "{key} expected `{expected}`, observed `{}`",
            observed
                .map(|value| value.to_string())
                .unwrap_or_else(|| "<missing or non-u64>".to_string())
        ));
    }
}

fn push_nested_json_string_mismatch(
    mismatches: &mut Vec<String>,
    value: &serde_json::Value,
    path: &[&str],
    expected: &str,
) {
    let observed = path
        .iter()
        .try_fold(value, |current, key| current.get(*key))
        .and_then(serde_json::Value::as_str);
    if observed != Some(expected) {
        mismatches.push(format!(
            "{} expected `{expected}`, observed `{}`",
            path.join("."),
            observed.unwrap_or("<missing or non-string>")
        ));
    }
}

fn parse_json_file<T>(
    path: &Path,
    label: &'static str,
    report: &mut AuditVerifyReport,
    success_detail: impl FnOnce(&T) -> String,
) -> Option<T>
where
    T: for<'de> Deserialize<'de>,
{
    match fs::read_to_string(path) {
        Ok(content) => match serde_json::from_str::<T>(&content) {
            Ok(value) => {
                report.pass(label, success_detail(&value));
                Some(value)
            }
            Err(error) => {
                report.fail(
                    label,
                    format!("failed to parse `{}`: {error}", path.display()),
                );
                None
            }
        },
        Err(error) => {
            report.fail(
                label,
                format!("failed to read `{}`: {error}", path.display()),
            );
            None
        }
    }
}

fn parse_rounds_file(path: &Path, report: &mut AuditVerifyReport) -> Option<Vec<RoundRecord>> {
    parse_tsv_file::<RoundRecord>(path, "rounds.tsv", report)
}

fn parse_candidate_calls_file(
    path: &Path,
    report: &mut AuditVerifyReport,
) -> Option<Vec<CandidateCiRecord>> {
    let mut reader = match ReaderBuilder::new().delimiter(b'\t').from_path(path) {
        Ok(reader) => reader,
        Err(error) => {
            report.fail(
                "candidate_calls.tsv",
                format!("failed to open `{}`: {error}", path.display()),
            );
            return None;
        }
    };

    let rows = match reader
        .deserialize::<CandidateCallTsvRow>()
        .collect::<std::result::Result<Vec<_>, _>>()
    {
        Ok(rows) => rows,
        Err(error) => {
            report.fail(
                "candidate_calls.tsv",
                format!("failed to parse TSV rows in `{}`: {error}", path.display()),
            );
            return None;
        }
    };

    let mut candidates = Vec::with_capacity(rows.len());
    for row in rows {
        let fraction_ci_95 = match serde_json::from_str::<[f64; 2]>(&row.fraction_ci_95) {
            Ok(value) => value,
            Err(error) => {
                report.fail(
                    "candidate_calls.tsv",
                    format!(
                        "candidate `{}` has invalid fraction_ci_95 JSON: {error}",
                        row.accession_or_group
                    ),
                );
                return None;
            }
        };
        let decision_reasons = match serde_json::from_str::<Vec<String>>(&row.decision_reasons) {
            Ok(value) => value,
            Err(error) => {
                report.fail(
                    "candidate_calls.tsv",
                    format!(
                        "candidate `{}` has invalid decision_reasons JSON: {error}",
                        row.accession_or_group
                    ),
                );
                return None;
            }
        };

        candidates.push(CandidateCiRecord {
            accession_or_group: row.accession_or_group,
            fraction_ci_95,
            decision_reasons,
        });
    }

    report.pass(
        "candidate_calls.tsv",
        format!("parsed {} candidate row(s)", candidates.len()),
    );
    Some(candidates)
}

fn parse_tsv_file<T>(
    path: &Path,
    label: &'static str,
    report: &mut AuditVerifyReport,
) -> Option<Vec<T>>
where
    T: for<'de> Deserialize<'de>,
{
    let mut reader = match ReaderBuilder::new().delimiter(b'\t').from_path(path) {
        Ok(reader) => reader,
        Err(error) => {
            report.fail(
                label,
                format!("failed to open `{}`: {error}", path.display()),
            );
            return None;
        }
    };

    match reader
        .deserialize::<T>()
        .collect::<std::result::Result<Vec<_>, _>>()
    {
        Ok(rows) => Some(rows),
        Err(error) => {
            report.fail(
                label,
                format!("failed to parse `{}`: {error}", path.display()),
            );
            None
        }
    }
}

fn verify_internal_bindings(
    summary: &SampleSummary,
    manifest: &RunManifest,
    report: &mut AuditVerifyReport,
) {
    let mut mismatches = Vec::new();

    if summary.reference_bundle_version != manifest.reference_bundle_version {
        mismatches.push(format!(
            "reference bundle version mismatch: sample_summary=`{}` run_manifest=`{}`",
            summary.reference_bundle_version, manifest.reference_bundle_version
        ));
    }
    if summary.calibration_profile_version != manifest.calibration_profile_version {
        mismatches.push(format!(
            "calibration profile version mismatch: sample_summary=`{}` run_manifest=`{}`",
            summary.calibration_profile_version, manifest.calibration_profile_version
        ));
    }
    if summary.backend != manifest.backend {
        mismatches.push(format!(
            "backend mismatch: sample_summary=`{}` run_manifest=`{}`",
            summary.backend, manifest.backend
        ));
    }
    if summary.seed != manifest.seed {
        mismatches.push(format!(
            "seed mismatch: sample_summary=`{}` run_manifest=`{}`",
            summary.seed, manifest.seed
        ));
    }

    if mismatches.is_empty() {
        report.pass(
            "report bundle bindings",
            "sample_summary.json and run_manifest.json agree on reference bundle, calibration profile, backend, and seed",
        );
    } else {
        report.fail("report bundle bindings", mismatches.join("; "));
    }
}

fn verify_candidate_ci_labels(candidates: &[CandidateCiRecord], report: &mut AuditVerifyReport) {
    let offenders = candidates
        .iter()
        .filter_map(|candidate| {
            if candidate.fraction_ci_95[0] > candidate.fraction_ci_95[1] {
                return Some(format!(
                    "candidate `{}` has invalid fraction_ci_95 bounds [{}, {}]",
                    candidate.accession_or_group,
                    candidate.fraction_ci_95[0],
                    candidate.fraction_ci_95[1]
                ));
            }
            sampling_only_ci_label_violation(&candidate.decision_reasons).map(|violation| {
                format!("candidate `{}` {violation}", candidate.accession_or_group)
            })
        })
        .collect::<Vec<_>>();

    if offenders.is_empty() {
        report.pass(
            "candidate CI labeling",
            format!(
                "{} candidate row(s) preserve the required sampling-only CI marker",
                candidates.len()
            ),
        );
    } else {
        report.fail("candidate CI labeling", offenders.join("; "));
    }
}

fn verify_reference_bundle_binding(
    reference_bundle_dir: &Path,
    summary: Option<&SampleSummary>,
    manifest: Option<&RunManifest>,
    report: &mut AuditVerifyReport,
) {
    let loaded_bundle = match load_reference_bundle(reference_bundle_dir) {
        Ok(bundle) => bundle,
        Err(error) => {
            report.fail(
                "reference bundle binding",
                format!(
                    "failed to load reference bundle `{}`: {error}",
                    reference_bundle_dir.display()
                ),
            );
            return;
        }
    };

    let Some(summary) = summary else {
        report.skip(
            "reference bundle binding",
            "skipped version comparison because sample_summary.json did not parse",
        );
        return;
    };
    let Some(manifest) = manifest else {
        report.skip(
            "reference bundle binding",
            "skipped version comparison because run_manifest.json did not parse",
        );
        return;
    };

    let expected = &loaded_bundle.bundle.version;
    let mut mismatches = Vec::new();
    if &summary.reference_bundle_version != expected {
        mismatches.push(format!(
            "sample_summary reference bundle version `{}` does not match loaded bundle `{expected}`",
            summary.reference_bundle_version
        ));
    }
    if &manifest.reference_bundle_version != expected {
        mismatches.push(format!(
            "run_manifest reference bundle version `{}` does not match loaded bundle `{expected}`",
            manifest.reference_bundle_version
        ));
    }

    if mismatches.is_empty() {
        report.pass(
            "reference bundle binding",
            format!("report references loaded reference bundle version `{expected}`"),
        );
    } else {
        report.fail("reference bundle binding", mismatches.join("; "));
    }
}

fn verify_calibration_profile_binding(
    calibration_profile_dir: &Path,
    summary: Option<&SampleSummary>,
    manifest: Option<&RunManifest>,
    report: &mut AuditVerifyReport,
) -> Option<ProfileToml> {
    let loaded_profile = match load_profile(calibration_profile_dir) {
        Ok(profile) => profile,
        Err(error) => {
            report.fail(
                "calibration profile binding",
                format!(
                    "failed to load calibration profile `{}`: {error}",
                    calibration_profile_dir.display()
                ),
            );
            return None;
        }
    };

    let Some(summary) = summary else {
        report.skip(
            "calibration profile binding",
            "skipped version comparison because sample_summary.json did not parse",
        );
        return Some(loaded_profile.profile);
    };
    let Some(manifest) = manifest else {
        report.skip(
            "calibration profile binding",
            "skipped version comparison because run_manifest.json did not parse",
        );
        return Some(loaded_profile.profile);
    };

    let mut mismatches = Vec::new();
    if summary.calibration_profile_version != loaded_profile.profile.profile_id {
        mismatches.push(format!(
            "sample_summary calibration profile version `{}` does not match loaded profile `{}`",
            summary.calibration_profile_version, loaded_profile.profile.profile_id
        ));
    }
    if manifest.calibration_profile_version != loaded_profile.profile.profile_id {
        mismatches.push(format!(
            "run_manifest calibration profile version `{}` does not match loaded profile `{}`",
            manifest.calibration_profile_version, loaded_profile.profile.profile_id
        ));
    }
    if summary.reference_bundle_version != loaded_profile.profile.reference_bundle {
        mismatches.push(format!(
            "sample_summary reference bundle version `{}` does not match calibration profile binding `{}`",
            summary.reference_bundle_version, loaded_profile.profile.reference_bundle
        ));
    }
    if manifest.reference_bundle_version != loaded_profile.profile.reference_bundle {
        mismatches.push(format!(
            "run_manifest reference bundle version `{}` does not match calibration profile binding `{}`",
            manifest.reference_bundle_version, loaded_profile.profile.reference_bundle
        ));
    }
    if summary.backend != loaded_profile.profile.backend {
        mismatches.push(format!(
            "sample_summary backend `{}` does not match calibration profile backend `{}`",
            summary.backend, loaded_profile.profile.backend
        ));
    }
    if manifest.backend != loaded_profile.profile.backend {
        mismatches.push(format!(
            "run_manifest backend `{}` does not match calibration profile backend `{}`",
            manifest.backend, loaded_profile.profile.backend
        ));
    }

    if mismatches.is_empty() {
        report.pass(
            "calibration profile binding",
            format!(
                "report references loaded calibration profile `{}` bound to reference bundle `{}`",
                loaded_profile.profile.profile_id, loaded_profile.profile.reference_bundle
            ),
        );
    } else {
        report.fail("calibration profile binding", mismatches.join("; "));
    }

    Some(loaded_profile.profile)
}

fn verify_release_status(
    summary: Option<&SampleSummary>,
    manifest: Option<&RunManifest>,
    profile: Option<&ProfileToml>,
    calibration_profile_dir: Option<&PathBuf>,
    report: &mut AuditVerifyReport,
) {
    let Some(summary) = summary else {
        report.skip(
            "release_status compliance",
            "skipped because sample_summary.json did not parse",
        );
        return;
    };
    let Some(profile) = profile else {
        report.skip(
            "release_status compliance",
            "skipped because --calibration-profile was not provided or did not load",
        );
        return;
    };

    let Some(manifest) = manifest else {
        report.skip(
            "release_status compliance",
            "skipped because run_manifest.json did not parse",
        );
        return;
    };

    let Some(calibration_profile_dir) = calibration_profile_dir else {
        report.skip(
            "release_status compliance",
            "skipped because --calibration-profile was not provided or did not load",
        );
        return;
    };

    if manifest.sampling_mode != profile.sampling.mode {
        report.fail(
            "release_status compliance",
            format!(
                "run_manifest sampling_mode `{}` does not match calibration profile sampling mode `{}`",
                manifest.sampling_mode, profile.sampling.mode
            ),
        );
        return;
    }

    if manifest.negative_control.required != profile.negative_control_required {
        report.fail(
            "release_status compliance",
            format!(
                "run_manifest negative_control.required `{}` does not match calibration profile negative_control_required `{}`",
                manifest.negative_control.required, profile.negative_control_required
            ),
        );
        return;
    }

    let benchmark_gates = match crate::decision::load_benchmark_gates(calibration_profile_dir) {
        Ok(gates) => gates,
        Err(error) => {
            report.fail(
                "release_status compliance",
                format!("failed to load calibration release_gate.json: {error}"),
            );
            return;
        }
    };

    let negative_control = negative_control_from_manifest(&manifest.negative_control);
    let expected = match expected_release_status(
        &manifest.sampling_mode,
        manifest.negative_control.required,
        &negative_control,
        benchmark_gates.as_ref(),
        &profile.status,
    ) {
        Ok(status) => status,
        Err(error) => {
            report.fail("release_status compliance", error.to_string());
            return;
        }
    };

    if summary.release_status != expected {
        report.fail(
            "release_status compliance",
            format!(
                "sample_summary release_status `{}` does not match recomputed release_status `{}` for sampling_mode `{}` / negative_control_status `{}` / negative_control_required={} / profile_status `{}`",
                release_status_label(&summary.release_status),
                release_status_label(&expected),
                manifest.sampling_mode,
                negative_control_status_label(&manifest.negative_control.status),
                manifest.negative_control.required,
                profile.status,
            ),
        );
    } else {
        report.pass(
            "release_status compliance",
            format!(
                "sample_summary release_status `{}` matches recomputed release gate for sampling_mode `{}` / negative_control_status `{}` / negative_control_required={} / profile_status `{}`",
                release_status_label(&summary.release_status),
                manifest.sampling_mode,
                negative_control_status_label(&manifest.negative_control.status),
                manifest.negative_control.required,
                profile.status,
            ),
        );
    }
}

fn negative_control_from_manifest(
    manifest: &crate::types::NegativeControlManifest,
) -> NegativeControlDecisionInput {
    match manifest.status {
        NegativeControlStatus::Missing => NegativeControlDecisionInput::Missing,
        NegativeControlStatus::Pass => NegativeControlDecisionInput::Passed {
            comparator: crate::decision::BackgroundComparator::new(&NegativeControlResult {
                control_id: manifest
                    .control_id
                    .clone()
                    .unwrap_or_else(|| "unknown".to_string()),
                control_status: manifest
                    .control_status
                    .clone()
                    .unwrap_or_else(|| "pass".to_string()),
                candidates: Vec::new(),
            }),
            result: NegativeControlResult {
                control_id: manifest
                    .control_id
                    .clone()
                    .unwrap_or_else(|| "unknown".to_string()),
                control_status: manifest
                    .control_status
                    .clone()
                    .unwrap_or_else(|| "pass".to_string()),
                candidates: Vec::new(),
            },
        },
        NegativeControlStatus::Fail => {
            NegativeControlDecisionInput::Failed(NegativeControlResult {
                control_id: manifest
                    .control_id
                    .clone()
                    .unwrap_or_else(|| "unknown".to_string()),
                control_status: manifest
                    .control_status
                    .clone()
                    .unwrap_or_else(|| "fail".to_string()),
                candidates: Vec::new(),
            })
        }
    }
}

fn negative_control_status_label(status: &NegativeControlStatus) -> &'static str {
    match status {
        NegativeControlStatus::Missing => "missing",
        NegativeControlStatus::Pass => "pass",
        NegativeControlStatus::Fail => "fail",
    }
}

fn decision_status_label(status: &DecisionStatus) -> &'static str {
    match status {
        DecisionStatus::Positive => "positive",
        DecisionStatus::Negative => "negative",
        DecisionStatus::Indeterminate => "indeterminate",
    }
}

fn release_status_label(release_status: &ReleaseStatus) -> &'static str {
    match release_status {
        ReleaseStatus::Final => "final",
        ReleaseStatus::Provisional => "provisional",
        ReleaseStatus::Blocked => "blocked",
    }
}

fn stop_reason_label(stop_reason: &StopReason) -> &'static str {
    match stop_reason {
        StopReason::PositiveBoundaryCrossed => "positive_boundary_crossed",
        StopReason::NegativeBoundaryConfirmed => "negative_boundary_confirmed",
        StopReason::MaxRoundsReached => "max_rounds_reached",
    }
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::decision::SAMPLING_ONLY_CI_LABEL;
    use crate::report::contract::FRACTION_CI_REASON_PREFIX;
    use crate::report::writer::{CANDIDATE_CALLS_HEADER, ROUNDS_HEADER};
    use crate::report::ReportWriter;
    use crate::types::{
        CandidateCall, DecisionStatus, EvidenceStrength, NegativeControlManifest,
        NegativeControlStatus, ProfileToml, SamplingConfig, StopReason,
    };
    use std::fs;
    use tempfile::tempdir;

    #[test]
    fn valid_report_bundle_passes() {
        let temp_dir = tempdir().expect("temp dir should be created");
        let fixture = write_audit_fixture(temp_dir.path(), "rvscreen_ref_2026.04.21-r1");

        let report = run_audit_verify(&AuditVerifyArgs {
            report_bundle: fixture.report_dir.clone(),
            reference_bundle: Some(fixture.reference_dir.clone()),
            calibration_profile: Some(fixture.profile_dir.clone()),
        })
        .expect("audit verify should run");

        assert!(report.passed(), "{}", report.render());
        assert!(report.render().contains("PASS rvscreen audit verify"));
        assert!(report.render().contains("candidate CI labeling"));
        assert!(report.render().contains("schema metadata contract"));
        assert!(report.render().contains("canonical sample summary facts"));
    }

    #[test]
    fn valid_legacy_only_report_bundle_passes() {
        let temp_dir = tempdir().expect("temp dir should be created");
        let fixture = write_audit_fixture(temp_dir.path(), "rvscreen_ref_2026.04.21-r1");
        make_legacy_only_bundle(&fixture.report_dir);

        let report = run_audit_verify(&AuditVerifyArgs {
            report_bundle: fixture.report_dir.clone(),
            reference_bundle: Some(fixture.reference_dir.clone()),
            calibration_profile: Some(fixture.profile_dir.clone()),
        })
        .expect("audit verify should run");

        assert!(report.passed(), "{}", report.render());
        assert!(report.render().contains("PASS rvscreen audit verify"));
        assert!(report.render().contains("legacy-only bundle"));
        assert!(report.render().contains("candidate CI labeling"));
    }

    #[test]
    fn candidate_tsv_mismatch_fails_with_specific_message() {
        let temp_dir = tempdir().expect("temp dir should be created");
        let fixture = write_audit_fixture(temp_dir.path(), "rvscreen_ref_2026.04.21-r1");
        let canonical_candidate_path = fixture
            .report_dir
            .join(CANONICAL_CANDIDATE_CALLS_ACCESSION_TSV);
        let candidate_content = fs::read_to_string(&canonical_candidate_path)
            .expect("canonical candidate calls should be readable")
            .replace("NC_SYNTHV1.1", "NC_SYNTHV1_DRIFT.1");
        fs::write(&canonical_candidate_path, candidate_content)
            .expect("canonical candidate calls should be rewritten");

        let report = run_audit_verify(&AuditVerifyArgs {
            report_bundle: fixture.report_dir,
            reference_bundle: Some(fixture.reference_dir),
            calibration_profile: Some(fixture.profile_dir),
        })
        .expect("audit verify should run");

        assert!(!report.passed(), "audit should fail");
        assert!(
            report.render().contains("candidate TSV mismatch"),
            "{}",
            report.render()
        );
    }

    #[test]
    fn run_manifest_mismatch_fails_with_specific_message() {
        let temp_dir = tempdir().expect("temp dir should be created");
        let fixture = write_audit_fixture(temp_dir.path(), "rvscreen_ref_2026.04.21-r1");
        let canonical_manifest_path = fixture.report_dir.join(CANONICAL_RUN_MANIFEST_JSON);
        let mut canonical_manifest: RunManifest = serde_json::from_str(
            &fs::read_to_string(&canonical_manifest_path)
                .expect("canonical run manifest should be readable"),
        )
        .expect("canonical run manifest should parse");
        canonical_manifest.seed += 1;
        fs::write(
            &canonical_manifest_path,
            serde_json::to_vec_pretty(&canonical_manifest)
                .expect("canonical run manifest should serialize"),
        )
        .expect("canonical run manifest should be rewritten");

        let report = run_audit_verify(&AuditVerifyArgs {
            report_bundle: fixture.report_dir,
            reference_bundle: Some(fixture.reference_dir),
            calibration_profile: Some(fixture.profile_dir),
        })
        .expect("audit verify should run");

        assert!(!report.passed(), "audit should fail");
        assert!(
            report.render().contains("run manifest mismatch"),
            "{}",
            report.render()
        );
    }

    #[test]
    fn missing_canonical_mandatory_artifact_fails_when_schema_says_v2() {
        let temp_dir = tempdir().expect("temp dir should be created");
        let fixture = write_audit_fixture(temp_dir.path(), "rvscreen_ref_2026.04.21-r1");
        fs::remove_file(fixture.report_dir.join(CANONICAL_SAMPLING_ROUNDS_TSV))
            .expect("canonical rounds.tsv should be removed");

        let report = run_audit_verify(&AuditVerifyArgs {
            report_bundle: fixture.report_dir,
            reference_bundle: Some(fixture.reference_dir),
            calibration_profile: Some(fixture.profile_dir),
        })
        .expect("audit verify should run");

        assert!(!report.passed(), "audit should fail");
        assert!(
            report
                .render()
                .contains("missing canonical mandatory artifact(s): results/sampling_rounds.tsv"),
            "{}",
            report.render()
        );
    }

    #[test]
    fn schema_metadata_mismatch_fails_with_specific_message() {
        let temp_dir = tempdir().expect("temp dir should be created");
        let fixture = write_audit_fixture(temp_dir.path(), "rvscreen_ref_2026.04.21-r1");
        let schema_path = fixture.report_dir.join(CANONICAL_SCHEMA_VERSIONS_JSON);
        let mut schema: serde_json::Value = serde_json::from_str(
            &fs::read_to_string(&schema_path).expect("schema versions should be readable"),
        )
        .expect("schema versions should parse");
        schema["canonical_outputs"]["run_manifest"] =
            serde_json::Value::String("provenance/run_manifest.drift.json".to_string());
        fs::write(
            &schema_path,
            serde_json::to_vec_pretty(&schema).expect("schema versions should serialize"),
        )
        .expect("schema versions should be rewritten");

        let report = run_audit_verify(&AuditVerifyArgs {
            report_bundle: fixture.report_dir,
            reference_bundle: Some(fixture.reference_dir),
            calibration_profile: Some(fixture.profile_dir),
        })
        .expect("audit verify should run");

        assert!(!report.passed(), "audit should fail");
        assert!(
            report
                .render()
                .contains("canonical_outputs does not match REPORT_ARTIFACT_CONTRACTS"),
            "{}",
            report.render()
        );
    }

    #[test]
    fn missing_required_file_fails_with_specific_message() {
        let temp_dir = tempdir().expect("temp dir should be created");
        let fixture = write_audit_fixture(temp_dir.path(), "rvscreen_ref_2026.04.21-r1");
        fs::remove_file(fixture.report_dir.join(ROUNDS_TSV)).expect("rounds.tsv should be removed");

        let report = run_audit_verify(&AuditVerifyArgs {
            report_bundle: fixture.report_dir,
            reference_bundle: Some(fixture.reference_dir),
            calibration_profile: Some(fixture.profile_dir),
        })
        .expect("audit verify should run");

        assert!(!report.passed(), "audit should fail");
        assert!(
            report.render().contains("rounds.tsv"),
            "{}",
            report.render()
        );
        assert!(
            report
                .render()
                .contains("missing required file(s): rounds.tsv"),
            "{}",
            report.render()
        );
    }

    #[test]
    fn version_mismatch_fails_with_specific_message() {
        let temp_dir = tempdir().expect("temp dir should be created");
        let fixture = write_audit_fixture(temp_dir.path(), "rvscreen_ref_2026.04.21-r1");
        write_reference_bundle_dir(&fixture.reference_dir, "rvscreen_ref_2099.01.01-r1");

        let report = run_audit_verify(&AuditVerifyArgs {
            report_bundle: fixture.report_dir,
            reference_bundle: Some(fixture.reference_dir),
            calibration_profile: Some(fixture.profile_dir),
        })
        .expect("audit verify should run");

        assert!(!report.passed(), "audit should fail");
        assert!(
            report.render().contains("reference bundle version"),
            "{}",
            report.render()
        );
        assert!(
            report
                .render()
                .contains("does not match loaded bundle `rvscreen_ref_2099.01.01-r1`"),
            "{}",
            report.render()
        );
    }

    #[test]
    fn duplicate_fraction_ci_labels_fail_audit() {
        let temp_dir = tempdir().expect("temp dir should be created");
        let fixture = write_audit_fixture(temp_dir.path(), "rvscreen_ref_2026.04.21-r1");

        let mut candidate_rows = candidate_fixtures();
        candidate_rows[0].decision_reasons.push(format!(
            "fraction_ci_95_label={}",
            crate::decision::SAMPLING_ONLY_CI_LABEL
        ));
        ReportWriter::write(
            &fixture.report_dir,
            &sample_summary_fixture("rvscreen_ref_2026.04.21-r1"),
            &candidate_rows,
            &round_fixtures(),
            &run_manifest_fixture("rvscreen_ref_2026.04.21-r1"),
        )
        .expect_err("writer should reject duplicate CI markers before audit");
    }

    #[test]
    fn forbidden_vocabulary_is_blocked_in_generated_artifacts_only() {
        let checked_vocab = checked_report_bundle_vocabulary();
        assert!(
            forbidden_terms_in_vocab(&checked_vocab).is_empty(),
            "forbidden boundary terms should not appear in generated report vocabulary"
        );

        let mut injected_vocab = checked_vocab.clone();
        injected_vocab.push("clinical_positive".to_string());
        let hits = forbidden_terms_in_vocab(&injected_vocab);
        assert!(
            hits.contains(&"clinical_positive".to_string()),
            "positive control should detect an injected forbidden term"
        );

        for allowed in [
            "decision",
            "evidence_strength",
            "decision_reasons",
            "entity_id",
            "aggregation_level",
            "rank",
            "display_priority",
        ] {
            assert!(
                checked_vocab.iter().any(|term| term == allowed),
                "neutral term `{allowed}` should remain allowed"
            );
        }
    }

    #[test]
    fn streaming_final_release_status_fails_audit() {
        let temp_dir = tempdir().expect("temp dir should be created");
        let fixture = write_audit_fixture(temp_dir.path(), "rvscreen_ref_2026.04.21-r1");
        write_profile_dir_with_mode_and_negative_control(
            &fixture.profile_dir,
            "rvscreen_ref_2026.04.21-r1",
            "streaming",
            false,
        );
        write_run_manifest_with_negative_control(&fixture.report_dir, |manifest| {
            manifest.sampling_mode = "streaming".to_string();
        });
        write_sample_summary_with_release_status(&fixture.report_dir, ReleaseStatus::Final);

        let report = run_audit_verify(&AuditVerifyArgs {
            report_bundle: fixture.report_dir,
            reference_bundle: Some(fixture.reference_dir),
            calibration_profile: Some(fixture.profile_dir),
        })
        .expect("audit verify should run");

        assert!(!report.passed(), "audit should fail");
        assert!(
            report
                .render()
                .contains("sample_summary release_status `final` does not match recomputed release_status `provisional`")
                && report.render().contains("sampling_mode `streaming`"),
            "{}",
            report.render()
        );
    }

    #[test]
    fn representative_required_negative_control_cannot_be_provisional() {
        let temp_dir = tempdir().expect("temp dir should be created");
        let fixture = write_audit_fixture(temp_dir.path(), "rvscreen_ref_2026.04.21-r1");
        write_profile_dir_with_mode_and_negative_control(
            &fixture.profile_dir,
            "rvscreen_ref_2026.04.21-r1",
            "representative",
            true,
        );
        write_run_manifest_with_negative_control(&fixture.report_dir, |manifest| {
            manifest.negative_control.required = true;
        });
        write_sample_summary_with_release_status(&fixture.report_dir, ReleaseStatus::Provisional);

        let report = run_audit_verify(&AuditVerifyArgs {
            report_bundle: fixture.report_dir,
            reference_bundle: Some(fixture.reference_dir),
            calibration_profile: Some(fixture.profile_dir),
        })
        .expect("audit verify should run");

        assert!(!report.passed(), "audit should fail");
        assert!(
            report
                .render()
                .contains("does not match recomputed release_status `final`"),
            "{}",
            report.render()
        );
    }

    #[test]
    fn missing_checksum_file_fails_audit() {
        let temp_dir = tempdir().expect("temp dir should be created");
        let fixture = write_audit_fixture(temp_dir.path(), "rvscreen_ref_2026.04.21-r1");
        fs::remove_file(fixture.report_dir.join(CHECKSUM_SHA256))
            .expect("checksum.sha256 should be removed");

        let report = run_audit_verify(&AuditVerifyArgs {
            report_bundle: fixture.report_dir,
            reference_bundle: Some(fixture.reference_dir),
            calibration_profile: Some(fixture.profile_dir),
        })
        .expect("audit verify should run");

        assert!(!report.passed(), "audit should fail");
        assert!(
            report
                .render()
                .contains("missing required file(s): checksum.sha256"),
            "{}",
            report.render()
        );
    }

    #[test]
    fn representative_required_negative_control_missing_is_blocked() {
        let temp_dir = tempdir().expect("temp dir should be created");
        let fixture = write_audit_fixture(temp_dir.path(), "rvscreen_ref_2026.04.21-r1");
        write_profile_dir_with_mode_and_negative_control(
            &fixture.profile_dir,
            "rvscreen_ref_2026.04.21-r1",
            "representative",
            true,
        );
        write_run_manifest_with_negative_control(&fixture.report_dir, |manifest| {
            manifest.negative_control.required = true;
            manifest.negative_control.status = NegativeControlStatus::Missing;
            manifest.negative_control.control_id = None;
            manifest.negative_control.control_status = None;
        });
        write_sample_summary_with_release_status(&fixture.report_dir, ReleaseStatus::Final);

        let report = run_audit_verify(&AuditVerifyArgs {
            report_bundle: fixture.report_dir,
            reference_bundle: Some(fixture.reference_dir),
            calibration_profile: Some(fixture.profile_dir),
        })
        .expect("audit verify should run");

        assert!(!report.passed(), "audit should fail");
        assert!(
            report
                .render()
                .contains("does not match recomputed release_status `blocked`"),
            "{}",
            report.render()
        );
    }

    fn make_legacy_only_bundle(report_dir: &Path) {
        remove_dir_if_exists(&report_dir.join("summary"));
        remove_dir_if_exists(&report_dir.join("results"));
        remove_dir_if_exists(&report_dir.join("evidence"));
        remove_dir_if_exists(&report_dir.join("provenance"));
        remove_dir_if_exists(&report_dir.join("audit"));
        remove_file_if_exists(&report_dir.join("logs").join("run.ndjson"));
        write_legacy_checksum(report_dir);
    }

    fn remove_dir_if_exists(path: &Path) {
        if path.exists() {
            fs::remove_dir_all(path).expect("directory should be removed");
        }
    }

    fn remove_file_if_exists(path: &Path) {
        if path.exists() {
            fs::remove_file(path).expect("file should be removed");
        }
    }

    fn write_legacy_checksum(report_dir: &Path) {
        let mut checksum_body = String::new();
        for relative_path in [
            SAMPLE_SUMMARY_JSON,
            CANDIDATE_CALLS_TSV,
            RUN_MANIFEST_JSON,
            ROUNDS_TSV,
        ] {
            let digest = sha256_file(&report_dir.join(relative_path))
                .expect("legacy checksum input should hash");
            checksum_body.push_str(&format!("{digest}  {relative_path}\n"));
        }
        fs::write(report_dir.join(CHECKSUM_SHA256), checksum_body)
            .expect("legacy checksum should be rewritten");
    }

    fn checked_report_bundle_vocabulary() -> Vec<String> {
        let schema_versions = serde_json::json!({
            "report_bundle_schema": "rvscreen.report_bundle.v2",
            "canonical_outputs": REPORT_ARTIFACT_CONTRACTS
                .iter()
                .map(|contract| (contract.key, contract.canonical_path))
                .collect::<std::collections::BTreeMap<_, _>>(),
            "legacy_outputs": REPORT_ARTIFACT_CONTRACTS
                .iter()
                .filter_map(|contract| contract
                    .legacy_alias_path
                    .map(|legacy_alias_path| (contract.key, legacy_alias_path)))
                .collect::<std::collections::BTreeMap<_, _>>(),
        });
        let release_gate = serde_json::json!({
            "schema_version": "rvscreen.release_gate.v1",
            "sample_id": "sample-1",
            "release_status": "provisional",
            "decision_status": "positive",
            "stop_reason": "positive_boundary_crossed",
            "rounds_run": 3,
            "reference_bundle_version": "rvscreen_ref_2026.04.21-r1",
            "calibration_profile_version": "rvscreen_ref_2026.04.21-r1",
            "backend": "minimap2",
            "seed": 1,
            "negative_control": {"status": "pass"},
            "input_file_count": 1,
        });
        let sample_run_summary = serde_json::json!({
            "schema_version": "rvscreen.sample_run_summary.v1",
            "sample_id": "sample-1",
            "decision_status": "positive",
            "release_status": "provisional",
            "stop_reason": "positive_boundary_crossed",
            "input_fragments": 10,
            "qc_passing_fragments": 10,
            "sampled_fragments": 10,
            "rounds_run": 3,
            "positive_entity_count": 1,
            "indeterminate_entity_count": 0,
            "negative_entity_count": 0,
        });
        let result_overview = serde_json::json!({
            "schema_version": "rvscreen.result_overview.v1",
            "sample": {
                "sample_id": "sample-1",
                "decision_status": "positive",
                "release_status": "provisional",
            },
            "run": {
                "reference_bundle_version": "rvscreen_ref_2026.04.21-r1",
                "calibration_profile_version": "rvscreen_ref_2026.04.21-r1",
                "backend": "minimap2",
                "sampling_mode": "representative",
                "negative_control_status": "pass",
            },
            "top_results": [{
                "rank": 1,
                "entity_id": "ebv",
                "entity_name": "Epstein-Barr virus",
                "aggregation_level": "accession",
                "decision": "positive",
                "evidence_strength": "high",
                "accepted_fragments": 10,
                "unique_fraction": 0.5,
                "breadth": 0.25,
                "background_ratio": 0.01,
            }],
        });

        let mut vocab = Vec::new();
        vocab.extend(collect_json_keys(&schema_versions));
        vocab.extend(collect_json_keys(&release_gate));
        vocab.extend(collect_json_keys(&sample_run_summary));
        vocab.extend(collect_json_keys(&result_overview));
        vocab.extend(
            REPORT_ARTIFACT_CONTRACTS
                .iter()
                .flat_map(|contract| [contract.key, contract.canonical_path])
                .map(|term| term.to_string()),
        );
        vocab.extend(VIRUS_SUMMARY_HEADER.iter().map(|term| term.to_string()));
        vocab.extend(CANDIDATE_CALLS_HEADER.iter().copied().map(String::from));
        vocab.extend(ROUNDS_HEADER.iter().copied().map(String::from));
        vocab.extend(
            [
                "decision",
                "evidence_strength",
                "decision_reasons",
                "entity_id",
                "aggregation_level",
                "rank",
                "display_priority",
            ]
            .into_iter()
            .map(String::from),
        );
        vocab
    }

    fn collect_json_keys(value: &serde_json::Value) -> Vec<String> {
        let mut keys = Vec::new();
        collect_json_keys_recursive(value, &mut keys);
        keys
    }

    fn collect_json_keys_recursive(value: &serde_json::Value, keys: &mut Vec<String>) {
        match value {
            serde_json::Value::Object(map) => {
                for (key, nested) in map {
                    keys.push(key.clone());
                    collect_json_keys_recursive(nested, keys);
                }
            }
            serde_json::Value::Array(items) => {
                for item in items {
                    collect_json_keys_recursive(item, keys);
                }
            }
            _ => {}
        }
    }

    fn forbidden_terms_in_vocab(vocab: &[String]) -> Vec<String> {
        const FORBIDDEN: &[&str] = &[
            "clinical_positive",
            "tumor_related",
            "report_conclusion",
            "cancer_matched",
            "should_report",
            "interpretation_text",
            "pathogenic",
            "diagnostic",
            "clinically_significant",
            "cancer",
            "tumor",
        ];

        vocab
            .iter()
            .filter(|term| FORBIDDEN.iter().any(|forbidden| term.contains(forbidden)))
            .cloned()
            .collect()
    }

    struct AuditFixture {
        report_dir: PathBuf,
        reference_dir: PathBuf,
        profile_dir: PathBuf,
    }

    fn write_audit_fixture(base_dir: &Path, reference_version: &str) -> AuditFixture {
        let report_dir = base_dir.join("report-bundle");
        let reference_dir = base_dir.join("reference-bundle");
        let profile_dir = base_dir.join("calibration-profile");
        let summary = sample_summary_fixture(reference_version);
        let manifest = run_manifest_fixture(reference_version);
        let rounds = round_fixtures();
        let candidates = candidate_fixtures();

        ReportWriter::write(&report_dir, &summary, &candidates, &rounds, &manifest)
            .expect("report bundle should be written");
        write_reference_bundle_dir(&reference_dir, reference_version);
        write_profile_dir(&profile_dir, reference_version);

        AuditFixture {
            report_dir,
            reference_dir,
            profile_dir,
        }
    }

    fn write_reference_bundle_dir(dir: &Path, version: &str) {
        fs::create_dir_all(dir).expect("reference bundle dir should exist");
        fs::write(
            dir.join("bundle.toml"),
            toml::to_string_pretty(&crate::types::BundleToml {
                version: version.to_string(),
                created_at: "2026-04-21T00:00:00Z".to_string(),
                included_layers: vec!["host_backbone".to_string(), "viral_panel".to_string()],
            })
            .expect("bundle.toml should serialize"),
        )
        .expect("bundle.toml should be written");
    }

    fn write_profile_dir(dir: &Path, reference_version: &str) {
        write_profile_dir_with_mode_and_negative_control(
            dir,
            reference_version,
            "representative",
            false,
        );
    }

    fn write_profile_dir_with_mode_and_negative_control(
        dir: &Path,
        reference_version: &str,
        sampling_mode: &str,
        negative_control_required: bool,
    ) {
        fs::create_dir_all(dir).expect("profile dir should exist");
        let profile = ProfileToml {
            profile_id: "rvscreen_calib_2026.04.21-r1".to_string(),
            status: "release_candidate".to_string(),
            reference_bundle: reference_version.to_string(),
            backend: "minimap2".to_string(),
            preset: "sr-conservative".to_string(),
            seed: 20_260_421,
            supported_input: vec!["fastq".to_string(), "fastq.gz".to_string()],
            supported_read_type: vec!["illumina_pe_shortread".to_string()],
            negative_control_required,
            sampling: SamplingConfig {
                mode: sampling_mode.to_string(),
                rounds: vec![50_000, 100_000],
                round_mode: None,
                round_proportions: None,
                max_rounds: 2,
            },
            fragment_rules: crate::types::FragmentRules {
                min_mapq: 20,
                min_as_diff: 12,
                max_nm: 8,
                require_pair_consistency: true,
            },
            candidate_rules: crate::types::CandidateRules {
                min_nonoverlap_fragments: 1,
                min_breadth: 0.0,
                max_background_ratio: 0.0,
            },
            decision_rules: crate::types::DecisionRules {
                theta_pos: 0.001,
                theta_neg: 0.0,
                allow_indeterminate: true,
            },
        };
        fs::write(
            dir.join("profile.toml"),
            toml::to_string_pretty(&profile).expect("profile.toml should serialize"),
        )
        .expect("profile.toml should be written");
        fs::write(
            dir.join("release_gate.json"),
            serde_json::to_vec_pretty(&crate::calibration::ReleaseGate {
                backend_gate: crate::calibration::BackendGate {
                    status: crate::calibration::GateStatus::Pass,
                    details: "audit fixture backend gate".to_string(),
                },
                reference_gate: crate::calibration::ReferenceGate {
                    status: crate::calibration::GateStatus::Pass,
                    details: "audit fixture reference gate".to_string(),
                },
                specificity_gate: crate::calibration::SpecificityGate {
                    status: crate::calibration::GateStatus::Pass,
                    details: "audit fixture specificity gate".to_string(),
                    negative_samples: 1,
                    false_positives: 0,
                },
                sensitivity_gate: crate::calibration::SensitivityGate {
                    status: crate::calibration::GateStatus::Pass,
                    details: "audit fixture sensitivity gate".to_string(),
                    spike_in_detected: 1,
                    spike_in_total: 1,
                },
            })
            .expect("release_gate.json should serialize"),
        )
        .expect("release_gate.json should be written");
    }

    fn write_sample_summary_with_release_status(report_dir: &Path, release_status: ReleaseStatus) {
        let mut summary = sample_summary_fixture("rvscreen_ref_2026.04.21-r1");
        summary.release_status = release_status;
        fs::write(
            report_dir.join(SAMPLE_SUMMARY_JSON),
            serde_json::to_vec_pretty(&summary).expect("sample summary should serialize"),
        )
        .expect("sample_summary.json should be rewritten");
        refresh_report_checksum(report_dir);
    }

    fn write_run_manifest_with_negative_control(
        report_dir: &Path,
        mutate: impl FnOnce(&mut RunManifest),
    ) {
        let manifest_path = report_dir.join(RUN_MANIFEST_JSON);
        let mut manifest: RunManifest = serde_json::from_str(
            &fs::read_to_string(&manifest_path).expect("run_manifest.json should be readable"),
        )
        .expect("run_manifest.json should parse");
        mutate(&mut manifest);
        fs::write(
            &manifest_path,
            serde_json::to_vec_pretty(&manifest).expect("run manifest should serialize"),
        )
        .expect("run_manifest.json should be rewritten");
        refresh_report_checksum(report_dir);
    }

    fn refresh_report_checksum(report_dir: &Path) {
        let summary: SampleSummary = serde_json::from_str(
            &fs::read_to_string(report_dir.join(SAMPLE_SUMMARY_JSON))
                .expect("sample_summary.json should be readable"),
        )
        .expect("sample_summary.json should parse");
        let manifest: RunManifest = serde_json::from_str(
            &fs::read_to_string(report_dir.join(RUN_MANIFEST_JSON))
                .expect("run_manifest.json should be readable"),
        )
        .expect("run_manifest.json should parse");

        ReportWriter::write(
            report_dir,
            &summary,
            &candidate_fixtures(),
            &round_fixtures(),
            &manifest,
        )
        .expect("report bundle should be rewritten with refreshed checksum");
    }

    fn sample_summary_fixture(reference_version: &str) -> SampleSummary {
        SampleSummary {
            sample_id: "S-AUDIT-001".to_string(),
            reference_bundle_version: reference_version.to_string(),
            calibration_profile_version: "rvscreen_calib_2026.04.21-r1".to_string(),
            backend: "minimap2".to_string(),
            seed: 20_260_421,
            input_fragments: 100_000,
            qc_passing_fragments: 95_000,
            sampled_fragments: 100_000,
            rounds_run: 2,
            stop_reason: StopReason::PositiveBoundaryCrossed,
            decision_status: DecisionStatus::Positive,
            release_status: ReleaseStatus::Final,
        }
    }

    fn run_manifest_fixture(reference_version: &str) -> RunManifest {
        RunManifest {
            reference_bundle_version: reference_version.to_string(),
            calibration_profile_version: "rvscreen_calib_2026.04.21-r1".to_string(),
            backend: "minimap2".to_string(),
            seed: 20_260_421,
            sampling_mode: "representative".to_string(),
            sampling_round_plan: None,
            negative_control: NegativeControlManifest {
                required: false,
                status: NegativeControlStatus::Pass,
                control_id: Some("neg-001".to_string()),
                control_status: Some("pass".to_string()),
            },
            input_files: vec![
                "reads_R1.fastq.gz".to_string(),
                "reads_R2.fastq.gz".to_string(),
            ],
        }
    }

    fn round_fixtures() -> Vec<RoundRecord> {
        vec![
            RoundRecord {
                sampled_fragments: 50_000,
                accepted_virus: 4,
                decision_status: DecisionStatus::Indeterminate,
            },
            RoundRecord {
                sampled_fragments: 100_000,
                accepted_virus: 12,
                decision_status: DecisionStatus::Positive,
            },
        ]
    }

    fn candidate_fixtures() -> Vec<CandidateCall> {
        vec![CandidateCall {
            virus_name: "Synthetic audit virus".to_string(),
            taxid: 9_999,
            accession_or_group: "NC_SYNTHV1.1".to_string(),
            accepted_fragments: 12,
            nonoverlap_fragments: 4,
            raw_fraction: 0.00012,
            unique_fraction: 12.0 / 100_000.0,
            fraction_ci_95: [0.00008, 0.00018],
            clopper_pearson_upper: 0.0002,
            breadth: 0.01,
            ambiguous_fragments: 1,
            background_ratio: 0.0,
            decision: DecisionStatus::Positive,
            decision_reasons: vec![
                format!("{FRACTION_CI_REASON_PREFIX}{SAMPLING_ONLY_CI_LABEL}"),
                "unique_fraction_above_theta_pos".to_string(),
            ],
            evidence_strength: EvidenceStrength::High,
        }]
    }
}
