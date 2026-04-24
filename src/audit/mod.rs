use crate::calibration::{load_profile, load_reference_bundle};
use crate::cli::AuditVerifyArgs;
use crate::decision::{
    expected_release_status, NegativeControlDecisionInput, NegativeControlResult,
};
use crate::error::{Result, RvScreenError};
use crate::report::contract::{
    sampling_only_ci_label_violation, CANDIDATE_CALLS_TSV, CHECKSUM_SHA256,
    REQUIRED_REPORT_BUNDLE_FILES, ROUNDS_TSV, RUN_MANIFEST_JSON, SAMPLE_SUMMARY_JSON,
};
use crate::types::{
    NegativeControlStatus, ProfileToml, ReleaseStatus, RoundRecord, RunManifest, SampleSummary,
};
use csv::ReaderBuilder;
use serde::Deserialize;
use sha2::{Digest, Sha256};
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

pub fn run_audit_verify(args: &AuditVerifyArgs) -> Result<AuditVerifyReport> {
    let mut report = AuditVerifyReport::new(args.report_bundle.clone());

    let missing_files = REQUIRED_REPORT_BUNDLE_FILES
        .iter()
        .filter(|name| !args.report_bundle.join(name).is_file())
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

    verify_checksum_if_present(&args.report_bundle, &mut report);

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
        &args.report_bundle.join(RUN_MANIFEST_JSON),
        "run_manifest.json",
        &mut report,
        |manifest| {
            format!(
                "parsed backend `{}` bound to reference bundle `{}`",
                manifest.backend, manifest.reference_bundle_version
            )
        },
    );
    let rounds = parse_rounds_file(&args.report_bundle.join(ROUNDS_TSV), &mut report);
    let candidates =
        parse_candidate_calls_file(&args.report_bundle.join(CANDIDATE_CALLS_TSV), &mut report);

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

fn verify_checksum_if_present(report_bundle: &Path, report: &mut AuditVerifyReport) {
    let checksum_path = report_bundle.join(CHECKSUM_SHA256);
    if !checksum_path.exists() {
        report.skip(
            "checksum.sha256",
            "checksum.sha256 not present; checksum validation skipped",
        );
        return;
    }

    let result = verify_checksum_file(&checksum_path, report_bundle);
    match result {
        Ok(verified_files) => report.pass(
            "checksum.sha256",
            format!(
                "verified {verified_files} checksum entr{}",
                if verified_files == 1 { "y" } else { "ies" }
            ),
        ),
        Err(error) => report.fail("checksum.sha256", error.to_string()),
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

fn release_status_label(release_status: &ReleaseStatus) -> &'static str {
    match release_status {
        ReleaseStatus::Final => "final",
        ReleaseStatus::Provisional => "provisional",
        ReleaseStatus::Blocked => "blocked",
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
