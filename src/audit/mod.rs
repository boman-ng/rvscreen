use crate::calibration::{load_profile, load_reference_bundle};
use crate::cli::AuditVerifyArgs;
use crate::decision::SAMPLING_ONLY_CI_LABEL;
use crate::error::{Result, RvScreenError};
use crate::types::{ProfileToml, ReleaseStatus, RoundRecord, RunManifest, SampleSummary};
use csv::ReaderBuilder;
use serde::Deserialize;
use sha2::{Digest, Sha256};
use std::fs::{self, File};
use std::io::{BufReader, Read};
use std::path::{Component, Path, PathBuf};

const SAMPLE_SUMMARY_JSON: &str = "sample_summary.json";
const CANDIDATE_CALLS_TSV: &str = "candidate_calls.tsv";
const RUN_MANIFEST_JSON: &str = "run_manifest.json";
const ROUNDS_TSV: &str = "rounds.tsv";
const CHECKSUM_SHA256: &str = "checksum.sha256";
const FRACTION_CI_REASON_PREFIX: &str = "fraction_ci_95_label=";

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

    let required_files = [
        SAMPLE_SUMMARY_JSON,
        CANDIDATE_CALLS_TSV,
        RUN_MANIFEST_JSON,
        ROUNDS_TSV,
    ];
    let missing_files = required_files
        .iter()
        .filter(|name| !args.report_bundle.join(name).is_file())
        .copied()
        .collect::<Vec<_>>();
    if missing_files.is_empty() {
        report.pass(
            "report bundle completeness",
            format!("found required files: {}", required_files.join(", ")),
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

    verify_release_status(summary.as_ref(), loaded_profile.as_ref(), &mut report);

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
            let observed = candidate
                .decision_reasons
                .iter()
                .find_map(|reason| reason.strip_prefix(FRACTION_CI_REASON_PREFIX));
            match observed {
                Some(label) if label == SAMPLING_ONLY_CI_LABEL => None,
                Some(label) => Some(format!(
                    "candidate `{}` carries non-sampling-only label `{label}`",
                    candidate.accession_or_group
                )),
                None => Some(format!(
                    "candidate `{}` is missing `{FRACTION_CI_REASON_PREFIX}{SAMPLING_ONLY_CI_LABEL}`",
                    candidate.accession_or_group
                )),
            }
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
    profile: Option<&ProfileToml>,
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

    let release_status = &summary.release_status;
    let sampling_mode = profile.sampling.mode.as_str();
    let mut violations = Vec::new();

    match sampling_mode {
        "representative" => {
            if profile.negative_control_required && *release_status == ReleaseStatus::Provisional {
                violations.push(
                    "release_status `provisional` is impossible for representative mode when negative_control_required=true under current task rules"
                        .to_string(),
                );
            }
        }
        "streaming" => {
            if *release_status == ReleaseStatus::Final {
                violations.push(
                    "release_status `final` is not allowed for streaming mode under current task rules"
                        .to_string(),
                );
            }
        }
        other => violations.push(format!(
            "unsupported calibration sampling mode `{other}` in loaded profile"
        )),
    }

    if violations.is_empty() {
        report.pass(
            "release_status compliance",
            format!(
                "no inferable release_status violations for mode `{}` / negative_control_required={} / release_status `{}`",
                profile.sampling.mode,
                profile.negative_control_required,
                release_status_label(release_status)
            ),
        );
    } else {
        report.fail("release_status compliance", violations.join("; "));
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
    use crate::report::ReportWriter;
    use crate::types::{
        CandidateCall, DecisionStatus, EvidenceStrength, ProfileToml, SamplingConfig, StopReason,
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
            negative_control_required: false,
            sampling: SamplingConfig {
                mode: "representative".to_string(),
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
            sampled_fragments: 50_000,
            rounds_run: 2,
            stop_reason: StopReason::PositiveBoundaryCrossed,
            decision_status: DecisionStatus::Positive,
            release_status: ReleaseStatus::Provisional,
        }
    }

    fn run_manifest_fixture(reference_version: &str) -> RunManifest {
        RunManifest {
            reference_bundle_version: reference_version.to_string(),
            calibration_profile_version: "rvscreen_calib_2026.04.21-r1".to_string(),
            backend: "minimap2".to_string(),
            seed: 20_260_421,
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
            unique_fraction: 0.00013,
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
