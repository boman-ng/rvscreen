use crate::decision::SAMPLING_ONLY_CI_LABEL;

pub const SAMPLE_SUMMARY_JSON: &str = "sample_summary.json";
pub const CANDIDATE_CALLS_TSV: &str = "candidate_calls.tsv";
pub const RUN_MANIFEST_JSON: &str = "run_manifest.json";
pub const ROUNDS_TSV: &str = "rounds.tsv";
pub const CHECKSUM_SHA256: &str = "checksum.sha256";
pub const COVERAGE_DIR: &str = "coverage";
pub const LOGS_DIR: &str = "logs";
pub const RUN_LOG: &str = "run.log";
pub const COVERAGE_HEADER: &str = "position\tdepth\n";
pub const FRACTION_CI_REASON_PREFIX: &str = "fraction_ci_95_label=";

pub const SUMMARY_DIR: &str = "summary";
pub const RESULTS_DIR: &str = "results";
pub const EVIDENCE_DIR: &str = "evidence";
pub const PROVENANCE_DIR: &str = "provenance";
pub const AUDIT_DIR: &str = "audit";
pub const POSITIVE_CANDIDATE_COVERAGE_DIR: &str = "positive_candidate_coverage";

pub const CANONICAL_SAMPLE_RUN_SUMMARY_JSON: &str = "summary/sample_run_summary.json";
pub const CANONICAL_RESULT_OVERVIEW_JSON: &str = "summary/result_overview.json";
pub const CANONICAL_VIRUS_SUMMARY_TSV: &str = "summary/virus_summary.tsv";
pub const CANONICAL_CANDIDATE_CALLS_ACCESSION_TSV: &str = "results/candidate_calls.accession.tsv";
pub const CANONICAL_SAMPLING_ROUNDS_TSV: &str = "results/sampling_rounds.tsv";
pub const CANONICAL_POSITIVE_CANDIDATE_COVERAGE_DIR: &str = "evidence/positive_candidate_coverage";
pub const CANONICAL_RUN_MANIFEST_JSON: &str = "provenance/run_manifest.json";
pub const CANONICAL_CHECKSUM_SHA256: &str = "audit/checksum.sha256";
pub const CANONICAL_SCHEMA_VERSIONS_JSON: &str = "audit/schema_versions.json";
pub const CANONICAL_RELEASE_GATE_JSON: &str = "audit/release_gate.json";
pub const CANONICAL_RUN_NDJSON: &str = "logs/run.ndjson";

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ReportArtifactKind {
    File,
    Directory,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum MandatoryStatus {
    Mandatory,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ChecksumInclusion {
    Include,
    Exclude,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ReportArtifactContract {
    pub key: &'static str,
    pub canonical_path: &'static str,
    pub legacy_alias_path: Option<&'static str>,
    pub source_type_or_helper: &'static str,
    pub mandatory_status: MandatoryStatus,
    pub checksum_inclusion: ChecksumInclusion,
    pub verifier_expectation: &'static str,
    pub kind: ReportArtifactKind,
}

pub const REPORT_ARTIFACT_CONTRACTS: &[ReportArtifactContract] = &[
    ReportArtifactContract {
        key: "sample_run_summary",
        canonical_path: CANONICAL_SAMPLE_RUN_SUMMARY_JSON,
        legacy_alias_path: Some(SAMPLE_SUMMARY_JSON),
        source_type_or_helper: "SampleSummary JSON writer",
        mandatory_status: MandatoryStatus::Mandatory,
        checksum_inclusion: ChecksumInclusion::Include,
        verifier_expectation: "JSON object parses and agrees with provenance/run_manifest.json",
        kind: ReportArtifactKind::File,
    },
    ReportArtifactContract {
        key: "result_overview",
        canonical_path: CANONICAL_RESULT_OVERVIEW_JSON,
        legacy_alias_path: None,
        source_type_or_helper: "result overview JSON writer",
        mandatory_status: MandatoryStatus::Mandatory,
        checksum_inclusion: ChecksumInclusion::Include,
        verifier_expectation: "JSON object parses as the compact result overview",
        kind: ReportArtifactKind::File,
    },
    ReportArtifactContract {
        key: "virus_summary",
        canonical_path: CANONICAL_VIRUS_SUMMARY_TSV,
        legacy_alias_path: None,
        source_type_or_helper: "virus summary TSV writer",
        mandatory_status: MandatoryStatus::Mandatory,
        checksum_inclusion: ChecksumInclusion::Include,
        verifier_expectation: "TSV header and rows parse as the virus summary table",
        kind: ReportArtifactKind::File,
    },
    ReportArtifactContract {
        key: "candidate_calls_accession",
        canonical_path: CANONICAL_CANDIDATE_CALLS_ACCESSION_TSV,
        legacy_alias_path: Some(CANDIDATE_CALLS_TSV),
        source_type_or_helper: "candidate calls accession TSV writer",
        mandatory_status: MandatoryStatus::Mandatory,
        checksum_inclusion: ChecksumInclusion::Include,
        verifier_expectation: "TSV header and rows parse as accession candidate calls",
        kind: ReportArtifactKind::File,
    },
    ReportArtifactContract {
        key: "sampling_rounds",
        canonical_path: CANONICAL_SAMPLING_ROUNDS_TSV,
        legacy_alias_path: Some(ROUNDS_TSV),
        source_type_or_helper: "sampling rounds TSV writer",
        mandatory_status: MandatoryStatus::Mandatory,
        checksum_inclusion: ChecksumInclusion::Include,
        verifier_expectation: "TSV header and rows parse as sampling rounds",
        kind: ReportArtifactKind::File,
    },
    ReportArtifactContract {
        key: "positive_candidate_coverage",
        canonical_path: CANONICAL_POSITIVE_CANDIDATE_COVERAGE_DIR,
        legacy_alias_path: Some(COVERAGE_DIR),
        source_type_or_helper: "positive candidate coverage TSV writer",
        mandatory_status: MandatoryStatus::Mandatory,
        checksum_inclusion: ChecksumInclusion::Include,
        verifier_expectation:
            "Directory exists and contains parseable coverage TSV files when calls are present",
        kind: ReportArtifactKind::Directory,
    },
    ReportArtifactContract {
        key: "run_manifest",
        canonical_path: CANONICAL_RUN_MANIFEST_JSON,
        legacy_alias_path: Some(RUN_MANIFEST_JSON),
        source_type_or_helper: "RunManifest JSON writer",
        mandatory_status: MandatoryStatus::Mandatory,
        checksum_inclusion: ChecksumInclusion::Include,
        verifier_expectation: "JSON object parses and agrees with summary/sample_run_summary.json",
        kind: ReportArtifactKind::File,
    },
    ReportArtifactContract {
        key: "checksum_sha256",
        canonical_path: CANONICAL_CHECKSUM_SHA256,
        legacy_alias_path: Some(CHECKSUM_SHA256),
        source_type_or_helper: "checksum writer",
        mandatory_status: MandatoryStatus::Mandatory,
        checksum_inclusion: ChecksumInclusion::Exclude,
        verifier_expectation: "Checksum file exists and is excluded from its own checksum set",
        kind: ReportArtifactKind::File,
    },
    ReportArtifactContract {
        key: "schema_versions",
        canonical_path: CANONICAL_SCHEMA_VERSIONS_JSON,
        legacy_alias_path: None,
        source_type_or_helper: "schema versions JSON writer",
        mandatory_status: MandatoryStatus::Mandatory,
        checksum_inclusion: ChecksumInclusion::Include,
        verifier_expectation: "JSON object parses and lists bundle schema versions",
        kind: ReportArtifactKind::File,
    },
    ReportArtifactContract {
        key: "release_gate",
        canonical_path: CANONICAL_RELEASE_GATE_JSON,
        legacy_alias_path: None,
        source_type_or_helper: "release gate JSON writer",
        mandatory_status: MandatoryStatus::Mandatory,
        checksum_inclusion: ChecksumInclusion::Include,
        verifier_expectation: "JSON object parses and records release gate status",
        kind: ReportArtifactKind::File,
    },
    ReportArtifactContract {
        key: "run_log",
        canonical_path: CANONICAL_RUN_NDJSON,
        legacy_alias_path: Some("logs/run.log"),
        source_type_or_helper: "structured run log writer",
        mandatory_status: MandatoryStatus::Mandatory,
        checksum_inclusion: ChecksumInclusion::Include,
        verifier_expectation: "NDJSON lines parse as structured run events",
        kind: ReportArtifactKind::File,
    },
];

const _: [&str; 6] = [
    SUMMARY_DIR,
    RESULTS_DIR,
    EVIDENCE_DIR,
    PROVENANCE_DIR,
    AUDIT_DIR,
    POSITIVE_CANDIDATE_COVERAGE_DIR,
];
const _: &[ReportArtifactContract] = REPORT_ARTIFACT_CONTRACTS;

pub const REQUIRED_REPORT_BUNDLE_FILES: [&str; 5] = [
    SAMPLE_SUMMARY_JSON,
    CANDIDATE_CALLS_TSV,
    RUN_MANIFEST_JSON,
    ROUNDS_TSV,
    CHECKSUM_SHA256,
];

pub fn sampling_only_ci_label(reasons: &[String]) -> Option<&str> {
    let mut labels = fraction_ci_labels(reasons);
    let first = labels.next()?;
    if labels.next().is_none() && first == SAMPLING_ONLY_CI_LABEL {
        Some(first)
    } else {
        None
    }
}

pub fn sampling_only_ci_label_violation(reasons: &[String]) -> Option<String> {
    let labels = fraction_ci_labels(reasons).collect::<Vec<_>>();

    match labels.as_slice() {
        [] => Some(format!(
            "is missing the required `{FRACTION_CI_REASON_PREFIX}{SAMPLING_ONLY_CI_LABEL}` marker"
        )),
        [label] if *label == SAMPLING_ONLY_CI_LABEL => None,
        [label] => Some(format!("carries non-sampling-only label `{label}`")),
        _ => Some(format!(
            "must carry exactly one `{FRACTION_CI_REASON_PREFIX}...` marker, found [{}]",
            labels.join(", ")
        )),
    }
}

fn fraction_ci_labels<'a>(reasons: &'a [String]) -> impl Iterator<Item = &'a str> + 'a {
    reasons
        .iter()
        .filter_map(|reason| reason.strip_prefix(FRACTION_CI_REASON_PREFIX))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn accepts_exactly_one_sampling_only_ci_label() {
        let reasons = vec![
            format!("{FRACTION_CI_REASON_PREFIX}{SAMPLING_ONLY_CI_LABEL}"),
            "unique_fraction_above_theta_pos".to_string(),
        ];

        assert_eq!(
            sampling_only_ci_label(&reasons),
            Some(SAMPLING_ONLY_CI_LABEL)
        );
        assert_eq!(sampling_only_ci_label_violation(&reasons), None);
    }

    #[test]
    fn rejects_duplicate_ci_labels() {
        let reasons = vec![
            format!("{FRACTION_CI_REASON_PREFIX}{SAMPLING_ONLY_CI_LABEL}"),
            format!("{FRACTION_CI_REASON_PREFIX}{SAMPLING_ONLY_CI_LABEL}"),
        ];

        assert_eq!(sampling_only_ci_label(&reasons), None);
        assert!(sampling_only_ci_label_violation(&reasons)
            .expect("duplicate labels should be rejected")
            .contains("exactly one"),);
    }

    #[test]
    fn canonical_constants_match_phase_one_paths() {
        assert_eq!(SUMMARY_DIR, "summary");
        assert_eq!(RESULTS_DIR, "results");
        assert_eq!(EVIDENCE_DIR, "evidence");
        assert_eq!(PROVENANCE_DIR, "provenance");
        assert_eq!(AUDIT_DIR, "audit");
        assert_eq!(
            POSITIVE_CANDIDATE_COVERAGE_DIR,
            "positive_candidate_coverage"
        );
        assert_eq!(
            CANONICAL_SAMPLE_RUN_SUMMARY_JSON,
            "summary/sample_run_summary.json"
        );
        assert_eq!(
            CANONICAL_RESULT_OVERVIEW_JSON,
            "summary/result_overview.json"
        );
        assert_eq!(CANONICAL_VIRUS_SUMMARY_TSV, "summary/virus_summary.tsv");
        assert_eq!(
            CANONICAL_CANDIDATE_CALLS_ACCESSION_TSV,
            "results/candidate_calls.accession.tsv"
        );
        assert_eq!(CANONICAL_SAMPLING_ROUNDS_TSV, "results/sampling_rounds.tsv");
        assert_eq!(
            CANONICAL_POSITIVE_CANDIDATE_COVERAGE_DIR,
            "evidence/positive_candidate_coverage"
        );
        assert_eq!(CANONICAL_RUN_MANIFEST_JSON, "provenance/run_manifest.json");
        assert_eq!(CANONICAL_CHECKSUM_SHA256, "audit/checksum.sha256");
        assert_eq!(CANONICAL_SCHEMA_VERSIONS_JSON, "audit/schema_versions.json");
        assert_eq!(CANONICAL_RELEASE_GATE_JSON, "audit/release_gate.json");
        assert_eq!(CANONICAL_RUN_NDJSON, "logs/run.ndjson");
    }

    #[test]
    fn legacy_constants_remain_compatible() {
        assert_eq!(SAMPLE_SUMMARY_JSON, "sample_summary.json");
        assert_eq!(CANDIDATE_CALLS_TSV, "candidate_calls.tsv");
        assert_eq!(RUN_MANIFEST_JSON, "run_manifest.json");
        assert_eq!(ROUNDS_TSV, "rounds.tsv");
        assert_eq!(CHECKSUM_SHA256, "checksum.sha256");
        assert_eq!(COVERAGE_DIR, "coverage");
        assert_eq!(LOGS_DIR, "logs");
        assert_eq!(RUN_LOG, "run.log");
        assert_eq!(COVERAGE_HEADER, "position\tdepth\n");
        assert_eq!(
            REQUIRED_REPORT_BUNDLE_FILES,
            [
                SAMPLE_SUMMARY_JSON,
                CANDIDATE_CALLS_TSV,
                RUN_MANIFEST_JSON,
                ROUNDS_TSV,
                CHECKSUM_SHA256,
            ]
        );
    }

    #[test]
    fn contract_matrix_covers_phase_one_artifacts() {
        let expected = [
            (
                "sample_run_summary",
                CANONICAL_SAMPLE_RUN_SUMMARY_JSON,
                Some(SAMPLE_SUMMARY_JSON),
                ChecksumInclusion::Include,
                ReportArtifactKind::File,
            ),
            (
                "result_overview",
                CANONICAL_RESULT_OVERVIEW_JSON,
                None,
                ChecksumInclusion::Include,
                ReportArtifactKind::File,
            ),
            (
                "virus_summary",
                CANONICAL_VIRUS_SUMMARY_TSV,
                None,
                ChecksumInclusion::Include,
                ReportArtifactKind::File,
            ),
            (
                "candidate_calls_accession",
                CANONICAL_CANDIDATE_CALLS_ACCESSION_TSV,
                Some(CANDIDATE_CALLS_TSV),
                ChecksumInclusion::Include,
                ReportArtifactKind::File,
            ),
            (
                "sampling_rounds",
                CANONICAL_SAMPLING_ROUNDS_TSV,
                Some(ROUNDS_TSV),
                ChecksumInclusion::Include,
                ReportArtifactKind::File,
            ),
            (
                "positive_candidate_coverage",
                CANONICAL_POSITIVE_CANDIDATE_COVERAGE_DIR,
                Some(COVERAGE_DIR),
                ChecksumInclusion::Include,
                ReportArtifactKind::Directory,
            ),
            (
                "run_manifest",
                CANONICAL_RUN_MANIFEST_JSON,
                Some(RUN_MANIFEST_JSON),
                ChecksumInclusion::Include,
                ReportArtifactKind::File,
            ),
            (
                "checksum_sha256",
                CANONICAL_CHECKSUM_SHA256,
                Some(CHECKSUM_SHA256),
                ChecksumInclusion::Exclude,
                ReportArtifactKind::File,
            ),
            (
                "schema_versions",
                CANONICAL_SCHEMA_VERSIONS_JSON,
                None,
                ChecksumInclusion::Include,
                ReportArtifactKind::File,
            ),
            (
                "release_gate",
                CANONICAL_RELEASE_GATE_JSON,
                None,
                ChecksumInclusion::Include,
                ReportArtifactKind::File,
            ),
            (
                "run_log",
                CANONICAL_RUN_NDJSON,
                Some("logs/run.log"),
                ChecksumInclusion::Include,
                ReportArtifactKind::File,
            ),
        ];

        assert_eq!(REPORT_ARTIFACT_CONTRACTS.len(), expected.len());

        for (contract, (key, canonical_path, legacy_alias_path, checksum_inclusion, kind)) in
            REPORT_ARTIFACT_CONTRACTS.iter().zip(expected)
        {
            assert_eq!(contract.key, key);
            assert_eq!(contract.canonical_path, canonical_path);
            assert_eq!(contract.legacy_alias_path, legacy_alias_path);
            assert_eq!(contract.mandatory_status, MandatoryStatus::Mandatory);
            assert_eq!(contract.checksum_inclusion, checksum_inclusion);
            assert_eq!(contract.kind, kind);
            assert!(!contract.source_type_or_helper.is_empty());
            assert!(!contract.verifier_expectation.is_empty());
        }

        let checksum_contract = REPORT_ARTIFACT_CONTRACTS
            .iter()
            .find(|contract| contract.key == "checksum_sha256")
            .expect("checksum contract should be present");
        assert_eq!(
            checksum_contract.checksum_inclusion,
            ChecksumInclusion::Exclude
        );
    }
}
