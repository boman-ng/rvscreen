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

        assert_eq!(sampling_only_ci_label(&reasons), Some(SAMPLING_ONLY_CI_LABEL));
        assert_eq!(sampling_only_ci_label_violation(&reasons), None);
    }

    #[test]
    fn rejects_duplicate_ci_labels() {
        let reasons = vec![
            format!("{FRACTION_CI_REASON_PREFIX}{SAMPLING_ONLY_CI_LABEL}"),
            format!("{FRACTION_CI_REASON_PREFIX}{SAMPLING_ONLY_CI_LABEL}"),
        ];

        assert_eq!(sampling_only_ci_label(&reasons), None);
        assert!(
            sampling_only_ci_label_violation(&reasons)
                .expect("duplicate labels should be rejected")
                .contains("exactly one"),
        );
    }
}
