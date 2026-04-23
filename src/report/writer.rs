use crate::error::{Result, RvScreenError};
use crate::report::contract::{
    sampling_only_ci_label, sampling_only_ci_label_violation, CANDIDATE_CALLS_TSV, CHECKSUM_SHA256,
    COVERAGE_DIR, COVERAGE_HEADER, LOGS_DIR, ROUNDS_TSV, RUN_LOG, RUN_MANIFEST_JSON,
    SAMPLE_SUMMARY_JSON,
};
use crate::types::{CandidateCall, DecisionStatus, RoundRecord, RunManifest, SampleSummary};
use csv::WriterBuilder;
use serde::Serialize;
use serde_json::json;
use std::collections::BTreeSet;
use std::fs::{self, File};
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::Path;
use sha2::{Digest, Sha256};
use std::time::{SystemTime, UNIX_EPOCH};

#[derive(Debug, Default, Clone, Copy)]
pub struct ReportWriter;

impl ReportWriter {
    pub fn write(
        output_dir: impl AsRef<Path>,
        summary: &SampleSummary,
        candidates: &[CandidateCall],
        rounds: &[RoundRecord],
        manifest: &RunManifest,
    ) -> Result<()> {
        let output_dir = output_dir.as_ref();
        validate_bundle_inputs(summary, candidates, rounds, manifest)?;

        create_dir_all(output_dir)?;
        create_dir_all(&output_dir.join(COVERAGE_DIR))?;
        create_dir_all(&output_dir.join(LOGS_DIR))?;

        write_json_pretty(&output_dir.join(SAMPLE_SUMMARY_JSON), summary)?;
        write_candidate_calls_tsv(&output_dir.join(CANDIDATE_CALLS_TSV), candidates)?;
        write_json_pretty(&output_dir.join(RUN_MANIFEST_JSON), manifest)?;
        write_rounds_tsv(&output_dir.join(ROUNDS_TSV), rounds)?;
        write_coverage_files(&output_dir.join(COVERAGE_DIR), candidates)?;
        write_run_log(
            &output_dir.join(LOGS_DIR).join(RUN_LOG),
            output_dir,
            summary,
            candidates,
        )?;
        write_checksum_file(
            &output_dir.join(CHECKSUM_SHA256),
            output_dir,
            &collect_report_bundle_files(output_dir)?,
        )?;

        Ok(())
    }
}

fn validate_bundle_inputs(
    summary: &SampleSummary,
    candidates: &[CandidateCall],
    rounds: &[RoundRecord],
    manifest: &RunManifest,
) -> Result<()> {
    if summary.reference_bundle_version != manifest.reference_bundle_version {
        return Err(RvScreenError::validation(
            "report.reference_bundle_version",
            format!(
                "sample summary `{}` does not match run manifest `{}`",
                summary.reference_bundle_version, manifest.reference_bundle_version
            ),
        ));
    }

    if summary.calibration_profile_version != manifest.calibration_profile_version {
        return Err(RvScreenError::validation(
            "report.calibration_profile_version",
            format!(
                "sample summary `{}` does not match run manifest `{}`",
                summary.calibration_profile_version, manifest.calibration_profile_version
            ),
        ));
    }

    if summary.backend != manifest.backend {
        return Err(RvScreenError::validation(
            "report.backend",
            format!(
                "sample summary backend `{}` does not match run manifest `{}`",
                summary.backend, manifest.backend
            ),
        ));
    }

    if summary.seed != manifest.seed {
        return Err(RvScreenError::validation(
            "report.seed",
            format!(
                "sample summary seed `{}` does not match run manifest `{}`",
                summary.seed, manifest.seed
            ),
        ));
    }

    if summary.rounds_run != rounds.len() as u64 {
        return Err(RvScreenError::validation(
            "report.rounds_run",
            format!(
                "sample summary rounds_run `{}` does not match rounds.tsv record count `{}`",
                summary.rounds_run,
                rounds.len()
            ),
        ));
    }

    for candidate in candidates {
        if let Some(violation) = sampling_only_ci_label_violation(&candidate.decision_reasons) {
            return Err(RvScreenError::validation(
                "report.candidate_calls.fraction_ci_95",
                format!(
                    "candidate `{}` {violation}",
                    candidate.accession_or_group,
                ),
            ));
        }
    }

    Ok(())
}

fn write_json_pretty(path: &Path, value: &impl Serialize) -> Result<()> {
    let bytes = serde_json::to_vec_pretty(value).map_err(|error| {
        RvScreenError::validation(
            format!("report.serialize.{}", path.display()),
            error.to_string(),
        )
    })?;
    fs::write(path, bytes).map_err(|error| RvScreenError::io(path, error))
}

fn write_candidate_calls_tsv(path: &Path, candidates: &[CandidateCall]) -> Result<()> {
    #[derive(Serialize)]
    struct CandidateCallRow<'a> {
        virus_name: &'a str,
        taxid: u64,
        accession_or_group: &'a str,
        accepted_fragments: u64,
        nonoverlap_fragments: u64,
        raw_fraction: f64,
        unique_fraction: f64,
        fraction_ci_95: String,
        clopper_pearson_upper: f64,
        breadth: f64,
        ambiguous_fragments: u64,
        background_ratio: f64,
        decision: String,
        decision_reasons: String,
        evidence_strength: String,
    }

    let mut writer = WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .map_err(|error| {
            RvScreenError::validation(
                "report.candidate_calls.tsv",
                format!(
                    "failed to open TSV writer for `{}`: {error}",
                    path.display()
                ),
            )
        })?;

    for candidate in candidates {
        let row = CandidateCallRow {
            virus_name: &candidate.virus_name,
            taxid: candidate.taxid,
            accession_or_group: &candidate.accession_or_group,
            accepted_fragments: candidate.accepted_fragments,
            nonoverlap_fragments: candidate.nonoverlap_fragments,
            raw_fraction: candidate.raw_fraction,
            unique_fraction: candidate.unique_fraction,
            fraction_ci_95: json_string(&candidate.fraction_ci_95)?,
            clopper_pearson_upper: candidate.clopper_pearson_upper,
            breadth: candidate.breadth,
            ambiguous_fragments: candidate.ambiguous_fragments,
            background_ratio: candidate.background_ratio,
            decision: enum_label(&candidate.decision)?,
            decision_reasons: json_string(&candidate.decision_reasons)?,
            evidence_strength: enum_label(&candidate.evidence_strength)?,
        };

        writer.serialize(row).map_err(|error| {
            RvScreenError::validation(
                "report.candidate_calls.tsv",
                format!(
                    "failed to serialize candidate `{}`: {error}",
                    candidate.accession_or_group
                ),
            )
        })?;
    }

    writer
        .flush()
        .map_err(|error| RvScreenError::io(path, error))
}

fn write_rounds_tsv(path: &Path, rounds: &[RoundRecord]) -> Result<()> {
    #[derive(Serialize)]
    struct RoundRow {
        sampled_fragments: u64,
        accepted_virus: u64,
        decision_status: String,
    }

    let mut writer = WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .map_err(|error| {
            RvScreenError::validation(
                "report.rounds.tsv",
                format!(
                    "failed to open TSV writer for `{}`: {error}",
                    path.display()
                ),
            )
        })?;

    for round in rounds {
        writer
            .serialize(RoundRow {
                sampled_fragments: round.sampled_fragments,
                accepted_virus: round.accepted_virus,
                decision_status: enum_label(&round.decision_status)?,
            })
            .map_err(|error| {
                RvScreenError::validation(
                    "report.rounds.tsv",
                    format!("failed to serialize round: {error}"),
                )
            })?;
    }

    writer
        .flush()
        .map_err(|error| RvScreenError::io(path, error))
}

fn write_coverage_files(coverage_dir: &Path, candidates: &[CandidateCall]) -> Result<()> {
    let mut seen_paths = BTreeSet::new();

    for candidate in candidates
        .iter()
        .filter(|candidate| candidate.decision == DecisionStatus::Positive)
    {
        let file_name = format!(
            "{}.coverage.tsv",
            sanitize_file_component(&candidate.accession_or_group)
        );
        let coverage_path = coverage_dir.join(file_name);
        if !seen_paths.insert(coverage_path.clone()) {
            return Err(RvScreenError::validation(
                "report.coverage",
                format!(
                    "duplicate positive candidate coverage path for `{}`",
                    candidate.accession_or_group
                ),
            ));
        }

        fs::write(&coverage_path, COVERAGE_HEADER)
            .map_err(|error| RvScreenError::io(&coverage_path, error))?;
    }

    Ok(())
}

fn write_run_log(
    path: &Path,
    output_dir: &Path,
    summary: &SampleSummary,
    candidates: &[CandidateCall],
) -> Result<()> {
    let positive_candidates = candidates
        .iter()
        .filter(|candidate| candidate.decision == DecisionStatus::Positive)
        .map(|candidate| candidate.accession_or_group.as_str())
        .collect::<Vec<_>>();
    let ci_labels = candidates
        .iter()
        .filter_map(candidate_sampling_only_ci_label)
        .collect::<Vec<_>>();

    let started = json!({
        "timestamp_unix_ms": unix_timestamp_millis(),
        "level": "info",
        "event": "report_bundle_write_started",
        "sample_id": summary.sample_id,
        "output_dir": output_dir.display().to_string(),
        "candidate_count": candidates.len(),
        "rounds_run": summary.rounds_run,
    });
    let completed = json!({
        "timestamp_unix_ms": unix_timestamp_millis(),
        "level": "info",
        "event": "report_bundle_write_completed",
        "decision_status": enum_label(&summary.decision_status)?,
        "release_status": enum_label(&summary.release_status)?,
        "positive_candidates": positive_candidates,
        "fraction_ci_95_labels": ci_labels,
    });

    let body = format!("{}\n{}\n", json_string(&started)?, json_string(&completed)?);
    fs::write(path, body).map_err(|error| RvScreenError::io(path, error))
}

fn create_dir_all(path: &Path) -> Result<()> {
    fs::create_dir_all(path).map_err(|error| RvScreenError::io(path, error))
}

fn collect_report_bundle_files(output_dir: &Path) -> Result<Vec<std::path::PathBuf>> {
    let mut files = Vec::new();
    collect_report_bundle_files_recursive(output_dir, &mut files)?;
    files.retain(|path| path.file_name().and_then(|name| name.to_str()) != Some(CHECKSUM_SHA256));
    files.sort();
    Ok(files)
}

fn collect_report_bundle_files_recursive(
    current: &Path,
    files: &mut Vec<std::path::PathBuf>,
) -> Result<()> {
    let mut entries = fs::read_dir(current)
        .map_err(|error| RvScreenError::io(current, error))?
        .collect::<std::result::Result<Vec<_>, _>>()
        .map_err(|error| RvScreenError::io(current, error))?;
    entries.sort_by_key(|entry| entry.path());

    for entry in entries {
        let path = entry.path();
        if path.is_dir() {
            collect_report_bundle_files_recursive(&path, files)?;
        } else {
            files.push(path);
        }
    }

    Ok(())
}

fn write_checksum_file(checksum_path: &Path, bundle_dir: &Path, files: &[std::path::PathBuf]) -> Result<()> {
    let file = File::create(checksum_path).map_err(|error| RvScreenError::io(checksum_path, error))?;
    let mut writer = BufWriter::new(file);

    for path in files {
        let digest = sha256_file(path)?;
        let relative = path.strip_prefix(bundle_dir).map_err(|_| {
            RvScreenError::validation(
                "report.checksum",
                format!(
                    "`{}` is not inside `{}`",
                    path.display(),
                    bundle_dir.display()
                ),
            )
        })?;
        writeln!(writer, "{digest}  {}", relative.display())
            .map_err(|error| RvScreenError::io(checksum_path, error))?;
    }

    writer.flush().map_err(|error| RvScreenError::io(checksum_path, error))
}

fn sha256_file(path: &Path) -> Result<String> {
    let file = File::open(path).map_err(|error| RvScreenError::io(path, error))?;
    let mut reader = BufReader::new(file);
    let mut hasher = Sha256::new();
    let mut buffer = [0u8; 8192];

    loop {
        let read = reader
            .read(&mut buffer)
            .map_err(|error| RvScreenError::io(path, error))?;
        if read == 0 {
            break;
        }
        hasher.update(&buffer[..read]);
    }

    Ok(format!("{:x}", hasher.finalize()))
}

fn json_string(value: &impl Serialize) -> Result<String> {
    serde_json::to_string(value)
        .map_err(|error| RvScreenError::validation("report.serialize", error.to_string()))
}

fn enum_label(value: &impl Serialize) -> Result<String> {
    let serialized = serde_json::to_value(value)
        .map_err(|error| RvScreenError::validation("report.enum_label", error.to_string()))?;
    serialized.as_str().map(ToOwned::to_owned).ok_or_else(|| {
        RvScreenError::validation(
            "report.enum_label",
            "expected enum to serialize to a string label",
        )
    })
}

fn candidate_sampling_only_ci_label(candidate: &CandidateCall) -> Option<&str> {
    sampling_only_ci_label(&candidate.decision_reasons)
}

fn sanitize_file_component(value: &str) -> String {
    let sanitized = value
        .chars()
        .map(|character| match character {
            'a'..='z' | 'A'..='Z' | '0'..='9' | '.' | '_' | '-' => character,
            _ => '_',
        })
        .collect::<String>();

    if sanitized.is_empty() {
        "candidate".to_string()
    } else {
        sanitized
    }
}

fn unix_timestamp_millis() -> u128 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|duration| duration.as_millis())
        .unwrap_or(0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::decision::SAMPLING_ONLY_CI_LABEL;
    use crate::report::contract::FRACTION_CI_REASON_PREFIX;
    use crate::types::{ReleaseStatus, StopReason};
    use csv::ReaderBuilder;
    use std::collections::BTreeSet;
    use std::fs;
    use tempfile::TempDir;

    #[test]
    fn writes_report_bundle_files_and_content() {
        let temp_dir = TempDir::new().expect("temp dir should be created");
        let output_dir = temp_dir.path().join("report_bundle");
        let summary = sample_summary_fixture();
        let manifest = run_manifest_fixture();
        let candidates = candidate_fixtures();
        let rounds = round_fixtures();

        ReportWriter::write(&output_dir, &summary, &candidates, &rounds, &manifest)
            .expect("report bundle should be written");

        let files = collect_relative_paths(&output_dir);

        assert_eq!(
            files.iter().cloned().collect::<BTreeSet<_>>(),
            BTreeSet::from([
                "candidate_calls.tsv".to_string(),
                "checksum.sha256".to_string(),
                "coverage/ebv.coverage.tsv".to_string(),
                "logs/run.log".to_string(),
                "rounds.tsv".to_string(),
                "run_manifest.json".to_string(),
                "sample_summary.json".to_string(),
            ])
        );

        let written_summary: SampleSummary = serde_json::from_str(
            &fs::read_to_string(output_dir.join(SAMPLE_SUMMARY_JSON))
                .expect("sample_summary.json should be readable"),
        )
        .expect("sample_summary.json should parse");
        assert_eq!(written_summary, summary);
        assert_eq!(written_summary.decision_status, DecisionStatus::Positive);
        assert_eq!(written_summary.release_status, ReleaseStatus::Provisional);

        let written_manifest: RunManifest = serde_json::from_str(
            &fs::read_to_string(output_dir.join(RUN_MANIFEST_JSON))
                .expect("run_manifest.json should be readable"),
        )
        .expect("run_manifest.json should parse");
        assert_eq!(written_manifest, manifest);

        let mut candidate_reader = ReaderBuilder::new()
            .delimiter(b'\t')
            .from_path(output_dir.join(CANDIDATE_CALLS_TSV))
            .expect("candidate_calls.tsv should be readable");
        let headers = candidate_reader
            .headers()
            .expect("candidate_calls.tsv should have headers")
            .clone();
        assert_eq!(
            headers.iter().collect::<Vec<_>>(),
            vec![
                "virus_name",
                "taxid",
                "accession_or_group",
                "accepted_fragments",
                "nonoverlap_fragments",
                "raw_fraction",
                "unique_fraction",
                "fraction_ci_95",
                "clopper_pearson_upper",
                "breadth",
                "ambiguous_fragments",
                "background_ratio",
                "decision",
                "decision_reasons",
                "evidence_strength",
            ]
        );

        let candidate_records = candidate_reader
            .records()
            .collect::<std::result::Result<Vec<_>, _>>()
            .expect("candidate_calls.tsv records should parse");
        assert_eq!(candidate_records.len(), 2);

        let positive_record = candidate_records
            .iter()
            .find(|record| record.get(2) == Some("ebv"))
            .expect("positive candidate should be present");
        let decision_reasons: Vec<String> = serde_json::from_str(
            positive_record
                .get(13)
                .expect("decision_reasons column should exist"),
        )
        .expect("decision_reasons should stay JSON encoded");
        assert!(decision_reasons.iter().any(|reason| {
            reason == &format!("{FRACTION_CI_REASON_PREFIX}{SAMPLING_ONLY_CI_LABEL}")
        }));

        let fraction_ci_95: [f64; 2] = serde_json::from_str(
            positive_record
                .get(7)
                .expect("fraction_ci_95 column should exist"),
        )
        .expect("fraction_ci_95 should stay JSON encoded");
        assert_eq!(fraction_ci_95, candidates[0].fraction_ci_95);

        let mut rounds_reader = ReaderBuilder::new()
            .delimiter(b'\t')
            .from_path(output_dir.join(ROUNDS_TSV))
            .expect("rounds.tsv should be readable");
        let round_headers = rounds_reader
            .headers()
            .expect("rounds.tsv should have headers")
            .clone();
        assert_eq!(
            round_headers.iter().collect::<Vec<_>>(),
            vec!["sampled_fragments", "accepted_virus", "decision_status"]
        );
        let round_records = rounds_reader
            .records()
            .collect::<std::result::Result<Vec<_>, _>>()
            .expect("round rows should parse");
        assert_eq!(round_records.len(), 3);
        assert_eq!(round_records[2].get(2), Some("positive"));

        assert_eq!(
            fs::read_to_string(output_dir.join(COVERAGE_DIR).join("ebv.coverage.tsv"))
                .expect("coverage file should be readable"),
            COVERAGE_HEADER
        );
        assert!(
            !output_dir
                .join(COVERAGE_DIR)
                .join("influenza-a.coverage.tsv")
                .exists(),
            "negative candidate should not emit coverage"
        );

        let run_log = fs::read_to_string(output_dir.join(LOGS_DIR).join(RUN_LOG))
            .expect("run.log should exist");
        assert!(run_log.contains("\"decision_status\":\"positive\""));
        assert!(run_log.contains("\"release_status\":\"provisional\""));
        assert!(run_log.contains(SAMPLING_ONLY_CI_LABEL));
    }

    #[test]
    fn preserves_sampling_only_ci_label_in_candidate_calls_tsv() {
        let temp_dir = TempDir::new().expect("temp dir should be created");
        let output_dir = temp_dir.path().join("report_bundle");

        ReportWriter::write(
            &output_dir,
            &sample_summary_fixture(),
            &candidate_fixtures(),
            &round_fixtures(),
            &run_manifest_fixture(),
        )
        .expect("report bundle should be written");

        let mut reader = ReaderBuilder::new()
            .delimiter(b'\t')
            .from_path(output_dir.join(CANDIDATE_CALLS_TSV))
            .expect("candidate_calls.tsv should be readable");
        let rows = reader
            .records()
            .collect::<std::result::Result<Vec<_>, _>>()
            .expect("candidate rows should parse");

        let reasons: Vec<String> = serde_json::from_str(
            rows[0]
                .get(13)
                .expect("decision_reasons column should exist"),
        )
        .expect("decision_reasons should stay JSON encoded");
        let label = reasons
            .into_iter()
            .find(|reason| reason.starts_with(FRACTION_CI_REASON_PREFIX))
            .expect("sampling-only CI label should be present");

        assert_eq!(
            label,
            format!("{FRACTION_CI_REASON_PREFIX}{SAMPLING_ONLY_CI_LABEL}")
        );
    }

    #[test]
    fn rejects_duplicate_fraction_ci_labels() {
        let temp_dir = TempDir::new().expect("temp dir should be created");
        let output_dir = temp_dir.path().join("report_bundle");
        let mut candidates = candidate_fixtures();
        candidates[0]
            .decision_reasons
            .push(format!("{FRACTION_CI_REASON_PREFIX}{SAMPLING_ONLY_CI_LABEL}"));

        let error = ReportWriter::write(
            &output_dir,
            &sample_summary_fixture(),
            &candidates,
            &round_fixtures(),
            &run_manifest_fixture(),
        )
        .expect_err("duplicate CI labels should fail validation");

        assert!(
            error
                .to_string()
                .contains("must carry exactly one `fraction_ci_95_label=...` marker"),
            "{error}"
        );
    }

    fn collect_relative_paths(root: &Path) -> Vec<String> {
        let mut paths = Vec::new();
        collect_relative_paths_recursive(root, root, &mut paths);
        paths.sort();
        paths
    }

    fn collect_relative_paths_recursive(root: &Path, current: &Path, paths: &mut Vec<String>) {
        let mut entries = fs::read_dir(current)
            .expect("directory should be readable")
            .map(|entry| entry.expect("directory entry should be readable").path())
            .collect::<Vec<_>>();
        entries.sort();

        for path in entries {
            if path.is_dir() {
                collect_relative_paths_recursive(root, &path, paths);
            } else {
                paths.push(
                    path.strip_prefix(root)
                        .expect("path should stay under root")
                        .to_string_lossy()
                        .replace('\\', "/"),
                );
            }
        }
    }

    fn sample_summary_fixture() -> SampleSummary {
        SampleSummary {
            sample_id: "S001".into(),
            reference_bundle_version: "rvscreen_ref_2026.04.20-r1".into(),
            calibration_profile_version: "rvscreen_calib_2026.04.20-r1".into(),
            backend: "minimap2".into(),
            seed: 20_260_420,
            input_fragments: 24_837_219,
            qc_passing_fragments: 24_100_988,
            sampled_fragments: 400_000,
            rounds_run: 3,
            stop_reason: StopReason::PositiveBoundaryCrossed,
            decision_status: DecisionStatus::Positive,
            release_status: ReleaseStatus::Provisional,
        }
    }

    fn run_manifest_fixture() -> RunManifest {
        RunManifest {
            reference_bundle_version: "rvscreen_ref_2026.04.20-r1".into(),
            calibration_profile_version: "rvscreen_calib_2026.04.20-r1".into(),
            backend: "minimap2".into(),
            seed: 20_260_420,
            sampling_mode: "representative".into(),
            negative_control: crate::types::NegativeControlManifest {
                required: false,
                status: crate::types::NegativeControlStatus::Pass,
                control_id: Some("neg-001".into()),
                control_status: Some("pass".into()),
            },
            input_files: vec!["reads_R1.fastq.gz".into(), "reads_R2.fastq.gz".into()],
        }
    }

    fn round_fixtures() -> Vec<RoundRecord> {
        vec![
            RoundRecord {
                sampled_fragments: 50_000,
                accepted_virus: 2,
                decision_status: DecisionStatus::Indeterminate,
            },
            RoundRecord {
                sampled_fragments: 100_000,
                accepted_virus: 4,
                decision_status: DecisionStatus::Indeterminate,
            },
            RoundRecord {
                sampled_fragments: 400_000,
                accepted_virus: 27,
                decision_status: DecisionStatus::Positive,
            },
        ]
    }

    fn candidate_fixtures() -> Vec<CandidateCall> {
        vec![
            CandidateCall {
                virus_name: "Epstein-Barr virus".into(),
                taxid: 10_376,
                accession_or_group: "ebv".into(),
                accepted_fragments: 27,
                nonoverlap_fragments: 6,
                raw_fraction: 0.00012,
                unique_fraction: 0.00008,
                fraction_ci_95: [0.00005, 0.00011],
                clopper_pearson_upper: 0.00014,
                breadth: 0.0025,
                ambiguous_fragments: 3,
                background_ratio: 2.1,
                decision: DecisionStatus::Positive,
                decision_reasons: vec![
                    format!("{FRACTION_CI_REASON_PREFIX}{SAMPLING_ONLY_CI_LABEL}"),
                    "unique_fraction_above_theta_pos".into(),
                    "background_ratio_meets_minimum".into(),
                ],
                evidence_strength: crate::types::EvidenceStrength::High,
            },
            CandidateCall {
                virus_name: "Influenza A virus".into(),
                taxid: 11_320,
                accession_or_group: "influenza-a".into(),
                accepted_fragments: 1,
                nonoverlap_fragments: 1,
                raw_fraction: 0.00001,
                unique_fraction: 0.000004,
                fraction_ci_95: [0.0, 0.00002],
                clopper_pearson_upper: 0.00003,
                breadth: 0.0001,
                ambiguous_fragments: 2,
                background_ratio: 0.8,
                decision: DecisionStatus::Negative,
                decision_reasons: vec![
                    format!("{FRACTION_CI_REASON_PREFIX}{SAMPLING_ONLY_CI_LABEL}"),
                    "unique_fraction_below_theta_neg".into(),
                    "background_negative_safe".into(),
                ],
                evidence_strength: crate::types::EvidenceStrength::Low,
            },
        ]
    }
}
