use crate::error::{Result, RvScreenError};
use crate::report::contract::{
    sampling_only_ci_label, sampling_only_ci_label_violation, ChecksumInclusion,
    ReportArtifactKind, CANDIDATE_CALLS_TSV, CANONICAL_CANDIDATE_CALLS_ACCESSION_TSV,
    CANONICAL_CHECKSUM_SHA256, CANONICAL_POSITIVE_CANDIDATE_COVERAGE_DIR,
    CANONICAL_RELEASE_GATE_JSON, CANONICAL_RESULT_OVERVIEW_JSON, CANONICAL_RUN_MANIFEST_JSON,
    CANONICAL_RUN_NDJSON, CANONICAL_SAMPLE_RUN_SUMMARY_JSON, CANONICAL_SAMPLING_ROUNDS_TSV,
    CANONICAL_SCHEMA_VERSIONS_JSON, CANONICAL_VIRUS_SUMMARY_TSV, CHECKSUM_SHA256, COVERAGE_DIR,
    COVERAGE_HEADER, LOGS_DIR, REPORT_ARTIFACT_CONTRACTS, ROUNDS_TSV, RUN_LOG, RUN_MANIFEST_JSON,
    SAMPLE_SUMMARY_JSON,
};
use crate::report::summary::{
    serialize_result_overview_json, serialize_sample_run_summary_json, serialize_virus_summary_tsv,
};
use crate::types::{CandidateCall, DecisionStatus, RoundRecord, RunManifest, SampleSummary};
use csv::WriterBuilder;
use serde::Serialize;
use serde_json::json;
use sha2::{Digest, Sha256};
use std::collections::{BTreeMap, BTreeSet};
use std::fs::{self, File};
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

pub(crate) const CANDIDATE_CALLS_HEADER: [&str; 15] = [
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
];

pub(crate) const ROUNDS_HEADER: [&str; 3] =
    ["sampled_fragments", "accepted_virus", "decision_status"];

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
        clean_managed_outputs(output_dir)?;
        create_dir_all(&output_dir.join(COVERAGE_DIR))?;
        create_dir_all(&output_dir.join(CANONICAL_POSITIVE_CANDIDATE_COVERAGE_DIR))?;
        create_dir_all(&output_dir.join(LOGS_DIR))?;
        create_parent_dir(&output_dir.join(CANONICAL_SAMPLE_RUN_SUMMARY_JSON))?;
        create_parent_dir(&output_dir.join(CANONICAL_CANDIDATE_CALLS_ACCESSION_TSV))?;
        create_parent_dir(&output_dir.join(CANONICAL_RUN_MANIFEST_JSON))?;
        create_parent_dir(&output_dir.join(CANONICAL_CHECKSUM_SHA256))?;
        create_parent_dir(&output_dir.join(CANONICAL_SCHEMA_VERSIONS_JSON))?;
        create_parent_dir(&output_dir.join(CANONICAL_RELEASE_GATE_JSON))?;

        write_json_pretty(&output_dir.join(SAMPLE_SUMMARY_JSON), summary)?;
        write_bytes(
            &output_dir.join(CANONICAL_SAMPLE_RUN_SUMMARY_JSON),
            serialize_sample_run_summary_json(summary, candidates, rounds)?,
        )?;
        write_bytes(
            &output_dir.join(CANONICAL_RESULT_OVERVIEW_JSON),
            serialize_result_overview_json(summary, candidates, manifest)?,
        )?;
        write_bytes(
            &output_dir.join(CANONICAL_VIRUS_SUMMARY_TSV),
            serialize_virus_summary_tsv(candidates)?,
        )?;
        write_candidate_calls_tsv(&output_dir.join(CANDIDATE_CALLS_TSV), candidates)?;
        write_candidate_calls_tsv(
            &output_dir.join(CANONICAL_CANDIDATE_CALLS_ACCESSION_TSV),
            candidates,
        )?;
        write_json_pretty(&output_dir.join(RUN_MANIFEST_JSON), manifest)?;
        write_json_pretty(&output_dir.join(CANONICAL_RUN_MANIFEST_JSON), manifest)?;
        write_json_pretty(
            &output_dir.join(CANONICAL_SCHEMA_VERSIONS_JSON),
            &schema_versions(),
        )?;
        write_json_pretty(
            &output_dir.join(CANONICAL_RELEASE_GATE_JSON),
            &release_gate(summary, manifest)?,
        )?;
        write_rounds_tsv(&output_dir.join(ROUNDS_TSV), rounds)?;
        write_rounds_tsv(&output_dir.join(CANONICAL_SAMPLING_ROUNDS_TSV), rounds)?;
        write_coverage_files(&output_dir.join(COVERAGE_DIR), candidates)?;
        write_canonical_coverage_files(
            &output_dir.join(CANONICAL_POSITIVE_CANDIDATE_COVERAGE_DIR),
            candidates,
        )?;
        write_run_log(
            &output_dir.join(LOGS_DIR).join(RUN_LOG),
            output_dir,
            summary,
            candidates,
        )?;
        write_run_log(
            &output_dir.join(CANONICAL_RUN_NDJSON),
            output_dir,
            summary,
            candidates,
        )?;
        let checksum_files = collect_report_bundle_files(output_dir)?;
        write_checksum_file(
            &output_dir.join(CANONICAL_CHECKSUM_SHA256),
            output_dir,
            &checksum_files,
        )?;
        write_checksum_file(
            &output_dir.join(CHECKSUM_SHA256),
            output_dir,
            &checksum_files,
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
                format!("candidate `{}` {violation}", candidate.accession_or_group,),
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

fn write_bytes(path: &Path, bytes: Vec<u8>) -> Result<()> {
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
        .has_headers(false)
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

    writer
        .write_record(CANDIDATE_CALLS_HEADER)
        .map_err(|error| {
            RvScreenError::validation(
                "report.candidate_calls.tsv",
                format!(
                    "failed to serialize candidate_calls.tsv header for `{}`: {error}",
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
        .has_headers(false)
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

    writer.write_record(ROUNDS_HEADER).map_err(|error| {
        RvScreenError::validation(
            "report.rounds.tsv",
            format!(
                "failed to serialize rounds.tsv header for `{}`: {error}",
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

fn write_canonical_coverage_files(coverage_dir: &Path, candidates: &[CandidateCall]) -> Result<()> {
    let mut seen_paths = BTreeSet::new();

    for candidate in candidates
        .iter()
        .filter(|candidate| candidate.decision == DecisionStatus::Positive)
    {
        let file_name = format!(
            "{}.{}.coverage.tsv",
            canonical_coverage_prefix(candidate),
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

fn canonical_coverage_prefix(candidate: &CandidateCall) -> &'static str {
    if candidate
        .decision_reasons
        .iter()
        .any(|reason| reason == "aggregation_level=virus_group")
    {
        "group"
    } else {
        "accession"
    }
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

fn schema_versions() -> impl Serialize {
    let canonical_outputs = REPORT_ARTIFACT_CONTRACTS
        .iter()
        .map(|contract| (contract.key, contract.canonical_path))
        .collect::<BTreeMap<_, _>>();
    let legacy_outputs = REPORT_ARTIFACT_CONTRACTS
        .iter()
        .filter_map(|contract| {
            contract
                .legacy_alias_path
                .map(|legacy_alias_path| (contract.key, legacy_alias_path))
        })
        .collect::<BTreeMap<_, _>>();

    json!({
        "report_bundle_schema": "rvscreen.report_bundle.v2",
        "canonical_outputs": canonical_outputs,
        "legacy_outputs": legacy_outputs,
    })
}

fn release_gate(summary: &SampleSummary, manifest: &RunManifest) -> Result<impl Serialize> {
    Ok(json!({
        "schema_version": "rvscreen.release_gate.v1",
        "sample_id": summary.sample_id,
        "release_status": enum_label(&summary.release_status)?,
        "decision_status": enum_label(&summary.decision_status)?,
        "stop_reason": enum_label(&summary.stop_reason)?,
        "rounds_run": summary.rounds_run,
        "reference_bundle_version": summary.reference_bundle_version,
        "calibration_profile_version": summary.calibration_profile_version,
        "backend": summary.backend,
        "seed": summary.seed,
        "negative_control": manifest.negative_control,
        "input_file_count": manifest.input_files.len(),
    }))
}

fn create_dir_all(path: &Path) -> Result<()> {
    fs::create_dir_all(path).map_err(|error| RvScreenError::io(path, error))
}

fn create_parent_dir(path: &Path) -> Result<()> {
    if let Some(parent) = path.parent() {
        create_dir_all(parent)
    } else {
        Ok(())
    }
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
struct ManagedOutputPath {
    relative_path: &'static str,
    kind: ReportArtifactKind,
}

fn clean_managed_outputs(output_dir: &Path) -> Result<()> {
    for managed_path in managed_output_paths() {
        clean_managed_output(output_dir, managed_path)?;
    }

    Ok(())
}

fn managed_output_paths() -> Vec<ManagedOutputPath> {
    let mut paths = Vec::new();

    for contract in REPORT_ARTIFACT_CONTRACTS {
        paths.push(ManagedOutputPath {
            relative_path: contract.canonical_path,
            kind: contract.kind,
        });

        if let Some(legacy_alias_path) = contract.legacy_alias_path {
            paths.push(ManagedOutputPath {
                relative_path: legacy_alias_path,
                kind: contract.kind,
            });
        }
    }

    paths.sort_by(|left, right| left.relative_path.cmp(right.relative_path));
    paths.dedup_by(|left, right| {
        left.relative_path == right.relative_path && left.kind == right.kind
    });
    paths
}

fn clean_managed_output(output_dir: &Path, managed_path: ManagedOutputPath) -> Result<()> {
    let relative_path = Path::new(managed_path.relative_path);
    ensure_managed_parent_dirs(output_dir, relative_path)?;
    let path = output_dir.join(relative_path);

    match managed_path.kind {
        ReportArtifactKind::File => remove_managed_file(&path, managed_path.relative_path),
        ReportArtifactKind::Directory => {
            remove_managed_directory(&path, managed_path.relative_path)
        }
    }
}

fn ensure_managed_parent_dirs(output_dir: &Path, relative_path: &Path) -> Result<()> {
    let Some(parent) = relative_path.parent() else {
        return Ok(());
    };

    let mut current = output_dir.to_path_buf();
    for component in parent.components() {
        current.push(component.as_os_str());
        match fs::symlink_metadata(&current) {
            Ok(metadata) if metadata.is_dir() => {}
            Ok(_) => {
                return Err(managed_output_conflict(
                    &current,
                    "parent directory",
                    relative_path,
                ));
            }
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => break,
            Err(error) => return Err(RvScreenError::io(&current, error)),
        }
    }

    Ok(())
}

fn remove_managed_file(path: &Path, relative_path: &str) -> Result<()> {
    match fs::symlink_metadata(path) {
        Ok(metadata) if metadata.is_file() || metadata.file_type().is_symlink() => {
            fs::remove_file(path).map_err(|error| RvScreenError::io(path, error))
        }
        Ok(metadata) if metadata.is_dir() => Err(managed_output_conflict(
            path,
            "file",
            Path::new(relative_path),
        )),
        Ok(_) => Err(managed_output_conflict(
            path,
            "file",
            Path::new(relative_path),
        )),
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => Ok(()),
        Err(error) => Err(RvScreenError::io(path, error)),
    }
}

fn remove_managed_directory(path: &Path, relative_path: &str) -> Result<()> {
    match fs::symlink_metadata(path) {
        Ok(metadata) if metadata.is_dir() => {
            fs::remove_dir_all(path).map_err(|error| RvScreenError::io(path, error))
        }
        Ok(_) => Err(managed_output_conflict(
            path,
            "directory",
            Path::new(relative_path),
        )),
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => Ok(()),
        Err(error) => Err(RvScreenError::io(path, error)),
    }
}

fn managed_output_conflict(path: &Path, expected: &str, relative_path: &Path) -> RvScreenError {
    RvScreenError::validation(
        "report.cleanup",
        format!(
            "managed output `{}` requires a {expected}, but `{}` is an unmanaged conflicting path; refusing to delete arbitrary data",
            relative_path.display(),
            path.display()
        ),
    )
}

fn collect_report_bundle_files(output_dir: &Path) -> Result<Vec<PathBuf>> {
    let mut files = Vec::new();

    for managed_path in managed_output_paths() {
        collect_managed_report_bundle_files(output_dir, managed_path, &mut files)?;
    }

    files.sort();
    files.dedup();
    Ok(files)
}

fn collect_managed_report_bundle_files(
    output_dir: &Path,
    managed_path: ManagedOutputPath,
    files: &mut Vec<PathBuf>,
) -> Result<()> {
    let path = output_dir.join(managed_path.relative_path);
    match fs::symlink_metadata(&path) {
        Ok(metadata) if metadata.is_file() || metadata.file_type().is_symlink() => {
            if checksum_inclusion(output_dir, &path)? == ChecksumInclusion::Include {
                files.push(path);
            }
        }
        Ok(metadata) if metadata.is_dir() && managed_path.kind == ReportArtifactKind::Directory => {
            collect_report_bundle_files_recursive(output_dir, &path, files)?;
        }
        Ok(_) => {}
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => {}
        Err(error) => return Err(RvScreenError::io(&path, error)),
    }

    Ok(())
}

fn collect_report_bundle_files_recursive(
    output_dir: &Path,
    current: &Path,
    files: &mut Vec<PathBuf>,
) -> Result<()> {
    let mut entries = fs::read_dir(current)
        .map_err(|error| RvScreenError::io(current, error))?
        .collect::<std::result::Result<Vec<_>, _>>()
        .map_err(|error| RvScreenError::io(current, error))?;
    entries.sort_by_key(|entry| entry.path());

    for entry in entries {
        let path = entry.path();
        if path.is_dir() {
            collect_report_bundle_files_recursive(output_dir, &path, files)?;
        } else if checksum_inclusion(output_dir, &path)? == ChecksumInclusion::Include {
            files.push(path);
        }
    }

    Ok(())
}

fn write_checksum_file(checksum_path: &Path, bundle_dir: &Path, files: &[PathBuf]) -> Result<()> {
    let file =
        File::create(checksum_path).map_err(|error| RvScreenError::io(checksum_path, error))?;
    let mut writer = BufWriter::new(file);

    for path in files {
        let digest = sha256_file(path)?;
        let relative = relative_report_path(bundle_dir, path)?;
        writeln!(writer, "{digest}  {}", relative.display())
            .map_err(|error| RvScreenError::io(checksum_path, error))?;
    }

    writer
        .flush()
        .map_err(|error| RvScreenError::io(checksum_path, error))
}

fn checksum_inclusion(bundle_dir: &Path, path: &Path) -> Result<ChecksumInclusion> {
    checksum_inclusion_for_relative_path(relative_report_path(bundle_dir, path)?)
}

fn checksum_inclusion_for_relative_path(relative_path: &Path) -> Result<ChecksumInclusion> {
    if relative_path.file_name().and_then(|name| name.to_str()) == Some(CHECKSUM_SHA256) {
        return Ok(ChecksumInclusion::Exclude);
    }

    for contract in REPORT_ARTIFACT_CONTRACTS {
        if relative_path == Path::new(contract.canonical_path) {
            return Ok(contract.checksum_inclusion);
        }

        if let Some(legacy_alias_path) = contract.legacy_alias_path {
            if relative_path == Path::new(legacy_alias_path) {
                return Ok(contract.checksum_inclusion);
            }
        }
    }

    Ok(ChecksumInclusion::Include)
}

fn relative_report_path<'a>(bundle_dir: &'a Path, path: &'a Path) -> Result<&'a Path> {
    path.strip_prefix(bundle_dir).map_err(|_| {
        RvScreenError::validation(
            "report.checksum",
            format!(
                "`{}` is not inside `{}`",
                path.display(),
                bundle_dir.display()
            ),
        )
    })
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
    use crate::audit::run_audit_verify;
    use crate::cli::AuditVerifyArgs;
    use crate::decision::SAMPLING_ONLY_CI_LABEL;
    use crate::report::contract::{
        CANONICAL_CANDIDATE_CALLS_ACCESSION_TSV, CANONICAL_CHECKSUM_SHA256,
        CANONICAL_POSITIVE_CANDIDATE_COVERAGE_DIR, CANONICAL_RELEASE_GATE_JSON,
        CANONICAL_RESULT_OVERVIEW_JSON, CANONICAL_RUN_MANIFEST_JSON, CANONICAL_RUN_NDJSON,
        CANONICAL_SAMPLE_RUN_SUMMARY_JSON, CANONICAL_SAMPLING_ROUNDS_TSV,
        CANONICAL_SCHEMA_VERSIONS_JSON, CANONICAL_VIRUS_SUMMARY_TSV, FRACTION_CI_REASON_PREFIX,
        SUMMARY_DIR,
    };
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
                "audit/checksum.sha256".to_string(),
                "audit/release_gate.json".to_string(),
                "audit/schema_versions.json".to_string(),
                "evidence/positive_candidate_coverage/accession.ebv.coverage.tsv".to_string(),
                "evidence/positive_candidate_coverage/group.herpes_viridae.coverage.tsv"
                    .to_string(),
                "coverage/ebv.coverage.tsv".to_string(),
                "coverage/herpes_viridae.coverage.tsv".to_string(),
                "logs/run.log".to_string(),
                "logs/run.ndjson".to_string(),
                "provenance/run_manifest.json".to_string(),
                "rounds.tsv".to_string(),
                "results/candidate_calls.accession.tsv".to_string(),
                "results/sampling_rounds.tsv".to_string(),
                "summary/result_overview.json".to_string(),
                "summary/sample_run_summary.json".to_string(),
                "summary/virus_summary.tsv".to_string(),
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

        let canonical_sample_run_summary: serde_json::Value = serde_json::from_str(
            &fs::read_to_string(output_dir.join(CANONICAL_SAMPLE_RUN_SUMMARY_JSON))
                .expect("sample_run_summary.json should be readable"),
        )
        .expect("sample_run_summary.json should parse");
        assert_eq!(
            canonical_sample_run_summary["schema_version"],
            "rvscreen.sample_run_summary.v1"
        );
        assert_eq!(canonical_sample_run_summary["sample_id"], summary.sample_id);
        assert_eq!(canonical_sample_run_summary["positive_entity_count"], 2);

        let result_overview: serde_json::Value = serde_json::from_str(
            &fs::read_to_string(output_dir.join(CANONICAL_RESULT_OVERVIEW_JSON))
                .expect("result_overview.json should be readable"),
        )
        .expect("result_overview.json should parse");
        assert_eq!(
            result_overview["schema_version"],
            "rvscreen.result_overview.v1"
        );
        assert_eq!(result_overview["sample"]["sample_id"], summary.sample_id);

        let mut virus_summary_reader = ReaderBuilder::new()
            .delimiter(b'\t')
            .from_path(output_dir.join(CANONICAL_VIRUS_SUMMARY_TSV))
            .expect("virus_summary.tsv should be readable");
        assert_eq!(
            virus_summary_reader
                .headers()
                .expect("virus_summary.tsv should have headers")
                .iter()
                .collect::<Vec<_>>(),
            vec![
                "rank",
                "entity_id",
                "entity_name",
                "entity_type",
                "taxid",
                "aggregation_level",
                "decision",
                "evidence_strength",
                "accepted_fragments",
                "nonoverlap_fragments",
                "ambiguous_fragments",
                "raw_fraction",
                "unique_fraction",
                "fraction_ci_95_low",
                "fraction_ci_95_high",
                "clopper_pearson_upper",
                "breadth",
                "background_ratio",
                "top_accession",
                "top_accession_name",
                "supporting_accession_count",
                "decision_reasons",
            ]
        );
        assert_eq!(
            virus_summary_reader
                .records()
                .collect::<std::result::Result<Vec<_>, _>>()
                .expect("virus summary rows should parse")
                .len(),
            3
        );

        let written_manifest: RunManifest = serde_json::from_str(
            &fs::read_to_string(output_dir.join(RUN_MANIFEST_JSON))
                .expect("run_manifest.json should be readable"),
        )
        .expect("run_manifest.json should parse");
        assert_eq!(written_manifest, manifest);

        let canonical_manifest: RunManifest = serde_json::from_str(
            &fs::read_to_string(output_dir.join(CANONICAL_RUN_MANIFEST_JSON))
                .expect("canonical run manifest should be readable"),
        )
        .expect("canonical run manifest should parse");
        assert_eq!(canonical_manifest, manifest);
        assert_eq!(
            fs::read(output_dir.join(RUN_MANIFEST_JSON)).expect("legacy manifest should read"),
            fs::read(output_dir.join(CANONICAL_RUN_MANIFEST_JSON))
                .expect("canonical manifest should read")
        );

        let schema_versions: serde_json::Value = serde_json::from_str(
            &fs::read_to_string(output_dir.join(CANONICAL_SCHEMA_VERSIONS_JSON))
                .expect("schema_versions.json should be readable"),
        )
        .expect("schema_versions.json should parse");
        assert_eq!(
            schema_versions["report_bundle_schema"],
            "rvscreen.report_bundle.v2"
        );
        for contract in REPORT_ARTIFACT_CONTRACTS {
            assert_eq!(
                schema_versions["canonical_outputs"][contract.key], contract.canonical_path,
                "canonical output for {} should come from contract",
                contract.key
            );
            if let Some(legacy_alias_path) = contract.legacy_alias_path {
                assert_eq!(
                    schema_versions["legacy_outputs"][contract.key], legacy_alias_path,
                    "legacy output for {} should come from contract",
                    contract.key
                );
            } else {
                assert!(
                    schema_versions["legacy_outputs"]
                        .get(contract.key)
                        .is_none(),
                    "{} should not have a legacy output",
                    contract.key
                );
            }
        }

        let release_gate: serde_json::Value = serde_json::from_str(
            &fs::read_to_string(output_dir.join(CANONICAL_RELEASE_GATE_JSON))
                .expect("release_gate.json should be readable"),
        )
        .expect("release_gate.json should parse");
        assert_eq!(release_gate["schema_version"], "rvscreen.release_gate.v1");
        assert_eq!(release_gate["sample_id"], summary.sample_id);
        assert_eq!(release_gate["release_status"], "provisional");
        assert_eq!(release_gate["decision_status"], "positive");
        assert_eq!(release_gate["stop_reason"], "positive_boundary_crossed");
        assert_eq!(release_gate["negative_control"]["status"], "pass");
        assert_eq!(release_gate["input_file_count"], manifest.input_files.len());

        assert_eq!(
            fs::read(output_dir.join(CANDIDATE_CALLS_TSV))
                .expect("legacy candidate calls should read"),
            fs::read(output_dir.join(CANONICAL_CANDIDATE_CALLS_ACCESSION_TSV))
                .expect("canonical candidate calls should read")
        );
        assert_eq!(
            fs::read(output_dir.join(ROUNDS_TSV)).expect("legacy rounds should read"),
            fs::read(output_dir.join(CANONICAL_SAMPLING_ROUNDS_TSV))
                .expect("canonical rounds should read")
        );

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
        assert_eq!(candidate_records.len(), 3);

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

        let accepted_fragments = positive_record
            .get(3)
            .expect("accepted_fragments column should exist")
            .parse::<u64>()
            .expect("accepted_fragments should parse");
        let raw_fraction = positive_record
            .get(5)
            .expect("raw_fraction column should exist")
            .parse::<f64>()
            .expect("raw_fraction should parse");
        let unique_fraction = positive_record
            .get(6)
            .expect("unique_fraction column should exist")
            .parse::<f64>()
            .expect("unique_fraction should parse");
        let sampled_fraction = accepted_fragments as f64 / written_summary.sampled_fragments as f64;
        assert!((raw_fraction - sampled_fraction).abs() < 1e-12);
        assert!((unique_fraction - sampled_fraction).abs() < 1e-12);

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
        assert_eq!(
            fs::read_to_string(
                output_dir
                    .join(COVERAGE_DIR)
                    .join("herpes_viridae.coverage.tsv")
            )
            .expect("group legacy coverage file should be readable"),
            COVERAGE_HEADER
        );
        assert_eq!(
            fs::read_to_string(
                output_dir
                    .join(CANONICAL_POSITIVE_CANDIDATE_COVERAGE_DIR)
                    .join("accession.ebv.coverage.tsv")
            )
            .expect("canonical accession coverage file should be readable"),
            COVERAGE_HEADER
        );
        assert_eq!(
            fs::read_to_string(
                output_dir
                    .join(CANONICAL_POSITIVE_CANDIDATE_COVERAGE_DIR)
                    .join("group.herpes_viridae.coverage.tsv")
            )
            .expect("canonical group coverage file should be readable"),
            COVERAGE_HEADER
        );
        assert_eq!(
            coverage_stems(&output_dir.join(COVERAGE_DIR), false),
            coverage_stems(
                &output_dir.join(CANONICAL_POSITIVE_CANDIDATE_COVERAGE_DIR),
                true
            )
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
        let run_ndjson = fs::read_to_string(output_dir.join(CANONICAL_RUN_NDJSON))
            .expect("run.ndjson should exist");
        let ndjson_events = run_ndjson
            .lines()
            .filter(|line| !line.is_empty())
            .map(|line| {
                let event = serde_json::from_str::<serde_json::Value>(line).expect("line parses");
                assert!(event.get("timestamp_unix_ms").is_some());
                assert!(event.get("level").is_some());
                assert!(event.get("event").is_some());
                event
            })
            .collect::<Vec<_>>();
        assert_eq!(ndjson_events.len(), 2);
        assert_eq!(ndjson_events[0]["event"], "report_bundle_write_started");
        assert_eq!(ndjson_events[1]["event"], "report_bundle_write_completed");

        let root_checksum = checksum_manifest_paths(
            &fs::read_to_string(output_dir.join(CHECKSUM_SHA256))
                .expect("root checksum should read"),
        );
        let canonical_checksum = checksum_manifest_paths(
            &fs::read_to_string(output_dir.join(CANONICAL_CHECKSUM_SHA256))
                .expect("canonical checksum should read"),
        );
        assert_eq!(root_checksum, canonical_checksum);
        assert!(root_checksum.contains(&CANONICAL_SAMPLE_RUN_SUMMARY_JSON.to_string()));
        assert!(root_checksum.contains(&CANONICAL_SCHEMA_VERSIONS_JSON.to_string()));
        assert!(root_checksum.contains(&CANONICAL_RELEASE_GATE_JSON.to_string()));
        assert!(root_checksum.contains(&CANONICAL_RUN_NDJSON.to_string()));
        assert!(!root_checksum.contains(&CHECKSUM_SHA256.to_string()));
        assert!(!root_checksum.contains(&CANONICAL_CHECKSUM_SHA256.to_string()));
    }

    #[test]
    fn writes_no_detection_bundle_with_valid_empty_summaries_and_audit() {
        let temp_dir = TempDir::new().expect("temp dir should be created");
        let output_dir = temp_dir.path().join("no_detection_bundle");
        let summary = no_detection_sample_summary_fixture();
        let manifest = run_manifest_fixture();
        let candidates = Vec::new();
        let rounds = Vec::new();

        ReportWriter::write(&output_dir, &summary, &candidates, &rounds, &manifest)
            .expect("no-detection report bundle should be written");

        let files = collect_relative_paths(&output_dir);
        assert!(files.contains(&CANONICAL_SAMPLE_RUN_SUMMARY_JSON.to_string()));
        assert!(files.contains(&CANONICAL_RESULT_OVERVIEW_JSON.to_string()));
        assert!(files.contains(&CANONICAL_VIRUS_SUMMARY_TSV.to_string()));
        assert!(files.contains(&CANONICAL_CANDIDATE_CALLS_ACCESSION_TSV.to_string()));
        assert!(files.contains(&CANONICAL_RUN_MANIFEST_JSON.to_string()));
        assert!(!files.iter().any(|path| path.ends_with(".coverage.tsv")));

        let sample_run_summary: serde_json::Value = serde_json::from_str(
            &fs::read_to_string(output_dir.join(CANONICAL_SAMPLE_RUN_SUMMARY_JSON))
                .expect("sample_run_summary.json should be readable"),
        )
        .expect("sample_run_summary.json should parse");
        assert_eq!(sample_run_summary["decision_status"], "negative");
        assert_eq!(sample_run_summary["rounds_run"], 0);
        assert_eq!(sample_run_summary["positive_entity_count"], 0);
        assert_eq!(sample_run_summary["indeterminate_entity_count"], 0);
        assert_eq!(sample_run_summary["negative_entity_count"], 0);

        let result_overview: serde_json::Value = serde_json::from_str(
            &fs::read_to_string(output_dir.join(CANONICAL_RESULT_OVERVIEW_JSON))
                .expect("result_overview.json should be readable"),
        )
        .expect("result_overview.json should parse");
        assert_eq!(
            result_overview["top_results"].as_array().map(Vec::len),
            Some(0)
        );

        let mut virus_summary_reader = ReaderBuilder::new()
            .delimiter(b'\t')
            .from_path(output_dir.join(CANONICAL_VIRUS_SUMMARY_TSV))
            .expect("virus_summary.tsv should be readable");
        assert_eq!(
            virus_summary_reader
                .headers()
                .expect("virus_summary.tsv should have headers")
                .iter()
                .collect::<Vec<_>>(),
            crate::report::summary::VIRUS_SUMMARY_HEADER
        );
        assert_eq!(
            virus_summary_reader
                .records()
                .collect::<std::result::Result<Vec<_>, _>>()
                .expect("virus summary records should parse")
                .len(),
            0
        );

        let mut candidate_reader = ReaderBuilder::new()
            .delimiter(b'\t')
            .from_path(output_dir.join(CANONICAL_CANDIDATE_CALLS_ACCESSION_TSV))
            .expect("candidate TSV should be readable");
        assert_eq!(
            candidate_reader
                .headers()
                .expect("candidate TSV should have headers")
                .iter()
                .collect::<Vec<_>>(),
            CANDIDATE_CALLS_HEADER
        );
        assert_eq!(
            candidate_reader
                .records()
                .collect::<std::result::Result<Vec<_>, _>>()
                .expect("candidate rows should parse")
                .len(),
            0
        );

        let audit_report = run_audit_verify(&AuditVerifyArgs {
            report_bundle: output_dir,
            reference_bundle: None,
            calibration_profile: None,
        })
        .expect("audit verify should run for no-detection bundle");
        assert!(audit_report.passed(), "{}", audit_report.render());
        assert!(audit_report.render().contains("parsed 0 candidate row(s)"));
        assert!(audit_report
            .render()
            .contains("parsed 0 virus summary row(s)"));
    }

    #[test]
    fn summary_entities_trace_to_candidate_rows_evidence_and_provenance() {
        let temp_dir = TempDir::new().expect("temp dir should be created");
        let output_dir = temp_dir.path().join("report_bundle");
        let summary = sample_summary_fixture();
        let manifest = run_manifest_fixture();
        let candidates = candidate_fixtures();

        ReportWriter::write(
            &output_dir,
            &summary,
            &candidates,
            &round_fixtures(),
            &manifest,
        )
        .expect("report bundle should be written");

        let mut candidate_reader = ReaderBuilder::new()
            .delimiter(b'\t')
            .from_path(output_dir.join(CANONICAL_CANDIDATE_CALLS_ACCESSION_TSV))
            .expect("canonical candidate TSV should be readable");
        let candidate_rows = candidate_reader
            .records()
            .collect::<std::result::Result<Vec<_>, _>>()
            .expect("candidate rows should parse");
        let candidate_ids = candidate_rows
            .iter()
            .map(|row| row.get(2).expect("candidate id column should exist"))
            .collect::<BTreeSet<_>>();

        let mut virus_summary_reader = ReaderBuilder::new()
            .delimiter(b'\t')
            .from_path(output_dir.join(CANONICAL_VIRUS_SUMMARY_TSV))
            .expect("virus summary should be readable");
        let virus_summary_rows = virus_summary_reader
            .records()
            .collect::<std::result::Result<Vec<_>, _>>()
            .expect("virus summary rows should parse");
        for row in &virus_summary_rows {
            assert!(
                candidate_ids.contains(row.get(1).expect("entity id column should exist")),
                "summary entity should be traceable to a candidate row: {row:?}"
            );
        }

        let ebv_summary = virus_summary_rows
            .iter()
            .find(|row| row.get(1) == Some("ebv"))
            .expect("accession summary row should exist");
        assert_eq!(ebv_summary.get(18), Some("ebv"));
        assert_eq!(ebv_summary.get(20), Some("1"));
        assert!(output_dir
            .join(COVERAGE_DIR)
            .join("ebv.coverage.tsv")
            .is_file());
        assert!(output_dir
            .join(CANONICAL_POSITIVE_CANDIDATE_COVERAGE_DIR)
            .join("accession.ebv.coverage.tsv")
            .is_file());

        let group_summary = virus_summary_rows
            .iter()
            .find(|row| row.get(1) == Some("herpes viridae"))
            .expect("group summary row should exist");
        assert_eq!(group_summary.get(5), Some("virus_group"));
        assert_eq!(group_summary.get(18), Some("NC_007605"));
        assert_eq!(group_summary.get(19), Some("Human herpesvirus 4"));
        assert_eq!(group_summary.get(20), Some("1"));
        assert!(output_dir
            .join(COVERAGE_DIR)
            .join("herpes_viridae.coverage.tsv")
            .is_file());
        assert!(output_dir
            .join(CANONICAL_POSITIVE_CANDIDATE_COVERAGE_DIR)
            .join("group.herpes_viridae.coverage.tsv")
            .is_file());

        let result_overview: serde_json::Value = serde_json::from_str(
            &fs::read_to_string(output_dir.join(CANONICAL_RESULT_OVERVIEW_JSON))
                .expect("result_overview.json should be readable"),
        )
        .expect("result_overview.json should parse");
        for top_result in result_overview["top_results"]
            .as_array()
            .expect("top_results should be an array")
        {
            let entity_id = top_result["entity_id"]
                .as_str()
                .expect("top result entity_id should be a string");
            assert!(candidate_ids.contains(entity_id));
        }

        let canonical_manifest: RunManifest = serde_json::from_str(
            &fs::read_to_string(output_dir.join(CANONICAL_RUN_MANIFEST_JSON))
                .expect("canonical run manifest should be readable"),
        )
        .expect("canonical run manifest should parse");
        assert_eq!(canonical_manifest, manifest);
        assert!(output_dir.join(CANONICAL_CHECKSUM_SHA256).is_file());
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
        candidates[0].decision_reasons.push(format!(
            "{FRACTION_CI_REASON_PREFIX}{SAMPLING_ONLY_CI_LABEL}"
        ));

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

    #[test]
    fn rerun_cleans_managed_outputs_and_preserves_unmanaged_sentinels() {
        let temp_dir = TempDir::new().expect("temp dir should be created");
        let output_dir = temp_dir.path().join("report_bundle");

        ReportWriter::write(
            &output_dir,
            &sample_summary_fixture(),
            &candidate_fixtures(),
            &round_fixtures(),
            &run_manifest_fixture(),
        )
        .expect("initial report bundle should be written");

        write_test_file(
            &output_dir.join(COVERAGE_DIR).join("stale.coverage.tsv"),
            "stale",
        );
        write_test_file(
            &output_dir.join(CANONICAL_SAMPLE_RUN_SUMMARY_JSON),
            "stale canonical summary",
        );
        write_test_file(
            &output_dir.join(CANONICAL_RESULT_OVERVIEW_JSON),
            "stale canonical overview",
        );
        write_test_file(
            &output_dir.join(CANONICAL_VIRUS_SUMMARY_TSV),
            "stale canonical virus summary",
        );
        write_test_file(
            &output_dir.join(CANONICAL_CANDIDATE_CALLS_ACCESSION_TSV),
            "stale canonical candidate calls",
        );
        write_test_file(
            &output_dir.join(CANONICAL_SAMPLING_ROUNDS_TSV),
            "stale canonical sampling rounds",
        );
        write_test_file(
            &output_dir
                .join(CANONICAL_POSITIVE_CANDIDATE_COVERAGE_DIR)
                .join("stale.coverage.tsv"),
            "stale canonical coverage",
        );
        write_test_file(
            &output_dir.join(CANONICAL_RUN_MANIFEST_JSON),
            "stale canonical manifest",
        );
        write_test_file(
            &output_dir.join(CANONICAL_CHECKSUM_SHA256),
            "stale canonical checksum",
        );
        write_test_file(
            &output_dir.join(CANONICAL_SCHEMA_VERSIONS_JSON),
            "stale canonical schema versions",
        );
        write_test_file(
            &output_dir.join(CANONICAL_RELEASE_GATE_JSON),
            "stale canonical release gate",
        );
        write_test_file(
            &output_dir.join(CANONICAL_RUN_NDJSON),
            "stale canonical log",
        );

        let root_sentinel = output_dir.join("notes.txt");
        let summary_sentinel = output_dir.join(SUMMARY_DIR).join("notes.txt");
        let log_sentinel = output_dir.join(LOGS_DIR).join("operator.log");
        write_test_file(&root_sentinel, "keep root sentinel");
        write_test_file(&summary_sentinel, "keep summary sentinel");
        write_test_file(&log_sentinel, "keep log sentinel");

        ReportWriter::write(
            &output_dir,
            &sample_summary_fixture(),
            &candidate_fixtures(),
            &round_fixtures(),
            &run_manifest_fixture(),
        )
        .expect("rerun report bundle should be written");

        assert!(output_dir.join(SAMPLE_SUMMARY_JSON).exists());
        assert!(output_dir.join(CANDIDATE_CALLS_TSV).exists());
        assert!(output_dir.join(RUN_MANIFEST_JSON).exists());
        assert!(output_dir.join(ROUNDS_TSV).exists());
        assert!(output_dir.join(CHECKSUM_SHA256).exists());
        assert!(output_dir
            .join(COVERAGE_DIR)
            .join("ebv.coverage.tsv")
            .exists());
        assert!(output_dir.join(LOGS_DIR).join(RUN_LOG).exists());

        assert!(!output_dir
            .join(COVERAGE_DIR)
            .join("stale.coverage.tsv")
            .exists());
        assert!(output_dir.join(CANONICAL_SAMPLE_RUN_SUMMARY_JSON).exists());
        assert!(output_dir.join(CANONICAL_RESULT_OVERVIEW_JSON).exists());
        assert!(output_dir.join(CANONICAL_VIRUS_SUMMARY_TSV).exists());
        assert!(output_dir
            .join(CANONICAL_CANDIDATE_CALLS_ACCESSION_TSV)
            .exists());
        assert!(output_dir.join(CANONICAL_SAMPLING_ROUNDS_TSV).exists());
        assert!(!output_dir
            .join(CANONICAL_POSITIVE_CANDIDATE_COVERAGE_DIR)
            .join("stale.coverage.tsv")
            .exists());
        assert!(output_dir
            .join(CANONICAL_POSITIVE_CANDIDATE_COVERAGE_DIR)
            .join("accession.ebv.coverage.tsv")
            .exists());
        assert!(output_dir
            .join(CANONICAL_POSITIVE_CANDIDATE_COVERAGE_DIR)
            .join("group.herpes_viridae.coverage.tsv")
            .exists());
        assert!(output_dir.join(CANONICAL_RUN_MANIFEST_JSON).exists());
        assert!(output_dir.join(CANONICAL_CHECKSUM_SHA256).exists());
        assert!(output_dir.join(CANONICAL_SCHEMA_VERSIONS_JSON).exists());
        assert!(output_dir.join(CANONICAL_RELEASE_GATE_JSON).exists());
        assert!(output_dir.join(CANONICAL_RUN_NDJSON).exists());

        assert_eq!(
            fs::read_to_string(root_sentinel).expect("root sentinel should remain"),
            "keep root sentinel"
        );
        assert_eq!(
            fs::read_to_string(summary_sentinel).expect("summary sentinel should remain"),
            "keep summary sentinel"
        );
        assert_eq!(
            fs::read_to_string(log_sentinel).expect("log sentinel should remain"),
            "keep log sentinel"
        );
    }

    #[test]
    fn refuses_unmanaged_file_that_blocks_managed_directory() {
        let temp_dir = TempDir::new().expect("temp dir should be created");
        let output_dir = temp_dir.path().join("report_bundle");
        fs::create_dir_all(&output_dir).expect("output dir should be created");
        let sentinel = output_dir.join("notes.txt");
        write_test_file(&sentinel, "keep sentinel");
        write_test_file(&output_dir.join(COVERAGE_DIR), "not a directory");

        let error = ReportWriter::write(
            &output_dir,
            &sample_summary_fixture(),
            &candidate_fixtures(),
            &round_fixtures(),
            &run_manifest_fixture(),
        )
        .expect_err("managed directory conflict should fail");

        let message = error.to_string();
        assert!(message.contains("report.cleanup"), "{message}");
        assert!(
            message.contains("managed output `coverage` requires a directory"),
            "{message}"
        );
        assert!(
            message.contains("refusing to delete arbitrary data"),
            "{message}"
        );
        assert_eq!(
            fs::read_to_string(&sentinel).expect("sentinel should remain"),
            "keep sentinel"
        );
        assert_eq!(
            fs::read_to_string(output_dir.join(COVERAGE_DIR))
                .expect("conflicting file should remain"),
            "not a directory"
        );
    }

    #[test]
    fn checksum_manifest_excludes_root_and_nested_checksum_files_by_name() {
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

        write_test_file(
            &output_dir.join(CANONICAL_SAMPLE_RUN_SUMMARY_JSON),
            "future canonical summary",
        );
        write_test_file(
            &output_dir.join(CANONICAL_CHECKSUM_SHA256),
            "future canonical checksum",
        );
        write_test_file(&output_dir.join("notes.txt"), "root sidecar");
        write_test_file(
            &output_dir.join(SUMMARY_DIR).join("notes.txt"),
            "summary sidecar",
        );
        write_test_file(
            &output_dir.join(LOGS_DIR).join("operator.log"),
            "log sidecar",
        );

        let files = collect_report_bundle_files(&output_dir)
            .expect("report bundle files should be collected");
        let relative_paths = files
            .iter()
            .map(|path| {
                path.strip_prefix(&output_dir)
                    .expect("collected path should be under output dir")
                    .to_string_lossy()
                    .replace('\\', "/")
            })
            .collect::<Vec<_>>();
        assert!(!relative_paths
            .iter()
            .any(|path| path.ends_with(CHECKSUM_SHA256)));

        write_checksum_file(&output_dir.join(CHECKSUM_SHA256), &output_dir, &files)
            .expect("checksum manifest should be rewritten");
        let checksum_manifest = fs::read_to_string(output_dir.join(CHECKSUM_SHA256))
            .expect("checksum manifest should be readable");
        let checksum_paths = checksum_manifest_paths(&checksum_manifest);

        assert!(!checksum_paths
            .iter()
            .any(|path| path.ends_with(CHECKSUM_SHA256)));
        assert!(checksum_paths.contains(&SAMPLE_SUMMARY_JSON.to_string()));
        assert!(checksum_paths.contains(&CANDIDATE_CALLS_TSV.to_string()));
        assert!(checksum_paths.contains(&RUN_MANIFEST_JSON.to_string()));
        assert!(checksum_paths.contains(&ROUNDS_TSV.to_string()));
        assert!(checksum_paths.contains(&format!("{COVERAGE_DIR}/ebv.coverage.tsv")));
        assert!(checksum_paths.contains(&format!("{LOGS_DIR}/{RUN_LOG}")));
        assert!(checksum_paths.contains(&CANONICAL_SAMPLE_RUN_SUMMARY_JSON.to_string()));
        assert!(checksum_paths.contains(&CANONICAL_SCHEMA_VERSIONS_JSON.to_string()));
        assert!(checksum_paths.contains(&CANONICAL_RELEASE_GATE_JSON.to_string()));
        assert!(!checksum_paths.contains(&CHECKSUM_SHA256.to_string()));
        assert!(!checksum_paths.contains(&CANONICAL_CHECKSUM_SHA256.to_string()));
        assert!(!checksum_paths.contains(&"notes.txt".to_string()));
        assert!(!checksum_paths.contains(&format!("{SUMMARY_DIR}/notes.txt")));
        assert!(!checksum_paths.contains(&format!("{LOGS_DIR}/operator.log")));
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

    fn write_test_file(path: &Path, contents: &str) {
        if let Some(parent) = path.parent() {
            fs::create_dir_all(parent).expect("test file parent should be created");
        }
        fs::write(path, contents).expect("test file should be written");
    }

    fn coverage_stems(coverage_dir: &Path, strip_canonical_prefix: bool) -> BTreeSet<String> {
        let mut entries = fs::read_dir(coverage_dir)
            .expect("coverage directory should be readable")
            .map(|entry| {
                entry
                    .expect("coverage directory entry should be readable")
                    .file_name()
                    .to_string_lossy()
                    .into_owned()
            })
            .collect::<Vec<_>>();
        entries.sort();

        entries
            .into_iter()
            .map(|file_name| {
                let without_suffix = file_name
                    .strip_suffix(".coverage.tsv")
                    .expect("coverage file should use expected suffix");
                if strip_canonical_prefix {
                    without_suffix
                        .split_once('.')
                        .expect("canonical coverage should use explicit prefix")
                        .1
                        .to_string()
                } else {
                    without_suffix.to_string()
                }
            })
            .collect()
    }

    fn checksum_manifest_paths(manifest: &str) -> Vec<String> {
        manifest
            .lines()
            .map(|line| {
                line.split_once("  ")
                    .expect("checksum line should contain digest/path separator")
                    .1
                    .to_string()
            })
            .collect()
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

    fn no_detection_sample_summary_fixture() -> SampleSummary {
        SampleSummary {
            sample_id: "S-no-detection".into(),
            reference_bundle_version: "rvscreen_ref_2026.04.20-r1".into(),
            calibration_profile_version: "rvscreen_calib_2026.04.20-r1".into(),
            backend: "minimap2".into(),
            seed: 20_260_420,
            input_fragments: 24_837_219,
            qc_passing_fragments: 24_100_988,
            sampled_fragments: 400_000,
            rounds_run: 0,
            stop_reason: StopReason::NegativeBoundaryConfirmed,
            decision_status: DecisionStatus::Negative,
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
            sampling_round_plan: None,
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
                raw_fraction: 27.0 / 400_000.0,
                unique_fraction: 27.0 / 400_000.0,
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
                raw_fraction: 1.0 / 400_000.0,
                unique_fraction: 1.0 / 400_000.0,
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
            CandidateCall {
                virus_name: "Herpesviridae".into(),
                taxid: 10_000,
                accession_or_group: "herpes viridae".into(),
                accepted_fragments: 14,
                nonoverlap_fragments: 5,
                raw_fraction: 14.0 / 400_000.0,
                unique_fraction: 14.0 / 400_000.0,
                fraction_ci_95: [0.00003, 0.00008],
                clopper_pearson_upper: 0.0001,
                breadth: 0.0015,
                ambiguous_fragments: 4,
                background_ratio: 1.9,
                decision: DecisionStatus::Positive,
                decision_reasons: vec![
                    format!("{FRACTION_CI_REASON_PREFIX}{SAMPLING_ONLY_CI_LABEL}"),
                    "aggregation_level=virus_group".into(),
                    "supporting_accession=NC_007605".into(),
                    "supporting_accession_name=Human herpesvirus 4".into(),
                ],
                evidence_strength: crate::types::EvidenceStrength::Medium,
            },
        ]
    }
}
