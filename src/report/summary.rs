use crate::error::{Result, RvScreenError};
use crate::types::{
    CandidateCall, DecisionStatus, EvidenceStrength, RoundRecord, RunManifest, SampleSummary,
};
use csv::WriterBuilder;
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};

pub(crate) const SAMPLE_RUN_SUMMARY_SCHEMA_VERSION: &str = "rvscreen.sample_run_summary.v1";
pub(crate) const RESULT_OVERVIEW_SCHEMA_VERSION: &str = "rvscreen.result_overview.v1";
pub(crate) const VIRUS_SUMMARY_HEADER: [&str; 22] = [
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
];

#[derive(Debug, Clone, PartialEq, Serialize)]
pub(crate) struct SampleRunSummary {
    schema_version: &'static str,
    sample_id: String,
    decision_status: String,
    release_status: String,
    stop_reason: String,
    input_fragments: u64,
    qc_passing_fragments: u64,
    sampled_fragments: u64,
    rounds_run: u64,
    positive_entity_count: usize,
    indeterminate_entity_count: usize,
    negative_entity_count: usize,
}

#[derive(Debug, Clone, PartialEq, Serialize)]
pub(crate) struct ResultOverview {
    schema_version: &'static str,
    sample: ResultOverviewSample,
    run: ResultOverviewRun,
    top_results: Vec<ResultOverviewTopResult>,
}

#[derive(Debug, Clone, PartialEq, Serialize)]
pub(crate) struct ResultOverviewSample {
    sample_id: String,
    decision_status: String,
    release_status: String,
}

#[derive(Debug, Clone, PartialEq, Serialize)]
pub(crate) struct ResultOverviewRun {
    reference_bundle_version: String,
    calibration_profile_version: String,
    backend: String,
    sampling_mode: String,
    negative_control_status: String,
}

#[derive(Debug, Clone, PartialEq, Serialize)]
pub(crate) struct ResultOverviewTopResult {
    rank: usize,
    entity_id: String,
    entity_name: String,
    aggregation_level: &'static str,
    decision: String,
    evidence_strength: String,
    accepted_fragments: u64,
    unique_fraction: f64,
    breadth: f64,
    background_ratio: f64,
}

#[derive(Debug, Clone, PartialEq, Serialize)]
pub(crate) struct VirusSummaryRow {
    rank: usize,
    entity_id: String,
    entity_name: String,
    entity_type: &'static str,
    taxid: u64,
    aggregation_level: &'static str,
    decision: String,
    evidence_strength: String,
    accepted_fragments: u64,
    nonoverlap_fragments: u64,
    ambiguous_fragments: u64,
    raw_fraction: f64,
    unique_fraction: f64,
    fraction_ci_95_low: f64,
    fraction_ci_95_high: f64,
    clopper_pearson_upper: f64,
    breadth: f64,
    background_ratio: f64,
    top_accession: String,
    top_accession_name: String,
    supporting_accession_count: u64,
    decision_reasons: String,
}

pub(crate) fn sample_run_summary(
    summary: &SampleSummary,
    candidates: &[CandidateCall],
    rounds: &[RoundRecord],
) -> Result<SampleRunSummary> {
    Ok(SampleRunSummary {
        schema_version: SAMPLE_RUN_SUMMARY_SCHEMA_VERSION,
        sample_id: summary.sample_id.clone(),
        decision_status: enum_label(&summary.decision_status)?,
        release_status: enum_label(&summary.release_status)?,
        stop_reason: enum_label(&summary.stop_reason)?,
        input_fragments: summary.input_fragments,
        qc_passing_fragments: summary.qc_passing_fragments,
        sampled_fragments: summary.sampled_fragments,
        rounds_run: rounds.len() as u64,
        positive_entity_count: decision_count(candidates, DecisionStatus::Positive),
        indeterminate_entity_count: decision_count(candidates, DecisionStatus::Indeterminate),
        negative_entity_count: decision_count(candidates, DecisionStatus::Negative),
    })
}

pub(crate) fn result_overview(
    summary: &SampleSummary,
    candidates: &[CandidateCall],
    manifest: &RunManifest,
) -> Result<ResultOverview> {
    Ok(ResultOverview {
        schema_version: RESULT_OVERVIEW_SCHEMA_VERSION,
        sample: ResultOverviewSample {
            sample_id: summary.sample_id.clone(),
            decision_status: enum_label(&summary.decision_status)?,
            release_status: enum_label(&summary.release_status)?,
        },
        run: ResultOverviewRun {
            reference_bundle_version: manifest.reference_bundle_version.clone(),
            calibration_profile_version: manifest.calibration_profile_version.clone(),
            backend: manifest.backend.clone(),
            sampling_mode: manifest.sampling_mode.clone(),
            negative_control_status: enum_label(&manifest.negative_control.status)?,
        },
        top_results: result_overview_top_results(candidates)?,
    })
}

pub(crate) fn virus_summary_rows(candidates: &[CandidateCall]) -> Result<Vec<VirusSummaryRow>> {
    let mut rows = candidates
        .iter()
        .cloned()
        .fold(
            BTreeMap::<VirusSummaryKey, VirusSummaryAccumulator>::new(),
            |mut accumulators, candidate| {
                let metadata = CandidateSummaryMetadata::from_candidate(&candidate);
                let key = VirusSummaryKey {
                    aggregation_level: metadata.aggregation_level,
                    entity_id: candidate.accession_or_group.clone(),
                };
                accumulators
                    .entry(key)
                    .or_insert_with(|| VirusSummaryAccumulator::new(metadata.aggregation_level))
                    .push(candidate, metadata);
                accumulators
            },
        )
        .into_values()
        .map(|accumulator| accumulator.into_row())
        .collect::<Result<Vec<_>>>()?;

    sort_virus_summary_rows(&mut rows);
    for (index, row) in rows.iter_mut().enumerate() {
        row.rank = index + 1;
    }
    Ok(rows)
}

pub(crate) fn serialize_sample_run_summary_json(
    summary: &SampleSummary,
    candidates: &[CandidateCall],
    rounds: &[RoundRecord],
) -> Result<Vec<u8>> {
    json_bytes(
        &sample_run_summary(summary, candidates, rounds)?,
        "sample_run_summary",
    )
}

pub(crate) fn serialize_result_overview_json(
    summary: &SampleSummary,
    candidates: &[CandidateCall],
    manifest: &RunManifest,
) -> Result<Vec<u8>> {
    json_bytes(
        &result_overview(summary, candidates, manifest)?,
        "result_overview",
    )
}

pub(crate) fn serialize_virus_summary_tsv(candidates: &[CandidateCall]) -> Result<Vec<u8>> {
    let mut buffer = Vec::new();
    {
        let mut writer = WriterBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_writer(&mut buffer);
        writer.write_record(VIRUS_SUMMARY_HEADER).map_err(|error| {
            RvScreenError::validation(
                "report.virus_summary.tsv",
                format!("failed to serialize virus summary header: {error}"),
            )
        })?;

        for row in virus_summary_rows(candidates)? {
            writer.serialize(row).map_err(|error| {
                RvScreenError::validation(
                    "report.virus_summary.tsv",
                    format!("failed to serialize virus summary row: {error}"),
                )
            })?;
        }

        writer.flush().map_err(|error| {
            RvScreenError::validation("report.virus_summary.tsv", error.to_string())
        })?;
    }
    Ok(buffer)
}

fn result_overview_top_results(
    candidates: &[CandidateCall],
) -> Result<Vec<ResultOverviewTopResult>> {
    ranked_candidates(candidates)
        .into_iter()
        .enumerate()
        .map(|(index, candidate)| {
            Ok(ResultOverviewTopResult {
                rank: index + 1,
                entity_id: candidate.accession_or_group.clone(),
                entity_name: candidate.virus_name.clone(),
                aggregation_level: "accession",
                decision: enum_label(&candidate.decision)?,
                evidence_strength: enum_label(&candidate.evidence_strength)?,
                accepted_fragments: candidate.accepted_fragments,
                unique_fraction: candidate.unique_fraction,
                breadth: candidate.breadth,
                background_ratio: candidate.background_ratio,
            })
        })
        .collect()
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct VirusSummaryKey {
    aggregation_level: SummaryAggregationLevel,
    entity_id: String,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum SummaryAggregationLevel {
    Accession,
    VirusGroup,
}

impl SummaryAggregationLevel {
    fn from_label(label: &str) -> Self {
        match label {
            "virus_group" => Self::VirusGroup,
            _ => Self::Accession,
        }
    }

    fn label(self) -> &'static str {
        match self {
            Self::Accession => "accession",
            Self::VirusGroup => "virus_group",
        }
    }

    fn entity_type(self) -> &'static str {
        self.label()
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct CandidateSummaryMetadata {
    aggregation_level: SummaryAggregationLevel,
    supporting_accession: Option<String>,
    supporting_accession_name: Option<String>,
}

impl CandidateSummaryMetadata {
    fn from_candidate(candidate: &CandidateCall) -> Self {
        let mut metadata = Self {
            aggregation_level: SummaryAggregationLevel::Accession,
            supporting_accession: None,
            supporting_accession_name: None,
        };

        for reason in &candidate.decision_reasons {
            if let Some((key, value)) = reason.split_once('=') {
                match key {
                    "aggregation_level" => {
                        metadata.aggregation_level = SummaryAggregationLevel::from_label(value)
                    }
                    "supporting_accession" if !value.is_empty() => {
                        metadata.supporting_accession = Some(value.to_owned())
                    }
                    "supporting_accession_name" if !value.is_empty() => {
                        metadata.supporting_accession_name = Some(value.to_owned())
                    }
                    _ => {}
                }
            }
        }

        metadata
    }
}

#[derive(Debug, Clone)]
struct CandidateSupportRow {
    candidate: CandidateCall,
    metadata: CandidateSummaryMetadata,
}

#[derive(Debug, Clone)]
struct VirusSummaryAccumulator {
    aggregation_level: SummaryAggregationLevel,
    support_rows: BTreeMap<String, CandidateSupportRow>,
    supporting_accessions: BTreeSet<String>,
    decision_reasons: BTreeSet<String>,
}

impl VirusSummaryAccumulator {
    fn new(aggregation_level: SummaryAggregationLevel) -> Self {
        Self {
            aggregation_level,
            support_rows: BTreeMap::new(),
            supporting_accessions: BTreeSet::new(),
            decision_reasons: BTreeSet::new(),
        }
    }

    fn push(&mut self, candidate: CandidateCall, metadata: CandidateSummaryMetadata) {
        self.decision_reasons
            .extend(candidate.decision_reasons.iter().cloned());

        if let Some(supporting_accession) = &metadata.supporting_accession {
            self.supporting_accessions
                .insert(supporting_accession.clone());
        }

        let support_key = match self.aggregation_level {
            SummaryAggregationLevel::Accession => candidate.accession_or_group.clone(),
            SummaryAggregationLevel::VirusGroup => metadata
                .supporting_accession
                .clone()
                .unwrap_or_else(|| "<missing-supporting-accession>".to_owned()),
        };
        let candidate_support_row = CandidateSupportRow {
            candidate,
            metadata,
        };

        match self.support_rows.get_mut(&support_key) {
            Some(existing)
                if is_stronger_candidate(&candidate_support_row.candidate, &existing.candidate) =>
            {
                *existing = candidate_support_row;
            }
            Some(_) => {}
            None => {
                self.support_rows.insert(support_key, candidate_support_row);
            }
        }
    }

    fn into_row(self) -> Result<VirusSummaryRow> {
        let selected_rows = self.support_rows.into_values().collect::<Vec<_>>();
        let representative = selected_rows
            .iter()
            .map(|row| &row.candidate)
            .max_by(|left, right| compare_candidate_strength(left, right))
            .expect("virus summary accumulator should contain at least one support row");
        let top_support = selected_rows
            .iter()
            .filter(|row| row.metadata.supporting_accession.is_some())
            .max_by(|left, right| compare_candidate_strength(&left.candidate, &right.candidate));

        let top_accession = match self.aggregation_level {
            SummaryAggregationLevel::Accession => representative.accession_or_group.clone(),
            SummaryAggregationLevel::VirusGroup => top_support
                .and_then(|row| row.metadata.supporting_accession.clone())
                .unwrap_or_default(),
        };
        let top_accession_name = match self.aggregation_level {
            SummaryAggregationLevel::Accession => representative.virus_name.clone(),
            SummaryAggregationLevel::VirusGroup => top_support
                .and_then(|row| row.metadata.supporting_accession_name.clone())
                .unwrap_or_default(),
        };

        Ok(VirusSummaryRow {
            rank: 0,
            entity_id: representative.accession_or_group.clone(),
            entity_name: representative.virus_name.clone(),
            entity_type: self.aggregation_level.entity_type(),
            taxid: representative.taxid,
            aggregation_level: self.aggregation_level.label(),
            decision: enum_label(&representative.decision)?,
            evidence_strength: enum_label(&representative.evidence_strength)?,
            accepted_fragments: selected_rows
                .iter()
                .map(|row| row.candidate.accepted_fragments)
                .sum(),
            nonoverlap_fragments: selected_rows
                .iter()
                .map(|row| row.candidate.nonoverlap_fragments)
                .sum(),
            ambiguous_fragments: selected_rows
                .iter()
                .map(|row| row.candidate.ambiguous_fragments)
                .sum(),
            raw_fraction: selected_rows
                .iter()
                .map(|row| row.candidate.raw_fraction)
                .sum(),
            unique_fraction: selected_rows
                .iter()
                .map(|row| row.candidate.unique_fraction)
                .sum(),
            fraction_ci_95_low: selected_rows
                .iter()
                .map(|row| row.candidate.fraction_ci_95[0])
                .fold(f64::INFINITY, f64::min),
            fraction_ci_95_high: selected_rows
                .iter()
                .map(|row| row.candidate.fraction_ci_95[1])
                .fold(f64::NEG_INFINITY, f64::max),
            clopper_pearson_upper: selected_rows
                .iter()
                .map(|row| row.candidate.clopper_pearson_upper)
                .fold(0.0, f64::max),
            breadth: selected_rows
                .iter()
                .map(|row| row.candidate.breadth)
                .fold(0.0, f64::max),
            background_ratio: selected_rows
                .iter()
                .map(|row| row.candidate.background_ratio)
                .fold(0.0, f64::max),
            top_accession,
            top_accession_name,
            supporting_accession_count: match self.aggregation_level {
                SummaryAggregationLevel::Accession => 1,
                SummaryAggregationLevel::VirusGroup => self.supporting_accessions.len() as u64,
            },
            decision_reasons: json_string(
                &self.decision_reasons.into_iter().collect::<Vec<String>>(),
            )?,
        })
    }
}

fn ranked_candidates(candidates: &[CandidateCall]) -> Vec<&CandidateCall> {
    let mut ranked = candidates.iter().collect::<Vec<_>>();
    ranked.sort_by(|left, right| {
        decision_rank(&right.decision)
            .cmp(&decision_rank(&left.decision))
            .then_with(|| {
                evidence_rank(&right.evidence_strength).cmp(&evidence_rank(&left.evidence_strength))
            })
            .then_with(|| right.unique_fraction.total_cmp(&left.unique_fraction))
            .then_with(|| right.breadth.total_cmp(&left.breadth))
            .then_with(|| right.accepted_fragments.cmp(&left.accepted_fragments))
            .then_with(|| left.accession_or_group.cmp(&right.accession_or_group))
    });
    ranked
}

fn sort_virus_summary_rows(rows: &mut [VirusSummaryRow]) {
    rows.sort_by(|left, right| {
        row_decision_rank(&right.decision)
            .cmp(&row_decision_rank(&left.decision))
            .then_with(|| {
                row_evidence_rank(&right.evidence_strength)
                    .cmp(&row_evidence_rank(&left.evidence_strength))
            })
            .then_with(|| right.unique_fraction.total_cmp(&left.unique_fraction))
            .then_with(|| right.breadth.total_cmp(&left.breadth))
            .then_with(|| right.accepted_fragments.cmp(&left.accepted_fragments))
            .then_with(|| left.entity_id.cmp(&right.entity_id))
    });
}

fn is_stronger_candidate(left: &CandidateCall, right: &CandidateCall) -> bool {
    compare_candidate_strength(left, right).is_gt()
}

fn compare_candidate_strength(left: &CandidateCall, right: &CandidateCall) -> std::cmp::Ordering {
    decision_rank(&left.decision)
        .cmp(&decision_rank(&right.decision))
        .then_with(|| {
            evidence_rank(&left.evidence_strength).cmp(&evidence_rank(&right.evidence_strength))
        })
        .then_with(|| left.unique_fraction.total_cmp(&right.unique_fraction))
        .then_with(|| left.breadth.total_cmp(&right.breadth))
        .then_with(|| left.accepted_fragments.cmp(&right.accepted_fragments))
        .then_with(|| right.accession_or_group.cmp(&left.accession_or_group))
}

fn row_decision_rank(decision: &str) -> u8 {
    match decision {
        "positive" => 3,
        "indeterminate" => 2,
        "negative" => 1,
        _ => 0,
    }
}

fn row_evidence_rank(evidence_strength: &str) -> u8 {
    match evidence_strength {
        "high" => 3,
        "medium" => 2,
        "low" => 1,
        _ => 0,
    }
}

fn decision_count(candidates: &[CandidateCall], decision: DecisionStatus) -> usize {
    candidates
        .iter()
        .filter(|candidate| candidate.decision == decision)
        .count()
}

fn decision_rank(decision: &DecisionStatus) -> u8 {
    match decision {
        DecisionStatus::Positive => 3,
        DecisionStatus::Indeterminate => 2,
        DecisionStatus::Negative => 1,
    }
}

fn evidence_rank(evidence_strength: &EvidenceStrength) -> u8 {
    match evidence_strength {
        EvidenceStrength::High => 3,
        EvidenceStrength::Medium => 2,
        EvidenceStrength::Low => 1,
    }
}

fn json_bytes(value: &impl Serialize, artifact_name: &str) -> Result<Vec<u8>> {
    serde_json::to_vec_pretty(value).map_err(|error| {
        RvScreenError::validation(
            format!("report.{artifact_name}.json"),
            format!("failed to serialize {artifact_name}: {error}"),
        )
    })
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::{NegativeControlManifest, NegativeControlStatus, ReleaseStatus, StopReason};
    use csv::ReaderBuilder;

    #[test]
    fn sample_run_summary_json_has_exact_keys_and_counts() {
        let bytes = serialize_sample_run_summary_json(
            &sample_summary_fixture(),
            &candidate_fixtures(),
            &round_fixtures(),
        )
        .expect("sample run summary should serialize");

        assert_eq!(
            String::from_utf8(bytes).expect("summary JSON should be UTF-8"),
            r#"{
  "schema_version": "rvscreen.sample_run_summary.v1",
  "sample_id": "S001",
  "decision_status": "positive",
  "release_status": "final",
  "stop_reason": "positive_boundary_crossed",
  "input_fragments": 24837219,
  "qc_passing_fragments": 24100988,
  "sampled_fragments": 400000,
  "rounds_run": 4,
  "positive_entity_count": 2,
  "indeterminate_entity_count": 1,
  "negative_entity_count": 1
}"#
        );
    }

    #[test]
    fn result_overview_json_has_exact_top_level_fields_and_empty_defaults() {
        let sample_summary_bytes =
            serialize_sample_run_summary_json(&empty_sample_summary_fixture(), &[], &[])
                .expect("empty sample run summary should serialize");
        assert_eq!(
            String::from_utf8(sample_summary_bytes).expect("summary JSON should be UTF-8"),
            r#"{
  "schema_version": "rvscreen.sample_run_summary.v1",
  "sample_id": "S-empty",
  "decision_status": "negative",
  "release_status": "provisional",
  "stop_reason": "negative_boundary_confirmed",
  "input_fragments": 0,
  "qc_passing_fragments": 0,
  "sampled_fragments": 0,
  "rounds_run": 0,
  "positive_entity_count": 0,
  "indeterminate_entity_count": 0,
  "negative_entity_count": 0
}"#
        );

        let bytes = serialize_result_overview_json(
            &empty_sample_summary_fixture(),
            &[],
            &run_manifest_fixture(),
        )
        .expect("result overview should serialize");

        assert_eq!(
            String::from_utf8(bytes).expect("result overview JSON should be UTF-8"),
            r#"{
  "schema_version": "rvscreen.result_overview.v1",
  "sample": {
    "sample_id": "S-empty",
    "decision_status": "negative",
    "release_status": "provisional"
  },
  "run": {
    "reference_bundle_version": "rvscreen_ref_2026.04.20-r1",
    "calibration_profile_version": "rvscreen_calib_2026.04.20-r1",
    "backend": "minimap2",
    "sampling_mode": "representative",
    "negative_control_status": "pass"
  },
  "top_results": []
}"#
        );
    }

    #[test]
    fn virus_summary_tsv_has_exact_header_and_empty_behavior() {
        let bytes = serialize_virus_summary_tsv(&[]).expect("empty virus summary should serialize");
        assert_eq!(
            String::from_utf8(bytes).expect("virus summary TSV should be UTF-8"),
            concat!(
                "rank\tentity_id\tentity_name\tentity_type\ttaxid\taggregation_level\tdecision\t",
                "evidence_strength\taccepted_fragments\tnonoverlap_fragments\tambiguous_fragments\t",
                "raw_fraction\tunique_fraction\tfraction_ci_95_low\tfraction_ci_95_high\t",
                "clopper_pearson_upper\tbreadth\tbackground_ratio\ttop_accession\t",
                "top_accession_name\tsupporting_accession_count\tdecision_reasons\n"
            )
        );
    }

    #[test]
    fn virus_summary_rows_use_phase_one_defaults_and_sort_order() {
        let bytes = serialize_virus_summary_tsv(&candidate_fixtures())
            .expect("virus summary should serialize");
        let mut reader = ReaderBuilder::new()
            .delimiter(b'\t')
            .from_reader(bytes.as_slice());
        let headers = reader.headers().expect("headers should parse").clone();
        assert_eq!(headers.iter().collect::<Vec<_>>(), VIRUS_SUMMARY_HEADER);

        let rows = reader
            .records()
            .collect::<std::result::Result<Vec<_>, _>>()
            .expect("rows should parse");
        let entity_ids = rows
            .iter()
            .map(|row| row.get(1).expect("entity_id column should exist"))
            .collect::<Vec<_>>();
        assert_eq!(
            entity_ids,
            vec![
                "positive-high-a",
                "positive-high-b",
                "indeterminate-mid",
                "negative-low"
            ]
        );

        let first = &rows[0];
        assert_eq!(first.get(0), Some("1"));
        assert_eq!(first.get(3), Some("accession"));
        assert_eq!(first.get(5), Some("accession"));
        assert_eq!(first.get(13), Some("0.00003"));
        assert_eq!(first.get(14), Some("0.00008"));
        assert_eq!(first.get(18), Some("positive-high-a"));
        assert_eq!(first.get(19), Some("Positive high A"));
        assert_eq!(first.get(20), Some("1"));
        assert_eq!(
            first.get(21),
            Some(r#"["unique_fraction_above_theta_pos"]"#)
        );
    }

    #[test]
    fn virus_summary_aggregates_multiple_accessions_for_same_virus_group() {
        let rows = virus_summary_rows(&[
            group_candidate_fixture(GroupCandidateFixtureSpec {
                group: "herpesviridae",
                group_name: "Herpesviridae",
                supporting_accession: Some("NC_001806"),
                supporting_accession_name: Some("Human herpesvirus 1"),
                accepted_fragments: 4,
                nonoverlap_fragments: 2,
                ambiguous_fragments: 1,
                unique_fraction: 0.04,
                breadth: 0.20,
                reason: "support_source=hsv1",
            }),
            group_candidate_fixture(GroupCandidateFixtureSpec {
                group: "herpesviridae",
                group_name: "Herpesviridae",
                supporting_accession: Some("NC_007605"),
                supporting_accession_name: Some("Human herpesvirus 4"),
                accepted_fragments: 7,
                nonoverlap_fragments: 3,
                ambiguous_fragments: 2,
                unique_fraction: 0.07,
                breadth: 0.35,
                reason: "support_source=ebv",
            }),
        ])
        .expect("virus summary rows should aggregate group support");

        assert_eq!(rows.len(), 1);
        let row = &rows[0];
        assert_eq!(row.entity_id, "herpesviridae");
        assert_eq!(row.entity_name, "Herpesviridae");
        assert_eq!(row.entity_type, "virus_group");
        assert_eq!(row.aggregation_level, "virus_group");
        assert_eq!(row.accepted_fragments, 11);
        assert_eq!(row.nonoverlap_fragments, 5);
        assert_eq!(row.ambiguous_fragments, 3);
        assert_eq!(row.top_accession, "NC_007605");
        assert_eq!(row.top_accession_name, "Human herpesvirus 4");
        assert_eq!(row.supporting_accession_count, 2);
        assert_decision_reasons(
            &row.decision_reasons,
            &[
                "aggregation_level=virus_group",
                "support_source=ebv",
                "support_source=hsv1",
                "supporting_accession=NC_001806",
                "supporting_accession=NC_007605",
                "supporting_accession_name=Human herpesvirus 1",
                "supporting_accession_name=Human herpesvirus 4",
            ],
        );
    }

    #[test]
    fn virus_summary_counts_ambiguous_support_without_unique_inflation() {
        let rows = virus_summary_rows(&[
            group_candidate_fixture(GroupCandidateFixtureSpec {
                group: "polyomavirus-group",
                group_name: "Polyomavirus group",
                supporting_accession: Some("NC_polyoma_A"),
                supporting_accession_name: Some("Polyomavirus A"),
                accepted_fragments: 8,
                nonoverlap_fragments: 4,
                ambiguous_fragments: 1,
                unique_fraction: 0.08,
                breadth: 0.30,
                reason: "support_class=unique",
            }),
            group_candidate_fixture(GroupCandidateFixtureSpec {
                group: "polyomavirus-group",
                group_name: "Polyomavirus group",
                supporting_accession: Some("NC_polyoma_B"),
                supporting_accession_name: Some("Polyomavirus B"),
                accepted_fragments: 0,
                nonoverlap_fragments: 0,
                ambiguous_fragments: 6,
                unique_fraction: 0.0,
                breadth: 0.0,
                reason: "support_class=ambiguous_multimap",
            }),
        ])
        .expect("ambiguous support should aggregate deterministically");

        assert_eq!(rows.len(), 1);
        let row = &rows[0];
        assert_eq!(row.accepted_fragments, 8);
        assert_eq!(row.nonoverlap_fragments, 4);
        assert_eq!(row.ambiguous_fragments, 7);
        assert_eq!(row.unique_fraction, 0.08);
        assert_eq!(row.top_accession, "NC_polyoma_A");
        assert_eq!(row.supporting_accession_count, 2);
    }

    #[test]
    fn virus_summary_deduplicates_group_support_units_before_summing_counts() {
        let rows = virus_summary_rows(&[
            group_candidate_fixture(GroupCandidateFixtureSpec {
                group: "adenovirus-group",
                group_name: "Adenovirus group",
                supporting_accession: Some("NC_adeno_A"),
                supporting_accession_name: Some("Adenovirus A"),
                accepted_fragments: 9,
                nonoverlap_fragments: 5,
                ambiguous_fragments: 2,
                unique_fraction: 0.09,
                breadth: 0.40,
                reason: "support_rank=strongest",
            }),
            group_candidate_fixture(GroupCandidateFixtureSpec {
                group: "adenovirus-group",
                group_name: "Adenovirus group",
                supporting_accession: Some("NC_adeno_A"),
                supporting_accession_name: Some("Adenovirus A duplicate"),
                accepted_fragments: 6,
                nonoverlap_fragments: 3,
                ambiguous_fragments: 5,
                unique_fraction: 0.06,
                breadth: 0.20,
                reason: "support_rank=duplicate_lower_score",
            }),
        ])
        .expect("duplicate support rows should not double-count");

        assert_eq!(rows.len(), 1);
        let row = &rows[0];
        assert_eq!(row.accepted_fragments, 9);
        assert_eq!(row.nonoverlap_fragments, 5);
        assert_eq!(row.ambiguous_fragments, 2);
        assert_eq!(row.top_accession, "NC_adeno_A");
        assert_eq!(row.top_accession_name, "Adenovirus A");
        assert_eq!(row.supporting_accession_count, 1);
        assert_decision_reasons(
            &row.decision_reasons,
            &[
                "aggregation_level=virus_group",
                "support_rank=duplicate_lower_score",
                "support_rank=strongest",
                "supporting_accession=NC_adeno_A",
                "supporting_accession_name=Adenovirus A",
                "supporting_accession_name=Adenovirus A duplicate",
            ],
        );
    }

    #[test]
    fn virus_summary_group_rows_with_missing_provenance_use_empty_traceability_defaults() {
        let rows = virus_summary_rows(&[group_candidate_fixture(GroupCandidateFixtureSpec {
            group: "unresolved-virus-group",
            group_name: "Unresolved virus group",
            supporting_accession: None,
            supporting_accession_name: None,
            accepted_fragments: 5,
            nonoverlap_fragments: 2,
            ambiguous_fragments: 3,
            unique_fraction: 0.05,
            breadth: 0.10,
            reason: "supporting_accession_unavailable=true",
        })])
        .expect("missing provenance defaults should serialize");

        assert_eq!(rows.len(), 1);
        let row = &rows[0];
        assert_eq!(row.entity_type, "virus_group");
        assert_eq!(row.aggregation_level, "virus_group");
        assert_eq!(row.top_accession, "");
        assert_eq!(row.top_accession_name, "");
        assert_eq!(row.supporting_accession_count, 0);
        assert_decision_reasons(
            &row.decision_reasons,
            &[
                "aggregation_level=virus_group",
                "supporting_accession_unavailable=true",
            ],
        );
    }

    #[test]
    fn virus_summary_aggregates_segmented_accession_set_without_double_counting_support() {
        let rows = virus_summary_rows(&[
            group_candidate_fixture(GroupCandidateFixtureSpec {
                group: "segmented-virus-set",
                group_name: "Segmented virus set",
                supporting_accession: Some("NC_segment_L"),
                supporting_accession_name: Some("Segment L"),
                accepted_fragments: 5,
                nonoverlap_fragments: 3,
                ambiguous_fragments: 1,
                unique_fraction: 0.05,
                breadth: 0.40,
                reason: "segment=L",
            }),
            group_candidate_fixture(GroupCandidateFixtureSpec {
                group: "segmented-virus-set",
                group_name: "Segmented virus set",
                supporting_accession: Some("NC_segment_S"),
                supporting_accession_name: Some("Segment S"),
                accepted_fragments: 4,
                nonoverlap_fragments: 2,
                ambiguous_fragments: 2,
                unique_fraction: 0.04,
                breadth: 0.30,
                reason: "segment=S",
            }),
            group_candidate_fixture(GroupCandidateFixtureSpec {
                group: "segmented-virus-set",
                group_name: "Segmented virus set",
                supporting_accession: Some("NC_segment_L"),
                supporting_accession_name: Some("Segment L lower support"),
                accepted_fragments: 2,
                nonoverlap_fragments: 1,
                ambiguous_fragments: 9,
                unique_fraction: 0.02,
                breadth: 0.10,
                reason: "segment=L_duplicate_lower_support",
            }),
        ])
        .expect("segmented support should aggregate");

        assert_eq!(rows.len(), 1);
        let row = &rows[0];
        assert_eq!(row.entity_id, "segmented-virus-set");
        assert_eq!(row.aggregation_level, "virus_group");
        assert_eq!(row.accepted_fragments, 9);
        assert_eq!(row.nonoverlap_fragments, 5);
        assert_eq!(row.ambiguous_fragments, 3);
        assert!((row.unique_fraction - 0.09).abs() < f64::EPSILON);
        assert_eq!(row.top_accession, "NC_segment_L");
        assert_eq!(row.top_accession_name, "Segment L");
        assert_eq!(row.supporting_accession_count, 2);
        assert_decision_reasons(
            &row.decision_reasons,
            &[
                "aggregation_level=virus_group",
                "segment=L",
                "segment=L_duplicate_lower_support",
                "segment=S",
                "supporting_accession=NC_segment_L",
                "supporting_accession=NC_segment_S",
                "supporting_accession_name=Segment L",
                "supporting_accession_name=Segment L lower support",
                "supporting_accession_name=Segment S",
            ],
        );
    }

    #[test]
    fn summary_serializers_cover_single_borderline_and_control_like_accession_defaults() {
        let candidates = vec![
            candidate_fixture(CandidateFixtureSpec {
                accession: "NC_single_hit.1",
                name: "Single accession hit",
                decision: DecisionStatus::Positive,
                evidence_strength: EvidenceStrength::Low,
                accepted_fragments: 1,
                unique_fraction: 0.00001,
                breadth: 0.0001,
                taxid: 11_111,
            }),
            CandidateCall {
                virus_name: "Borderline low coverage".into(),
                taxid: 22_222,
                accession_or_group: "NC_borderline.1".into(),
                accepted_fragments: 2,
                nonoverlap_fragments: 1,
                raw_fraction: 0.00002,
                unique_fraction: 0.00001,
                fraction_ci_95: [0.0, 0.00004],
                clopper_pearson_upper: 0.00005,
                breadth: 0.0002,
                ambiguous_fragments: 5,
                background_ratio: 1.0,
                decision: DecisionStatus::Indeterminate,
                decision_reasons: vec![
                    "low_coverage_borderline=true".into(),
                    "threshold_margin=epsilon".into(),
                ],
                evidence_strength: EvidenceStrength::Low,
            },
            CandidateCall {
                virus_name: "Control-like environmental hit".into(),
                taxid: 33_333,
                accession_or_group: "NC_control_like.1".into(),
                accepted_fragments: 3,
                nonoverlap_fragments: 1,
                raw_fraction: 0.00003,
                unique_fraction: 0.00001,
                fraction_ci_95: [0.0, 0.00006],
                clopper_pearson_upper: 0.00007,
                breadth: 0.0001,
                ambiguous_fragments: 1,
                background_ratio: 9.0,
                decision: DecisionStatus::Negative,
                decision_reasons: vec![
                    "control_like_hit=true".into(),
                    "contaminant_or_control_specific_type_unavailable=true".into(),
                ],
                evidence_strength: EvidenceStrength::Low,
            },
        ];

        let summary = sample_run_summary(&sample_summary_fixture(), &candidates, &round_fixtures())
            .expect("sample run summary should derive counts");
        assert_eq!(summary.positive_entity_count, 1);
        assert_eq!(summary.indeterminate_entity_count, 1);
        assert_eq!(summary.negative_entity_count, 1);

        let overview = result_overview(
            &sample_summary_fixture(),
            &candidates,
            &run_manifest_fixture(),
        )
        .expect("result overview should serialize edge cases");
        assert_eq!(overview.top_results[0].entity_id, "NC_single_hit.1");
        assert_eq!(overview.top_results[1].entity_id, "NC_borderline.1");
        assert_eq!(overview.top_results[2].entity_id, "NC_control_like.1");

        let rows = virus_summary_rows(&candidates).expect("virus summary rows should serialize");
        let borderline = rows
            .iter()
            .find(|row| row.entity_id == "NC_borderline.1")
            .expect("borderline row should exist");
        assert_eq!(borderline.decision, "indeterminate");
        assert_eq!(borderline.ambiguous_fragments, 5);
        assert_eq!(borderline.unique_fraction, 0.00001);

        let control_like = rows
            .iter()
            .find(|row| row.entity_id == "NC_control_like.1")
            .expect("control-like row should exist");
        assert_eq!(control_like.entity_type, "accession");
        assert_eq!(control_like.decision, "negative");
        assert_eq!(control_like.top_accession, "NC_control_like.1");
        assert_decision_reasons(
            &control_like.decision_reasons,
            &[
                "contaminant_or_control_specific_type_unavailable=true",
                "control_like_hit=true",
            ],
        );
    }

    #[test]
    fn result_overview_top_results_follow_the_same_stable_sort() {
        let bytes = serialize_result_overview_json(
            &sample_summary_fixture(),
            &candidate_fixtures(),
            &run_manifest_fixture(),
        )
        .expect("result overview should serialize");
        let value: serde_json::Value = serde_json::from_slice(&bytes).expect("JSON should parse");
        let top_results = value["top_results"]
            .as_array()
            .expect("top_results should be an array");

        let ids = top_results
            .iter()
            .map(|result| {
                result["entity_id"]
                    .as_str()
                    .expect("entity_id should be string")
            })
            .collect::<Vec<_>>();
        assert_eq!(
            ids,
            vec![
                "positive-high-a",
                "positive-high-b",
                "indeterminate-mid",
                "negative-low"
            ]
        );
        assert_eq!(top_results[0]["rank"], 1);
        assert_eq!(top_results[0]["aggregation_level"], "accession");
    }

    #[test]
    fn serializers_are_byte_stable_across_repeated_calls() {
        let summary = sample_summary_fixture();
        let candidates = candidate_fixtures();
        let rounds = round_fixtures();
        let manifest = run_manifest_fixture();

        assert_eq!(
            serialize_sample_run_summary_json(&summary, &candidates, &rounds)
                .expect("first sample run summary should serialize"),
            serialize_sample_run_summary_json(&summary, &candidates, &rounds)
                .expect("second sample run summary should serialize")
        );
        assert_eq!(
            serialize_result_overview_json(&summary, &candidates, &manifest)
                .expect("first result overview should serialize"),
            serialize_result_overview_json(&summary, &candidates, &manifest)
                .expect("second result overview should serialize")
        );
        assert_eq!(
            serialize_virus_summary_tsv(&candidates).expect("first virus summary should serialize"),
            serialize_virus_summary_tsv(&candidates)
                .expect("second virus summary should serialize")
        );
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
            rounds_run: 4,
            stop_reason: StopReason::PositiveBoundaryCrossed,
            decision_status: DecisionStatus::Positive,
            release_status: ReleaseStatus::Final,
        }
    }

    fn empty_sample_summary_fixture() -> SampleSummary {
        SampleSummary {
            sample_id: "S-empty".into(),
            reference_bundle_version: "rvscreen_ref_2026.04.20-r1".into(),
            calibration_profile_version: "rvscreen_calib_2026.04.20-r1".into(),
            backend: "minimap2".into(),
            seed: 20_260_421,
            input_fragments: 0,
            qc_passing_fragments: 0,
            sampled_fragments: 0,
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
            negative_control: NegativeControlManifest {
                required: false,
                status: NegativeControlStatus::Pass,
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
                sampled_fragments: 200_000,
                accepted_virus: 12,
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
            candidate_fixture(CandidateFixtureSpec {
                accession: "negative-low",
                name: "Negative low",
                decision: DecisionStatus::Negative,
                evidence_strength: EvidenceStrength::Low,
                accepted_fragments: 200,
                unique_fraction: 0.99,
                breadth: 0.99,
                taxid: 900,
            }),
            candidate_fixture(CandidateFixtureSpec {
                accession: "positive-high-a",
                name: "Positive high A",
                decision: DecisionStatus::Positive,
                evidence_strength: EvidenceStrength::High,
                accepted_fragments: 50,
                unique_fraction: 0.40,
                breadth: 0.25,
                taxid: 100,
            }),
            candidate_fixture(CandidateFixtureSpec {
                accession: "indeterminate-mid",
                name: "Indeterminate mid",
                decision: DecisionStatus::Indeterminate,
                evidence_strength: EvidenceStrength::Medium,
                accepted_fragments: 300,
                unique_fraction: 0.80,
                breadth: 0.80,
                taxid: 300,
            }),
            candidate_fixture(CandidateFixtureSpec {
                accession: "positive-high-b",
                name: "Positive high B",
                decision: DecisionStatus::Positive,
                evidence_strength: EvidenceStrength::High,
                accepted_fragments: 40,
                unique_fraction: 0.40,
                breadth: 0.25,
                taxid: 101,
            }),
        ]
    }

    struct CandidateFixtureSpec {
        accession: &'static str,
        name: &'static str,
        decision: DecisionStatus,
        evidence_strength: EvidenceStrength,
        accepted_fragments: u64,
        unique_fraction: f64,
        breadth: f64,
        taxid: u64,
    }

    struct GroupCandidateFixtureSpec {
        group: &'static str,
        group_name: &'static str,
        supporting_accession: Option<&'static str>,
        supporting_accession_name: Option<&'static str>,
        accepted_fragments: u64,
        nonoverlap_fragments: u64,
        ambiguous_fragments: u64,
        unique_fraction: f64,
        breadth: f64,
        reason: &'static str,
    }

    fn group_candidate_fixture(spec: GroupCandidateFixtureSpec) -> CandidateCall {
        let mut decision_reasons = vec!["aggregation_level=virus_group".to_string()];
        if let Some(supporting_accession) = spec.supporting_accession {
            decision_reasons.push(format!("supporting_accession={supporting_accession}"));
        }
        if let Some(supporting_accession_name) = spec.supporting_accession_name {
            decision_reasons.push(format!(
                "supporting_accession_name={supporting_accession_name}"
            ));
        }
        decision_reasons.push(spec.reason.to_string());

        CandidateCall {
            virus_name: spec.group_name.into(),
            taxid: 10_001,
            accession_or_group: spec.group.into(),
            accepted_fragments: spec.accepted_fragments,
            nonoverlap_fragments: spec.nonoverlap_fragments,
            raw_fraction: spec.unique_fraction + 0.01,
            unique_fraction: spec.unique_fraction,
            fraction_ci_95: [0.00003, 0.00008],
            clopper_pearson_upper: 0.0001,
            breadth: spec.breadth,
            ambiguous_fragments: spec.ambiguous_fragments,
            background_ratio: 2.4,
            decision: DecisionStatus::Positive,
            decision_reasons,
            evidence_strength: EvidenceStrength::High,
        }
    }

    fn assert_decision_reasons(actual: &str, expected: &[&str]) {
        let parsed = serde_json::from_str::<Vec<String>>(actual)
            .expect("decision_reasons should remain a JSON array string");
        assert_eq!(parsed, expected);
    }

    fn candidate_fixture(spec: CandidateFixtureSpec) -> CandidateCall {
        CandidateCall {
            virus_name: spec.name.into(),
            taxid: spec.taxid,
            accession_or_group: spec.accession.into(),
            accepted_fragments: spec.accepted_fragments,
            nonoverlap_fragments: spec.accepted_fragments / 2,
            raw_fraction: spec.unique_fraction + 0.01,
            unique_fraction: spec.unique_fraction,
            fraction_ci_95: [0.00003, 0.00008],
            clopper_pearson_upper: 0.0001,
            breadth: spec.breadth,
            ambiguous_fragments: 2,
            background_ratio: 2.4,
            decision: spec.decision,
            decision_reasons: vec!["unique_fraction_above_theta_pos".into()],
            evidence_strength: spec.evidence_strength,
        }
    }
}
