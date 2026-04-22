use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
#[serde(transparent)]
pub struct BundleManifest(pub Vec<ContigEntry>);

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct ContigEntry {
    pub contig: String,
    pub accession: String,
    pub taxid: u64,
    pub virus_name: String,
    pub segment: Option<String>,
    pub group: String,
    pub genome_length: u64,
    pub source_release: String,
    pub source_type: String,
    pub masked_regions: Vec<[u64; 2]>,
}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct BundleToml {
    pub version: String,
    pub created_at: String,
    pub included_layers: Vec<String>,
}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct ProfileToml {
    pub profile_id: String,
    pub status: String,
    pub reference_bundle: String,
    pub backend: String,
    pub preset: String,
    pub seed: u64,
    pub supported_input: Vec<String>,
    pub supported_read_type: Vec<String>,
    pub negative_control_required: bool,
    pub sampling: SamplingConfig,
    pub fragment_rules: FragmentRules,
    pub candidate_rules: CandidateRules,
    pub decision_rules: DecisionRules,
}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct SamplingConfig {
    pub mode: String,
    pub rounds: Vec<u64>,
    pub max_rounds: u64,
}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct FragmentRules {
    pub min_mapq: u32,
    pub min_as_diff: i32,
    pub max_nm: u32,
    pub require_pair_consistency: bool,
}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct CandidateRules {
    pub min_nonoverlap_fragments: u64,
    pub min_breadth: f64,
    pub max_background_ratio: f64,
}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct DecisionRules {
    pub theta_pos: f64,
    pub theta_neg: f64,
    pub allow_indeterminate: bool,
}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
#[serde(rename_all = "snake_case")]
pub enum DecisionStatus {
    Positive,
    Negative,
    Indeterminate,
}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
#[serde(rename_all = "snake_case")]
pub enum ReleaseStatus {
    Final,
    Provisional,
    Blocked,
}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
#[serde(rename_all = "snake_case")]
pub enum FragmentClass {
    Host,
    Virus,
    Ambiguous,
    Unmapped,
}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
#[serde(rename_all = "snake_case")]
pub enum EvidenceStrength {
    Low,
    Medium,
    High,
}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
#[serde(rename_all = "snake_case")]
pub enum StopReason {
    PositiveBoundaryCrossed,
    NegativeBoundaryConfirmed,
    MaxRoundsReached,
}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct SampleSummary {
    pub sample_id: String,
    pub reference_bundle_version: String,
    pub calibration_profile_version: String,
    pub backend: String,
    pub seed: u64,
    pub input_fragments: u64,
    pub qc_passing_fragments: u64,
    pub sampled_fragments: u64,
    pub rounds_run: u64,
    pub stop_reason: StopReason,
    pub decision_status: DecisionStatus,
    pub release_status: ReleaseStatus,
}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct CandidateCall {
    pub virus_name: String,
    pub taxid: u64,
    pub accession_or_group: String,
    pub accepted_fragments: u64,
    pub nonoverlap_fragments: u64,
    pub raw_fraction: f64,
    pub unique_fraction: f64,
    pub fraction_ci_95: [f64; 2],
    pub clopper_pearson_upper: f64,
    pub breadth: f64,
    pub ambiguous_fragments: u64,
    pub background_ratio: f64,
    pub decision: DecisionStatus,
    pub decision_reasons: Vec<String>,
    pub evidence_strength: EvidenceStrength,
}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct RunManifest {
    pub reference_bundle_version: String,
    pub calibration_profile_version: String,
    pub backend: String,
    pub seed: u64,
    pub input_files: Vec<String>,
}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct RoundRecord {
    pub sampled_fragments: u64,
    pub accepted_virus: u64,
    pub decision_status: DecisionStatus,
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeSet;

    #[test]
    fn test_profile_toml_roundtrip() {
        let input = r#"
profile_id = "rvscreen_calib_2026.04.20-r1"
status = "release_candidate"
reference_bundle = "rvscreen_ref_2026.04.20-r1"
backend = "minimap2"
preset = "sr-conservative"
seed = 20260420
supported_input = ["fastq.gz", "ubam", "bam", "cram"]
supported_read_type = ["illumina_pe_shortread"]
negative_control_required = true

[sampling]
mode = "representative"
rounds = [50000, 100000, 200000, 400000, 800000]
max_rounds = 5

[fragment_rules]
min_mapq = 20
min_as_diff = 12
max_nm = 8
require_pair_consistency = true

[candidate_rules]
min_nonoverlap_fragments = 3
min_breadth = 0.001
max_background_ratio = 1.5

[decision_rules]
theta_pos = 0.00005
theta_neg = 0.000005
allow_indeterminate = true
"#;

        let parsed: ProfileToml = toml::from_str(input).expect("profile TOML should parse");
        assert_eq!(parsed.profile_id, "rvscreen_calib_2026.04.20-r1");
        assert_eq!(parsed.fragment_rules.min_mapq, 20);
        assert_eq!(parsed.fragment_rules.min_as_diff, 12);
        assert_eq!(
            parsed.sampling.rounds,
            vec![50_000, 100_000, 200_000, 400_000, 800_000]
        );

        let roundtrip = toml::to_string(&parsed).expect("profile TOML should serialize");
        let reparsed: ProfileToml =
            toml::from_str(&roundtrip).expect("round-tripped profile TOML should parse");
        assert_eq!(reparsed, parsed);
    }

    #[test]
    fn test_sample_summary_json_fields() {
        let summary = sample_summary_fixture();
        let json = serde_json::to_value(&summary).expect("sample summary should serialize");
        let object = json
            .as_object()
            .expect("sample summary JSON should be an object");

        assert_json_keys_eq(
            object.keys().map(String::as_str),
            [
                "sample_id",
                "reference_bundle_version",
                "calibration_profile_version",
                "backend",
                "seed",
                "input_fragments",
                "qc_passing_fragments",
                "sampled_fragments",
                "rounds_run",
                "stop_reason",
                "decision_status",
                "release_status",
            ],
        );
        assert_eq!(object.len(), 12);

        let roundtrip: SampleSummary =
            serde_json::from_value(json).expect("sample summary should deserialize");
        assert_eq!(roundtrip, summary);
    }

    #[test]
    fn test_candidate_call_json() {
        let candidate = candidate_call_fixture();
        let json = serde_json::to_value(&candidate).expect("candidate call should serialize");
        let object = json
            .as_object()
            .expect("candidate call JSON should be an object");

        assert_json_keys_eq(
            object.keys().map(String::as_str),
            [
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
            ],
        );
        assert_eq!(object.len(), 15);

        let roundtrip: CandidateCall =
            serde_json::from_value(json).expect("candidate call should deserialize");
        assert_eq!(roundtrip, candidate);
    }

    #[test]
    fn test_core_enum_serialization() {
        assert_json_roundtrip(DecisionStatus::Positive, "\"positive\"");
        assert_json_roundtrip(DecisionStatus::Negative, "\"negative\"");
        assert_json_roundtrip(DecisionStatus::Indeterminate, "\"indeterminate\"");

        assert_json_roundtrip(ReleaseStatus::Final, "\"final\"");
        assert_json_roundtrip(ReleaseStatus::Provisional, "\"provisional\"");
        assert_json_roundtrip(ReleaseStatus::Blocked, "\"blocked\"");

        assert_json_roundtrip(FragmentClass::Host, "\"host\"");
        assert_json_roundtrip(FragmentClass::Virus, "\"virus\"");
        assert_json_roundtrip(FragmentClass::Ambiguous, "\"ambiguous\"");
        assert_json_roundtrip(FragmentClass::Unmapped, "\"unmapped\"");

        assert_json_roundtrip(EvidenceStrength::Low, "\"low\"");
        assert_json_roundtrip(EvidenceStrength::Medium, "\"medium\"");
        assert_json_roundtrip(EvidenceStrength::High, "\"high\"");

        assert_json_roundtrip(
            StopReason::PositiveBoundaryCrossed,
            "\"positive_boundary_crossed\"",
        );
        assert_json_roundtrip(
            StopReason::NegativeBoundaryConfirmed,
            "\"negative_boundary_confirmed\"",
        );
        assert_json_roundtrip(StopReason::MaxRoundsReached, "\"max_rounds_reached\"");
    }

    #[test]
    fn test_reference_and_report_metadata_roundtrip() {
        let manifest = BundleManifest(vec![ContigEntry {
            contig: "NC_007605.1".into(),
            accession: "NC_007605.1".into(),
            taxid: 10_376,
            virus_name: "Epstein-Barr virus".into(),
            segment: None,
            group: "ebv".into(),
            genome_length: 171_823,
            source_release: "RefSeq_2026_04".into(),
            source_type: "virus_panel".into(),
            masked_regions: Vec::new(),
        }]);
        let manifest_json = serde_json::to_string(&manifest).expect("manifest should serialize");
        let manifest_roundtrip: BundleManifest =
            serde_json::from_str(&manifest_json).expect("manifest should deserialize");
        assert_eq!(manifest_roundtrip, manifest);

        let bundle = BundleToml {
            version: "rvscreen_ref_2026.04.20-r1".into(),
            created_at: "2026-04-20T00:00:00Z".into(),
            included_layers: vec![
                "host_backbone".into(),
                "viral_panel".into(),
                "decoy_panel".into(),
            ],
        };
        let bundle_toml = toml::to_string(&bundle).expect("bundle TOML should serialize");
        let bundle_roundtrip: BundleToml =
            toml::from_str(&bundle_toml).expect("bundle TOML should deserialize");
        assert_eq!(bundle_roundtrip, bundle);

        let run_manifest = RunManifest {
            reference_bundle_version: "rvscreen_ref_2026.04.20-r1".into(),
            calibration_profile_version: "rvscreen_calib_2026.04.20-r1".into(),
            backend: "minimap2".into(),
            seed: 20_260_420,
            input_files: vec!["reads_R1.fastq.gz".into(), "reads_R2.fastq.gz".into()],
        };
        let run_manifest_json =
            serde_json::to_string(&run_manifest).expect("run manifest should serialize");
        let run_manifest_roundtrip: RunManifest =
            serde_json::from_str(&run_manifest_json).expect("run manifest should deserialize");
        assert_eq!(run_manifest_roundtrip, run_manifest);

        let round = RoundRecord {
            sampled_fragments: 400_000,
            accepted_virus: 32,
            decision_status: DecisionStatus::Positive,
        };
        let round_json = serde_json::to_string(&round).expect("round record should serialize");
        let round_roundtrip: RoundRecord =
            serde_json::from_str(&round_json).expect("round record should deserialize");
        assert_eq!(round_roundtrip, round);
    }

    fn assert_json_keys_eq<'a, const N: usize>(
        actual: impl IntoIterator<Item = &'a str>,
        expected: [&str; N],
    ) {
        let actual: BTreeSet<&str> = actual.into_iter().collect();
        let expected: BTreeSet<&str> = expected.into_iter().collect();
        assert_eq!(actual, expected);
    }

    fn assert_json_roundtrip<T>(value: T, expected_json: &str)
    where
        T: Serialize + for<'de> Deserialize<'de> + PartialEq + std::fmt::Debug,
    {
        let serialized = serde_json::to_string(&value).expect("enum should serialize");
        assert_eq!(serialized, expected_json);

        let deserialized: T = serde_json::from_str(expected_json).expect("enum should deserialize");
        assert_eq!(deserialized, value);
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

    fn candidate_call_fixture() -> CandidateCall {
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
                "unique_fraction_above_theta_pos".into(),
                "background_ratio_above_threshold".into(),
            ],
            evidence_strength: EvidenceStrength::High,
        }
    }
}
