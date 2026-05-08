use crate::calibration::{GateStatus, ReleaseGate as BenchmarkGates};
use crate::cli::ScreenMode;
use crate::decision::NegativeControlDecisionInput;
use crate::error::{Result, RvScreenError};
use crate::types::{
    CandidateRules, DecisionRules, EvidenceStrength, FragmentRules, ProfileToml, ReleaseStatus,
    SamplingConfig,
};
use std::fs;
use std::path::Path;

const RELEASE_GATE_JSON: &str = "release_gate.json";
const FAILED_PROFILE_STATUS: &str = "benchmark_failed";

#[derive(Debug, Clone, PartialEq)]
pub struct ReleaseContext<'a> {
    pub sampling_mode: ScreenMode,
    pub negative_control: &'a NegativeControlDecisionInput,
    pub calibration_profile: &'a ProfileToml,
    pub benchmark_gates: Option<&'a BenchmarkGates>,
}

#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct ReleaseGateDecision {
    pub status: ReleaseStatus,
    #[serde(default)]
    pub primary_blockers: Vec<String>,
    #[serde(default)]
    pub secondary_blockers: Vec<String>,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct ReleaseGate;

impl ReleaseGate {
    pub fn evaluate(context: &ReleaseContext) -> ReleaseStatus {
        Self::evaluate_decision(context).status
    }

    pub fn evaluate_decision(context: &ReleaseContext) -> ReleaseGateDecision {
        let mut primary_blockers = Vec::new();
        let mut secondary_blockers = Vec::new();

        if profile_gate_blocked(context) {
            primary_blockers.push("benchmark_failed_profile".to_string());
        }
        if negative_control_gate_blocked(context) {
            primary_blockers.push("failed_negative_control".to_string());
        }
        match context.benchmark_gates {
            None => primary_blockers.push("missing_calibration_release_gate".to_string()),
            Some(gates) => {
                if gate_failed(&gates.backend_gate.status) {
                    primary_blockers.push("backend_gate_failed".to_string());
                }
                if gate_failed(&gates.reference_gate.status) {
                    primary_blockers.push("reference_gate_failed".to_string());
                }
                if gate_failed(&gates.specificity_gate.status) {
                    primary_blockers.push("specificity_gate_failed".to_string());
                }
                if gate_failed(&gates.sensitivity_gate.status) {
                    primary_blockers.push("sensitivity_gate_failed".to_string());
                }
            }
        }
        if context.sampling_mode != ScreenMode::Representative {
            secondary_blockers.push("non_representative_sampling_mode".to_string());
        }

        let status = if !primary_blockers.is_empty() {
            ReleaseStatus::Blocked
        } else if !secondary_blockers.is_empty() {
            ReleaseStatus::Provisional
        } else {
            ReleaseStatus::Final
        };

        ReleaseGateDecision {
            status,
            primary_blockers,
            secondary_blockers,
        }
    }
}

pub fn expected_release_status(
    sampling_mode: &str,
    negative_control_required: bool,
    negative_control: &NegativeControlDecisionInput,
    benchmark_gates: Option<&BenchmarkGates>,
    profile_status: &str,
) -> Result<ReleaseStatus> {
    Ok(expected_release_gate_decision(
        sampling_mode,
        negative_control_required,
        negative_control,
        benchmark_gates,
        profile_status,
    )?
    .status)
}

pub fn expected_release_gate_decision(
    sampling_mode: &str,
    negative_control_required: bool,
    negative_control: &NegativeControlDecisionInput,
    benchmark_gates: Option<&BenchmarkGates>,
    profile_status: &str,
) -> Result<ReleaseGateDecision> {
    let sampling_mode = parse_sampling_mode(sampling_mode)?;
    let calibration_profile =
        synthetic_profile(negative_control_required, profile_status, sampling_mode);
    Ok(ReleaseGate::evaluate_decision(&ReleaseContext {
        sampling_mode,
        negative_control,
        calibration_profile: &calibration_profile,
        benchmark_gates,
    }))
}

pub fn load_benchmark_gates(profile_dir: impl AsRef<Path>) -> Result<Option<BenchmarkGates>> {
    let path = profile_dir.as_ref().join(RELEASE_GATE_JSON);
    if !path.exists() {
        return Ok(None);
    }

    let body = fs::read_to_string(&path).map_err(|source| RvScreenError::io(&path, source))?;
    let gates = serde_json::from_str(&body).map_err(|error| {
        RvScreenError::parse(
            &path,
            error.line() as u64,
            format!("column {}: {error}", error.column()),
        )
    })?;

    Ok(Some(gates))
}

fn negative_control_gate_blocked(context: &ReleaseContext) -> bool {
    match context.negative_control {
        NegativeControlDecisionInput::Missing => false,
        NegativeControlDecisionInput::Failed(_) => true,
        NegativeControlDecisionInput::Passed { .. } => false,
    }
}

fn profile_gate_blocked(context: &ReleaseContext) -> bool {
    context
        .calibration_profile
        .status
        .trim()
        .eq_ignore_ascii_case(FAILED_PROFILE_STATUS)
}

fn gate_failed(status: &GateStatus) -> bool {
    *status == GateStatus::Fail
}

fn parse_sampling_mode(mode: &str) -> Result<ScreenMode> {
    match mode {
        "representative" => Ok(ScreenMode::Representative),
        "streaming" => Ok(ScreenMode::Streaming),
        other => Err(RvScreenError::validation(
            "run_manifest.sampling_mode",
            format!("unsupported sampling mode `{other}` in report bundle"),
        )),
    }
}

fn synthetic_profile(
    negative_control_required: bool,
    status: &str,
    sampling_mode: ScreenMode,
) -> ProfileToml {
    ProfileToml {
        profile_id: "audit-profile".to_string(),
        status: status.to_string(),
        reference_bundle: "audit-reference".to_string(),
        backend: "minimap2".to_string(),
        preset: "sr-conservative".to_string(),
        seed: 0,
        supported_input: vec!["fastq".to_string()],
        supported_read_type: vec!["illumina_pe_shortread".to_string()],
        negative_control_required,
        sampling: SamplingConfig {
            mode: match sampling_mode {
                ScreenMode::Representative => "representative".to_string(),
                ScreenMode::Streaming => "streaming".to_string(),
            },
            rounds: vec![1],
            round_mode: None,
            round_proportions: None,
            max_rounds: 1,
        },
        fragment_rules: FragmentRules {
            min_mapq: 0,
            min_as_diff: 0,
            max_nm: 0,
            require_pair_consistency: true,
        },
        candidate_rules: CandidateRules {
            min_nonoverlap_fragments: 1,
            min_breadth: 0.0,
            max_background_ratio: 0.0,
            theta_pos_absolute: 1,
            max_ambiguous_fraction_for_positive: 0.20,
            min_positive_evidence_strength: EvidenceStrength::Moderate,
            weak_positive_enabled: true,
        },
        decision_rules: DecisionRules {
            theta_pos: 0.0,
            theta_neg: 0.0,
            allow_indeterminate: true,
            theta_neg_fixed: 0.00001,
            positive_alpha_global: 0.05,
            positive_statistic_method: "exact_binomial_survival_bonferroni".to_string(),
            negative_cp_lod_max: 0.00001,
            noise_floor_fraction: 0.000001,
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::calibration::{BackendGate, ReferenceGate, SensitivityGate, SpecificityGate};
    use crate::decision::{BackgroundComparator, NegativeControlResult};
    use crate::types::{CandidateRules, DecisionRules, FragmentRules, ProfileToml, SamplingConfig};
    use tempfile::tempdir;

    #[test]
    fn representative_passed_negative_control_and_passing_gates_are_final() {
        let benchmark_gates = benchmark_gates(
            GateStatus::Pass,
            GateStatus::Pass,
            GateStatus::Pass,
            GateStatus::Pass,
        );
        let negative_control = passing_negative_control();

        let status = ReleaseGate::evaluate(&ReleaseContext {
            sampling_mode: ScreenMode::Representative,
            negative_control: &negative_control,
            calibration_profile: &profile(true, "release_candidate"),
            benchmark_gates: Some(&benchmark_gates),
        });

        assert_eq!(status, ReleaseStatus::Final);
    }

    #[test]
    fn missing_negative_control_even_when_profile_requires_is_metadata_only() {
        let benchmark_gates = benchmark_gates(
            GateStatus::Pass,
            GateStatus::Pass,
            GateStatus::Pass,
            GateStatus::Pass,
        );

        let status = ReleaseGate::evaluate(&ReleaseContext {
            sampling_mode: ScreenMode::Representative,
            negative_control: &NegativeControlDecisionInput::Missing,
            calibration_profile: &profile(true, "release_candidate"),
            benchmark_gates: Some(&benchmark_gates),
        });

        assert_eq!(status, ReleaseStatus::Final);
    }

    #[test]
    fn streaming_mode_is_provisional_even_with_passing_gates() {
        let benchmark_gates = benchmark_gates(
            GateStatus::Pass,
            GateStatus::Pass,
            GateStatus::Pass,
            GateStatus::Pass,
        );
        let negative_control = passing_negative_control();

        let status = ReleaseGate::evaluate(&ReleaseContext {
            sampling_mode: ScreenMode::Streaming,
            negative_control: &negative_control,
            calibration_profile: &profile(true, "release_candidate"),
            benchmark_gates: Some(&benchmark_gates),
        });

        assert_eq!(status, ReleaseStatus::Provisional);
    }

    #[test]
    fn optional_missing_negative_control_with_passing_gates_is_final() {
        let benchmark_gates = benchmark_gates(
            GateStatus::Pass,
            GateStatus::Pass,
            GateStatus::Pass,
            GateStatus::Pass,
        );

        let status = ReleaseGate::evaluate(&ReleaseContext {
            sampling_mode: ScreenMode::Representative,
            negative_control: &NegativeControlDecisionInput::Missing,
            calibration_profile: &profile(false, "release_candidate"),
            benchmark_gates: Some(&benchmark_gates),
        });

        assert_eq!(status, ReleaseStatus::Final);
    }

    #[test]
    fn failed_negative_control_is_blocked() {
        let benchmark_gates = benchmark_gates(
            GateStatus::Pass,
            GateStatus::Pass,
            GateStatus::Pass,
            GateStatus::Pass,
        );
        let negative_control = NegativeControlDecisionInput::Failed(NegativeControlResult {
            control_id: "neg-001".to_string(),
            control_status: "fail".to_string(),
            candidates: Vec::new(),
        });

        let status = ReleaseGate::evaluate(&ReleaseContext {
            sampling_mode: ScreenMode::Representative,
            negative_control: &negative_control,
            calibration_profile: &profile(false, "release_candidate"),
            benchmark_gates: Some(&benchmark_gates),
        });

        assert_eq!(status, ReleaseStatus::Blocked);
    }

    #[test]
    fn failed_backend_gate_is_blocked() {
        let benchmark_gates = benchmark_gates(
            GateStatus::Fail,
            GateStatus::Pass,
            GateStatus::Pass,
            GateStatus::Pass,
        );
        let negative_control = passing_negative_control();

        let status = ReleaseGate::evaluate(&ReleaseContext {
            sampling_mode: ScreenMode::Representative,
            negative_control: &negative_control,
            calibration_profile: &profile(true, "release_candidate"),
            benchmark_gates: Some(&benchmark_gates),
        });

        assert_eq!(status, ReleaseStatus::Blocked);
    }

    #[test]
    fn failed_reference_gate_is_blocked() {
        let benchmark_gates = benchmark_gates(
            GateStatus::Pass,
            GateStatus::Fail,
            GateStatus::Pass,
            GateStatus::Pass,
        );
        let negative_control = passing_negative_control();

        let status = ReleaseGate::evaluate(&ReleaseContext {
            sampling_mode: ScreenMode::Representative,
            negative_control: &negative_control,
            calibration_profile: &profile(true, "release_candidate"),
            benchmark_gates: Some(&benchmark_gates),
        });

        assert_eq!(status, ReleaseStatus::Blocked);
    }

    #[test]
    fn failed_statistical_gate_is_blocked() {
        let benchmark_gates = benchmark_gates(
            GateStatus::Pass,
            GateStatus::Pass,
            GateStatus::Fail,
            GateStatus::Pass,
        );
        let negative_control = passing_negative_control();

        let status = ReleaseGate::evaluate(&ReleaseContext {
            sampling_mode: ScreenMode::Representative,
            negative_control: &negative_control,
            calibration_profile: &profile(true, "release_candidate"),
            benchmark_gates: Some(&benchmark_gates),
        });

        assert_eq!(status, ReleaseStatus::Blocked);
    }

    #[test]
    fn failed_sensitivity_gate_is_blocked() {
        let benchmark_gates = benchmark_gates(
            GateStatus::Pass,
            GateStatus::Pass,
            GateStatus::Pass,
            GateStatus::Fail,
        );
        let negative_control = passing_negative_control();

        let status = ReleaseGate::evaluate(&ReleaseContext {
            sampling_mode: ScreenMode::Representative,
            negative_control: &negative_control,
            calibration_profile: &profile(true, "release_candidate"),
            benchmark_gates: Some(&benchmark_gates),
        });

        assert_eq!(status, ReleaseStatus::Blocked);
    }

    #[test]
    fn benchmark_failed_profile_is_blocked() {
        let benchmark_gates = benchmark_gates(
            GateStatus::Pass,
            GateStatus::Pass,
            GateStatus::Pass,
            GateStatus::Pass,
        );
        let negative_control = passing_negative_control();

        let status = ReleaseGate::evaluate(&ReleaseContext {
            sampling_mode: ScreenMode::Representative,
            negative_control: &negative_control,
            calibration_profile: &profile(true, "benchmark_failed"),
            benchmark_gates: Some(&benchmark_gates),
        });

        assert_eq!(status, ReleaseStatus::Blocked);
    }

    #[test]
    fn missing_release_gate_artifact_blocks_final_status() {
        let negative_control = passing_negative_control();

        let status = ReleaseGate::evaluate(&ReleaseContext {
            sampling_mode: ScreenMode::Representative,
            negative_control: &negative_control,
            calibration_profile: &profile(true, "release_candidate"),
            benchmark_gates: None,
        });

        assert_eq!(status, ReleaseStatus::Blocked);
    }

    #[test]
    fn load_benchmark_gates_reads_release_gate_json() {
        let tempdir = tempdir().expect("tempdir should be created");
        let gates = benchmark_gates(
            GateStatus::Pass,
            GateStatus::Pass,
            GateStatus::Pass,
            GateStatus::Pass,
        );
        let artifact_path = tempdir.path().join(RELEASE_GATE_JSON);
        fs::write(
            &artifact_path,
            serde_json::to_string_pretty(&gates).expect("release gate should serialize"),
        )
        .expect("release gate artifact should be written");

        let loaded = load_benchmark_gates(tempdir.path())
            .expect("release gate artifact should load")
            .expect("release gate artifact should be present");

        assert_eq!(loaded, gates);
    }

    fn profile(negative_control_required: bool, status: &str) -> ProfileToml {
        ProfileToml {
            profile_id: "rvscreen_calib_test".to_string(),
            status: status.to_string(),
            reference_bundle: "rvscreen_ref_test".to_string(),
            backend: "minimap2".to_string(),
            preset: "sr-conservative".to_string(),
            seed: 20_260_420,
            supported_input: vec!["fastq".to_string()],
            supported_read_type: vec!["illumina_pe_shortread".to_string()],
            negative_control_required,
            sampling: SamplingConfig {
                mode: "representative".to_string(),
                rounds: vec![50, 100],
                round_mode: None,
                round_proportions: None,
                max_rounds: 2,
            },
            fragment_rules: FragmentRules {
                min_mapq: 0,
                min_as_diff: 0,
                max_nm: 100,
                require_pair_consistency: true,
            },
            candidate_rules: CandidateRules {
                min_nonoverlap_fragments: 1,
                min_breadth: 0.0,
                max_background_ratio: 1.5,
                theta_pos_absolute: 1,
                max_ambiguous_fraction_for_positive: 0.20,
                min_positive_evidence_strength: EvidenceStrength::Moderate,
                weak_positive_enabled: true,
            },
            decision_rules: DecisionRules {
                theta_pos: 0.01,
                theta_neg: 0.0001,
                allow_indeterminate: true,
                theta_neg_fixed: 0.00001,
                positive_alpha_global: 0.05,
                positive_statistic_method: "exact_binomial_survival_bonferroni".to_string(),
                negative_cp_lod_max: 0.00001,
                noise_floor_fraction: 0.000001,
            },
        }
    }

    fn benchmark_gates(
        backend: GateStatus,
        reference: GateStatus,
        specificity: GateStatus,
        sensitivity: GateStatus,
    ) -> BenchmarkGates {
        let specificity_failed = specificity == GateStatus::Fail;
        let sensitivity_passed = sensitivity == GateStatus::Pass;

        BenchmarkGates {
            backend_gate: BackendGate {
                status: backend,
                details: "backend gate".to_string(),
            },
            reference_gate: ReferenceGate {
                status: reference,
                details: "reference gate".to_string(),
            },
            specificity_gate: SpecificityGate {
                status: specificity,
                details: "specificity gate".to_string(),
                negative_samples: 1,
                false_positives: u64::from(specificity_failed),
            },
            sensitivity_gate: SensitivityGate {
                status: sensitivity,
                details: "sensitivity gate".to_string(),
                spike_in_detected: u64::from(sensitivity_passed),
                spike_in_total: 1,
            },
        }
    }

    fn passing_negative_control() -> NegativeControlDecisionInput {
        let result = NegativeControlResult {
            control_id: "neg-001".to_string(),
            control_status: "pass".to_string(),
            candidates: Vec::new(),
        };

        NegativeControlDecisionInput::Passed {
            comparator: BackgroundComparator::new(&result),
            result,
        }
    }
}
