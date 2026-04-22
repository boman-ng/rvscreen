use crate::calibration::{GateStatus, ReleaseGate as BenchmarkGates};
use crate::cli::ScreenMode;
use crate::decision::NegativeControlDecisionInput;
use crate::error::{Result, RvScreenError};
use crate::types::{ProfileToml, ReleaseStatus};
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

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct ReleaseGate;

impl ReleaseGate {
    pub fn evaluate(context: &ReleaseContext) -> ReleaseStatus {
        if negative_control_gate_blocked(context)
            || profile_gate_blocked(context)
            || benchmark_gate_blocked(context)
        {
            return ReleaseStatus::Blocked;
        }

        if !has_qualified_negative_control(context)
            || context.sampling_mode != ScreenMode::Representative
        {
            return ReleaseStatus::Provisional;
        }

        ReleaseStatus::Final
    }
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

fn has_qualified_negative_control(context: &ReleaseContext) -> bool {
    matches!(
        context.negative_control,
        NegativeControlDecisionInput::Passed { .. }
    )
}

fn negative_control_gate_blocked(context: &ReleaseContext) -> bool {
    match context.negative_control {
        NegativeControlDecisionInput::Missing => {
            context.calibration_profile.negative_control_required
        }
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

fn benchmark_gate_blocked(context: &ReleaseContext) -> bool {
    let Some(gates) = context.benchmark_gates else {
        return true;
    };

    gate_failed(&gates.backend_gate.status)
        || gate_failed(&gates.reference_gate.status)
        || gate_failed(&gates.specificity_gate.status)
        || gate_failed(&gates.sensitivity_gate.status)
}

fn gate_failed(status: &GateStatus) -> bool {
    *status == GateStatus::Fail
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
    fn missing_negative_control_when_required_is_blocked() {
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

        assert_eq!(status, ReleaseStatus::Blocked);
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
    fn optional_missing_negative_control_is_provisional() {
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

        assert_eq!(status, ReleaseStatus::Provisional);
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
            },
            decision_rules: DecisionRules {
                theta_pos: 0.01,
                theta_neg: 0.0001,
                allow_indeterminate: true,
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
