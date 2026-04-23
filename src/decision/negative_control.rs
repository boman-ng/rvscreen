use crate::error::{Result, RvScreenError};
use crate::types::{CandidateCall, NegativeControlManifest, NegativeControlStatus};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::fs;
use std::path::Path;

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct NegativeControlResult {
    pub control_id: String,
    pub control_status: String,
    #[serde(default)]
    pub candidates: Vec<NegativeControlCandidate>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct NegativeControlCandidate {
    pub accession_or_group: String,
    pub unique_fraction: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct BackgroundComparator {
    background_by_candidate: BTreeMap<String, f64>,
}

#[derive(Debug, Clone, PartialEq)]
pub enum NegativeControlDecisionInput {
    Missing,
    Failed(NegativeControlResult),
    Passed {
        result: NegativeControlResult,
        comparator: BackgroundComparator,
    },
}

impl BackgroundComparator {
    pub fn new(control: &NegativeControlResult) -> Self {
        let mut background_by_candidate = BTreeMap::new();

        for candidate in &control.candidates {
            background_by_candidate
                .entry(candidate.accession_or_group.clone())
                .and_modify(|existing| {
                    if candidate.unique_fraction > *existing {
                        *existing = candidate.unique_fraction;
                    }
                })
                .or_insert(candidate.unique_fraction);
        }

        Self {
            background_by_candidate,
        }
    }

    pub fn background_fraction(&self, accession_or_group: &str) -> Option<f64> {
        self.background_by_candidate
            .get(accession_or_group)
            .copied()
    }

    pub fn compute_background_ratio(&self, candidate: &CandidateCall) -> f64 {
        match self.background_fraction(&candidate.accession_or_group) {
            Some(control_background_fraction) if control_background_fraction > 0.0 => {
                candidate.unique_fraction / control_background_fraction
            }
            Some(_) | None if candidate.unique_fraction > 0.0 => f64::INFINITY,
            Some(_) | None => 0.0,
        }
    }
}

impl NegativeControlDecisionInput {
    pub fn from_optional_path(path: Option<&Path>) -> Result<Self> {
        let Some(path) = path else {
            return Ok(Self::Missing);
        };

        validate_negative_control_file(path)?;
        let result = load_negative_control(path)?;

        if result.control_status == "pass" {
            Ok(Self::Passed {
                comparator: BackgroundComparator::new(&result),
                result,
            })
        } else {
            Ok(Self::Failed(result))
        }
    }

    pub fn manifest(&self, required: bool) -> NegativeControlManifest {
        match self {
            Self::Missing => NegativeControlManifest {
                required,
                status: NegativeControlStatus::Missing,
                control_id: None,
                control_status: None,
            },
            Self::Failed(result) => NegativeControlManifest {
                required,
                status: NegativeControlStatus::Fail,
                control_id: Some(result.control_id.clone()),
                control_status: Some(result.control_status.clone()),
            },
            Self::Passed { result, .. } => NegativeControlManifest {
                required,
                status: NegativeControlStatus::Pass,
                control_id: Some(result.control_id.clone()),
                control_status: Some(result.control_status.clone()),
            },
        }
    }
}

pub fn load_negative_control(path: &Path) -> Result<NegativeControlResult> {
    let body = fs::read_to_string(path).map_err(|source| RvScreenError::io(path, source))?;
    let result: NegativeControlResult = serde_json::from_str(&body).map_err(|error| {
        RvScreenError::parse(
            path,
            error.line() as u64,
            format!("column {}: {error}", error.column()),
        )
    })?;
    validate_negative_control_result(&result)?;
    Ok(result)
}

pub fn apply_negative_control(
    candidates: &[CandidateCall],
    negative_control: &NegativeControlDecisionInput,
) -> Vec<CandidateCall> {
    candidates
        .iter()
        .cloned()
        .map(|candidate| annotate_candidate(candidate, negative_control))
        .collect()
}

fn validate_negative_control_file(path: &Path) -> Result<()> {
    let metadata = fs::metadata(path).map_err(|source| RvScreenError::io(path, source))?;
    if !metadata.is_file() {
        return Err(RvScreenError::validation(
            "negative_control",
            format!("`{}` is not a regular file", path.display()),
        ));
    }

    if metadata.len() == 0 {
        return Err(RvScreenError::validation(
            "negative_control",
            format!("`{}` is empty", path.display()),
        ));
    }

    Ok(())
}

fn validate_negative_control_result(result: &NegativeControlResult) -> Result<()> {
    if result.control_id.trim().is_empty() {
        return Err(RvScreenError::validation(
            "negative_control.control_id",
            "control_id must not be empty",
        ));
    }

    if result.control_status.trim().is_empty() {
        return Err(RvScreenError::validation(
            "negative_control.control_status",
            "control_status must not be empty",
        ));
    }

    for candidate in &result.candidates {
        if candidate.accession_or_group.trim().is_empty() {
            return Err(RvScreenError::validation(
                "negative_control.candidates.accession_or_group",
                "accession_or_group must not be empty",
            ));
        }

        if !candidate.unique_fraction.is_finite() || candidate.unique_fraction < 0.0 {
            return Err(RvScreenError::validation(
                "negative_control.candidates.unique_fraction",
                format!(
                    "unique_fraction for `{}` must be a finite non-negative value",
                    candidate.accession_or_group
                ),
            ));
        }
    }

    Ok(())
}

fn annotate_candidate(
    mut candidate: CandidateCall,
    negative_control: &NegativeControlDecisionInput,
) -> CandidateCall {
    let mut reasons = std::mem::take(&mut candidate.decision_reasons)
        .into_iter()
        .filter(|reason| {
            reason != "background_ratio_pending_task_19"
                && reason != "background_ratio_computed_from_negative_control"
                && !reason.starts_with("background_ratio_unavailable_")
                && !reason.starts_with("negative_control_")
        })
        .collect::<Vec<_>>();

    match negative_control {
        NegativeControlDecisionInput::Missing => {
            candidate.background_ratio = 0.0;
            reasons.push("background_ratio_unavailable_missing_negative_control".to_string());
        }
        NegativeControlDecisionInput::Failed(result) => {
            candidate.background_ratio = 0.0;
            reasons.push(format!("negative_control_id={}", result.control_id));
            reasons.push(format!("negative_control_status={}", result.control_status));
            reasons.push("background_ratio_unavailable_negative_control_not_pass".to_string());
        }
        NegativeControlDecisionInput::Passed { result, comparator } => {
            candidate.background_ratio = comparator.compute_background_ratio(&candidate);
            reasons.push(format!("negative_control_id={}", result.control_id));
            reasons.push("negative_control_status=pass".to_string());

            match comparator.background_fraction(&candidate.accession_or_group) {
                Some(control_background_fraction) if control_background_fraction > 0.0 => {
                    reasons.push(format!(
                        "negative_control_background_fraction={control_background_fraction:.12e}"
                    ));
                    reasons.push("background_ratio_computed_from_negative_control".to_string());
                }
                Some(_) => {
                    reasons.push("negative_control_background_fraction_zero".to_string());
                    reasons.push("background_ratio_computed_from_negative_control".to_string());
                }
                None => {
                    reasons.push("negative_control_background_fraction_absent".to_string());
                    reasons.push("background_ratio_computed_from_negative_control".to_string());
                }
            }
        }
    }

    candidate.decision_reasons = reasons;
    candidate
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::decision::DecisionEngine;
    use crate::types::{CandidateRules, DecisionRules, DecisionStatus, EvidenceStrength};

    #[test]
    fn parses_negative_control_json() {
        let parsed: NegativeControlResult = serde_json::from_str(
            r#"{
                "control_id": "batch_neg_001",
                "control_status": "pass",
                "candidates": [
                    {
                        "accession_or_group": "NC_007605.1",
                        "unique_fraction": 0.000002
                    }
                ]
            }"#,
        )
        .expect("negative-control JSON should parse");

        assert_eq!(parsed.control_id, "batch_neg_001");
        assert_eq!(parsed.control_status, "pass");
        assert_eq!(parsed.candidates.len(), 1);
        assert_eq!(parsed.candidates[0].accession_or_group, "NC_007605.1");
        assert_eq!(parsed.candidates[0].unique_fraction, 0.000002);
    }

    #[test]
    fn ratio_below_threshold_blocks_positive() {
        let engine = decision_engine();
        let control = NegativeControlResult {
            control_id: "neg-001".to_string(),
            control_status: "pass".to_string(),
            candidates: vec![NegativeControlCandidate {
                accession_or_group: "ACC-A".to_string(),
                unique_fraction: 0.00008,
            }],
        };
        let compared = apply_negative_control(
            &[candidate_call()
                .with_unique_fraction(0.0001)
                .with_nonoverlap_fragments(10)
                .with_breadth(0.01)
                .with_counts(10, 1)
                .build()],
            &NegativeControlDecisionInput::Passed {
                comparator: BackgroundComparator::new(&control),
                result: control,
            },
        );

        assert_eq!(compared[0].background_ratio, 1.25);

        let evaluated = engine.evaluate_candidates(&compared);
        assert_eq!(evaluated[0].decision, DecisionStatus::Indeterminate);
        assert!(evaluated[0]
            .decision_reasons
            .iter()
            .any(|reason| reason == "background_ratio_below_positive_threshold"));
    }

    #[test]
    fn high_ratio_allows_positive_path() {
        let engine = decision_engine();
        let control = NegativeControlResult {
            control_id: "neg-002".to_string(),
            control_status: "pass".to_string(),
            candidates: vec![NegativeControlCandidate {
                accession_or_group: "ACC-A".to_string(),
                unique_fraction: 0.000001,
            }],
        };
        let compared = apply_negative_control(
            &[candidate_call()
                .with_unique_fraction(0.001)
                .with_nonoverlap_fragments(10)
                .with_breadth(0.01)
                .with_counts(10, 1)
                .build()],
            &NegativeControlDecisionInput::Passed {
                comparator: BackgroundComparator::new(&control),
                result: control,
            },
        );

        assert!((compared[0].background_ratio - 1000.0).abs() < 1e-9);

        let evaluated = engine.evaluate_candidates(&compared);
        assert_eq!(evaluated[0].decision, DecisionStatus::Positive);
        assert!(evaluated[0]
            .decision_reasons
            .iter()
            .any(|reason| reason == "background_ratio_meets_minimum"));
    }

    #[test]
    fn absent_or_zero_control_background_is_explicit() {
        let control = NegativeControlResult {
            control_id: "neg-003".to_string(),
            control_status: "pass".to_string(),
            candidates: vec![NegativeControlCandidate {
                accession_or_group: "OTHER".to_string(),
                unique_fraction: 0.0,
            }],
        };
        let compared = apply_negative_control(
            &[candidate_call().with_unique_fraction(0.001).build()],
            &NegativeControlDecisionInput::Passed {
                comparator: BackgroundComparator::new(&control),
                result: control,
            },
        );

        assert!(compared[0].background_ratio.is_infinite());
        assert!(compared[0]
            .decision_reasons
            .iter()
            .any(|reason| reason == "negative_control_background_fraction_absent"));
    }

    fn decision_engine() -> DecisionEngine {
        DecisionEngine::new(
            &DecisionRules {
                theta_pos: 0.00005,
                theta_neg: 0.000005,
                allow_indeterminate: true,
            },
            &CandidateRules {
                min_nonoverlap_fragments: 3,
                min_breadth: 0.001,
                max_background_ratio: 1.5,
            },
        )
    }

    fn candidate_call() -> CandidateCallBuilder {
        CandidateCallBuilder {
            inner: CandidateCall {
                virus_name: "virus-a".to_string(),
                taxid: 10_376,
                accession_or_group: "ACC-A".to_string(),
                accepted_fragments: 4,
                nonoverlap_fragments: 4,
                raw_fraction: 0.0,
                unique_fraction: 0.0,
                fraction_ci_95: [0.0, 0.0],
                clopper_pearson_upper: 0.0,
                breadth: 0.0,
                ambiguous_fragments: 0,
                background_ratio: 0.0,
                decision: DecisionStatus::Indeterminate,
                decision_reasons: vec!["background_ratio_pending_task_19".to_string()],
                evidence_strength: EvidenceStrength::Low,
            },
        }
    }

    struct CandidateCallBuilder {
        inner: CandidateCall,
    }

    impl CandidateCallBuilder {
        fn build(self) -> CandidateCall {
            self.inner
        }

        fn with_unique_fraction(mut self, unique_fraction: f64) -> Self {
            self.inner.unique_fraction = unique_fraction;
            self
        }

        fn with_nonoverlap_fragments(mut self, nonoverlap_fragments: u64) -> Self {
            self.inner.nonoverlap_fragments = nonoverlap_fragments;
            self
        }

        fn with_breadth(mut self, breadth: f64) -> Self {
            self.inner.breadth = breadth;
            self
        }

        fn with_counts(mut self, accepted_fragments: u64, ambiguous_fragments: u64) -> Self {
            self.inner.accepted_fragments = accepted_fragments;
            self.inner.ambiguous_fragments = ambiguous_fragments;
            self
        }
    }
}
