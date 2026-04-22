use crate::types::{CandidateCall, CandidateRules, DecisionRules, DecisionStatus, StopReason};

const MIN_REPRESENTATIVE_ROUNDS: u64 = 2;
const MAX_AMBIGUOUS_FRACTION_FOR_POSITIVE: f64 = 0.5;

#[derive(Debug, Clone, PartialEq)]
pub struct DecisionEngine {
    rules: DecisionRules,
    candidate_rules: CandidateRules,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DecisionContext {
    pub rounds_completed: u64,
    pub max_rounds: u64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct DecisionOutcome {
    pub status: DecisionStatus,
    pub stop_reason: Option<StopReason>,
}

impl DecisionEngine {
    pub fn new(rules: &DecisionRules, candidate_rules: &CandidateRules) -> Self {
        Self {
            rules: rules.clone(),
            candidate_rules: candidate_rules.clone(),
        }
    }

    pub fn rules(&self) -> &DecisionRules {
        &self.rules
    }

    pub fn candidate_rules(&self) -> &CandidateRules {
        &self.candidate_rules
    }

    pub fn decide(
        &self,
        candidates: &[CandidateCall],
        context: &DecisionContext,
    ) -> DecisionOutcome {
        if candidates
            .iter()
            .any(|candidate| self.is_positive_candidate(candidate))
        {
            return DecisionOutcome::stop(
                DecisionStatus::Positive,
                StopReason::PositiveBoundaryCrossed,
            );
        }

        if self.can_confirm_negative(candidates, context) {
            return DecisionOutcome::stop(
                DecisionStatus::Negative,
                StopReason::NegativeBoundaryConfirmed,
            );
        }

        if context.max_rounds_reached() {
            return DecisionOutcome::stop(
                DecisionStatus::Indeterminate,
                StopReason::MaxRoundsReached,
            );
        }

        debug_assert!(
            self.rules.allow_indeterminate,
            "Task 16 assumes indeterminate remains an allowed outcome"
        );
        DecisionOutcome::continue_with(DecisionStatus::Indeterminate)
    }

    pub fn is_positive_candidate(&self, candidate: &CandidateCall) -> bool {
        candidate.unique_fraction > self.rules.theta_pos
            && candidate.nonoverlap_fragments >= self.candidate_rules.min_nonoverlap_fragments
            && candidate.breadth >= self.candidate_rules.min_breadth
            && candidate.background_ratio >= self.candidate_rules.max_background_ratio
            && self.has_no_strong_host_explanation(candidate)
    }

    pub fn evaluate_candidates(&self, candidates: &[CandidateCall]) -> Vec<CandidateCall> {
        candidates
            .iter()
            .cloned()
            .map(|candidate| self.evaluate_candidate(candidate))
            .collect()
    }

    pub fn can_confirm_negative(
        &self,
        candidates: &[CandidateCall],
        context: &DecisionContext,
    ) -> bool {
        context.has_minimum_representative_rounds()
            && candidates
                .iter()
                .all(|candidate| candidate.unique_fraction < self.rules.theta_neg)
            && candidates
                .iter()
                .all(|candidate| !self.has_abnormal_background(candidate))
    }

    fn has_no_strong_host_explanation(&self, candidate: &CandidateCall) -> bool {
        let Some(accepted_fragments) =
            (candidate.accepted_fragments > 0).then_some(candidate.accepted_fragments as f64)
        else {
            return false;
        };

        (candidate.ambiguous_fragments as f64 / accepted_fragments)
            < MAX_AMBIGUOUS_FRACTION_FOR_POSITIVE
    }

    fn has_abnormal_background(&self, candidate: &CandidateCall) -> bool {
        candidate.background_ratio >= self.candidate_rules.max_background_ratio
    }

    fn evaluate_candidate(&self, mut candidate: CandidateCall) -> CandidateCall {
        let strong_host_explanation = !self.has_no_strong_host_explanation(&candidate);
        let abnormal_background = self.has_abnormal_background(&candidate);
        let mut reasons = std::mem::take(&mut candidate.decision_reasons)
            .into_iter()
            .filter(|reason| reason != "decision_pending_task_16")
            .collect::<Vec<_>>();

        if candidate.unique_fraction > self.rules.theta_pos {
            reasons.push("unique_fraction_above_theta_pos".to_string());
        } else if candidate.unique_fraction < self.rules.theta_neg {
            reasons.push("unique_fraction_below_theta_neg".to_string());
        } else {
            reasons.push("unique_fraction_between_boundaries".to_string());
        }

        if candidate.nonoverlap_fragments >= self.candidate_rules.min_nonoverlap_fragments {
            reasons.push("nonoverlap_fragments_meet_minimum".to_string());
        }
        if candidate.breadth >= self.candidate_rules.min_breadth {
            reasons.push("breadth_meets_minimum".to_string());
        }
        if candidate.background_ratio >= self.candidate_rules.max_background_ratio {
            reasons.push("background_ratio_meets_minimum".to_string());
        } else {
            reasons.push("background_ratio_below_positive_threshold".to_string());
        }
        if strong_host_explanation {
            reasons.push("strong_host_explanation_present".to_string());
        } else {
            reasons.push("no_strong_host_explanation".to_string());
        }
        if abnormal_background {
            reasons.push("background_not_negative_safe".to_string());
        } else {
            reasons.push("background_negative_safe".to_string());
        }

        candidate.decision = if self.is_positive_candidate(&candidate) {
            DecisionStatus::Positive
        } else if candidate.unique_fraction < self.rules.theta_neg && !abnormal_background {
            DecisionStatus::Negative
        } else {
            DecisionStatus::Indeterminate
        };
        candidate.decision_reasons = reasons;
        candidate
    }
}

impl DecisionContext {
    pub fn new(rounds_completed: u64, max_rounds: u64) -> Self {
        Self {
            rounds_completed,
            max_rounds,
        }
    }

    pub fn has_minimum_representative_rounds(&self) -> bool {
        self.rounds_completed >= MIN_REPRESENTATIVE_ROUNDS
    }

    pub fn max_rounds_reached(&self) -> bool {
        self.max_rounds > 0 && self.rounds_completed >= self.max_rounds
    }
}

impl DecisionOutcome {
    fn stop(status: DecisionStatus, stop_reason: StopReason) -> Self {
        Self {
            status,
            stop_reason: Some(stop_reason),
        }
    }

    fn continue_with(status: DecisionStatus) -> Self {
        Self {
            status,
            stop_reason: None,
        }
    }

    pub fn should_stop(&self) -> bool {
        self.stop_reason.is_some()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::{DecisionStatus, EvidenceStrength};

    #[test]
    fn positive_candidate_crosses_boundary_and_stops() {
        let engine = decision_engine();
        let outcome = engine.decide(
            &[candidate_call()
                .with_unique_fraction(0.001)
                .with_nonoverlap_fragments(10)
                .with_breadth(0.01)
                .with_background_ratio(5.0)
                .with_counts(10, 1)
                .build()],
            &DecisionContext::new(1, 5),
        );

        assert_eq!(outcome.status, DecisionStatus::Positive);
        assert_eq!(
            outcome.stop_reason,
            Some(StopReason::PositiveBoundaryCrossed)
        );
        assert!(outcome.should_stop());
    }

    #[test]
    fn positive_requires_all_conditions() {
        let engine = decision_engine();
        let outcome = engine.decide(
            &[candidate_call()
                .with_unique_fraction(0.001)
                .with_nonoverlap_fragments(10)
                .with_breadth(0.01)
                .with_background_ratio(5.0)
                .with_counts(10, 5)
                .build()],
            &DecisionContext::new(1, 5),
        );

        assert_eq!(outcome.status, DecisionStatus::Indeterminate);
        assert_eq!(outcome.stop_reason, None);
        assert!(!outcome.should_stop());
    }

    #[test]
    fn theta_neg_minus_epsilon_confirms_negative_after_two_rounds() {
        let engine = decision_engine();
        let epsilon = 1e-9;
        let outcome = engine.decide(
            &[
                candidate_call()
                    .with_unique_fraction(engine.rules().theta_neg - epsilon)
                    .with_background_ratio(0.0)
                    .build(),
                candidate_call()
                    .with_unique_fraction(0.0)
                    .with_background_ratio(0.1)
                    .build(),
            ],
            &DecisionContext::new(3, 5),
        );

        assert_eq!(outcome.status, DecisionStatus::Negative);
        assert_eq!(
            outcome.stop_reason,
            Some(StopReason::NegativeBoundaryConfirmed)
        );
        assert!(outcome.should_stop());
    }

    #[test]
    fn negative_requires_minimum_representative_rounds() {
        let engine = decision_engine();
        let outcome = engine.decide(
            &[candidate_call().with_unique_fraction(0.0).build()],
            &DecisionContext::new(1, 5),
        );

        assert_eq!(outcome.status, DecisionStatus::Indeterminate);
        assert_eq!(outcome.stop_reason, None);
        assert!(!outcome.should_stop());
    }

    #[test]
    fn abnormal_background_blocks_negative_confirmation() {
        let engine = decision_engine();
        let outcome = engine.decide(
            &[candidate_call()
                .with_unique_fraction(0.0)
                .with_background_ratio(engine.candidate_rules().max_background_ratio)
                .build()],
            &DecisionContext::new(3, 5),
        );

        assert_eq!(outcome.status, DecisionStatus::Indeterminate);
        assert_eq!(outcome.stop_reason, None);
    }

    #[test]
    fn gray_zone_candidate_is_indeterminate_before_max_rounds() {
        let engine = decision_engine();
        let outcome = engine.decide(
            &[candidate_call().with_unique_fraction(0.00003).build()],
            &DecisionContext::new(1, 5),
        );

        assert_eq!(outcome.status, DecisionStatus::Indeterminate);
        assert_eq!(outcome.stop_reason, None);
        assert!(!outcome.should_stop());
    }

    #[test]
    fn max_rounds_reached_stops_with_indeterminate() {
        let engine = decision_engine();
        let outcome = engine.decide(
            &[candidate_call().with_unique_fraction(0.00003).build()],
            &DecisionContext::new(5, 5),
        );

        assert_eq!(outcome.status, DecisionStatus::Indeterminate);
        assert_eq!(outcome.stop_reason, Some(StopReason::MaxRoundsReached));
        assert!(outcome.should_stop());
    }

    #[test]
    fn theta_pos_is_strictly_greater_than_boundary() {
        let engine = decision_engine();
        let outcome = engine.decide(
            &[candidate_call()
                .with_unique_fraction(engine.rules().theta_pos)
                .with_nonoverlap_fragments(10)
                .with_breadth(0.01)
                .with_background_ratio(5.0)
                .with_counts(10, 0)
                .build()],
            &DecisionContext::new(1, 5),
        );

        assert_eq!(outcome.status, DecisionStatus::Indeterminate);
        assert_eq!(outcome.stop_reason, None);
    }

    #[test]
    fn evaluate_candidates_sets_positive_negative_and_indeterminate_labels() {
        let engine = decision_engine();
        let candidates = engine.evaluate_candidates(&[
            candidate_call()
                .with_unique_fraction(0.001)
                .with_nonoverlap_fragments(10)
                .with_breadth(0.01)
                .with_background_ratio(5.0)
                .with_counts(10, 1)
                .build(),
            candidate_call()
                .with_unique_fraction(0.0)
                .with_background_ratio(0.0)
                .build(),
            candidate_call()
                .with_unique_fraction(0.00003)
                .with_background_ratio(0.0)
                .build(),
        ]);

        assert_eq!(candidates[0].decision, DecisionStatus::Positive);
        assert!(candidates[0]
            .decision_reasons
            .iter()
            .any(|reason| reason == "no_strong_host_explanation"));

        assert_eq!(candidates[1].decision, DecisionStatus::Negative);
        assert!(candidates[1]
            .decision_reasons
            .iter()
            .any(|reason| reason == "background_negative_safe"));

        assert_eq!(candidates[2].decision, DecisionStatus::Indeterminate);
        assert!(candidates[2]
            .decision_reasons
            .iter()
            .any(|reason| reason == "unique_fraction_between_boundaries"));
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
                decision_reasons: vec!["decision_pending_task_16".to_string()],
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

        fn with_background_ratio(mut self, background_ratio: f64) -> Self {
            self.inner.background_ratio = background_ratio;
            self
        }

        fn with_counts(mut self, accepted_fragments: u64, ambiguous_fragments: u64) -> Self {
            self.inner.accepted_fragments = accepted_fragments;
            self.inner.ambiguous_fragments = ambiguous_fragments;
            self
        }
    }
}
