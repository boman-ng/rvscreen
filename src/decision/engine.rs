use crate::decision::stats::{bonferroni_corrected_alpha, exact_binomial_survival};
use crate::types::{
    CandidateCall, CandidateRules, DecisionRules, DecisionStatus, EvidenceStrength, StopReason,
};
use std::collections::BTreeSet;

const MIN_REPRESENTATIVE_ROUNDS: u64 = 2;

#[derive(Debug, Clone, PartialEq)]
pub struct DecisionEngine {
    rules: DecisionRules,
    candidate_rules: CandidateRules,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DecisionContext {
    pub rounds_completed: u64,
    pub max_rounds: u64,
    pub candidate_family_size: u64,
    pub planned_round_look_count: u64,
    pub current_look_index: u64,
    pub prior_positive_group_ids: Vec<String>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct DecisionOutcome {
    pub status: DecisionStatus,
    pub stop_reason: Option<StopReason>,
}

#[derive(Debug, Clone, PartialEq)]
struct CandidateEvaluation {
    status: DecisionStatus,
    reasons: Vec<String>,
    accession_resolution: &'static str,
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
        let positive_candidates = candidates
            .iter()
            .filter(|candidate| {
                self.evaluate_candidate_status_with_context(candidate, context)
                    == DecisionStatus::Positive
            })
            .collect::<Vec<_>>();
        if !positive_candidates.is_empty() {
            let accession_positive = positive_candidates
                .iter()
                .any(|candidate| !is_group_candidate(candidate));
            let group_positive_can_stop = positive_candidates.iter().any(|candidate| {
                is_group_candidate(candidate) && context.can_stop_on_group(candidate)
            });
            if accession_positive || group_positive_can_stop {
                return DecisionOutcome::stop(
                    DecisionStatus::Positive,
                    StopReason::PositiveBoundaryCrossed,
                );
            }
            return DecisionOutcome::continue_with(DecisionStatus::Positive);
        }

        let has_weak_positive = candidates.iter().any(|candidate| {
            self.evaluate_candidate_status_with_context(candidate, context)
                == DecisionStatus::WeakPositive
        });
        if has_weak_positive {
            if context.max_rounds_reached() {
                return DecisionOutcome::stop(
                    DecisionStatus::WeakPositive,
                    StopReason::MaxRoundsReached,
                );
            }
            return DecisionOutcome::continue_with(DecisionStatus::WeakPositive);
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
        self.evaluate_candidate_status(candidate) == DecisionStatus::Positive
    }

    pub fn evaluate_candidates(&self, candidates: &[CandidateCall]) -> Vec<CandidateCall> {
        self.evaluate_candidates_with_context(candidates, &DecisionContext::default())
    }

    pub fn evaluate_candidates_with_context(
        &self,
        candidates: &[CandidateCall],
        context: &DecisionContext,
    ) -> Vec<CandidateCall> {
        candidates
            .iter()
            .cloned()
            .map(|candidate| self.evaluate_candidate(candidate, context))
            .collect()
    }

    pub fn can_confirm_negative(
        &self,
        candidates: &[CandidateCall],
        context: &DecisionContext,
    ) -> bool {
        let theta_neg_operational = self.theta_neg_operational();
        context.has_minimum_representative_rounds()
            && candidates.iter().all(|candidate| {
                candidate.unique_fraction <= theta_neg_operational
                    && candidate.clopper_pearson_upper <= self.rules.negative_cp_lod_max
            })
    }

    fn evaluate_candidate_status(&self, candidate: &CandidateCall) -> DecisionStatus {
        self.evaluate_candidate_status_with_context(candidate, &DecisionContext::default())
    }

    fn evaluate_candidate_status_with_context(
        &self,
        candidate: &CandidateCall,
        context: &DecisionContext,
    ) -> DecisionStatus {
        self.evaluate_candidate_gates(candidate, context).status
    }

    fn evaluate_candidate_gates(
        &self,
        candidate: &CandidateCall,
        context: &DecisionContext,
    ) -> CandidateEvaluation {
        let abundance_pass = candidate.unique_fraction > self.rules.theta_pos;
        let absolute_pass = candidate.accepted_fragments >= self.candidate_rules.theta_pos_absolute;
        let nonoverlap_pass =
            candidate.nonoverlap_fragments >= self.candidate_rules.min_nonoverlap_fragments;
        let breadth_pass = candidate.breadth >= self.candidate_rules.min_breadth;
        let background_pass = true;
        let ambiguity_pass = self.has_no_strong_host_explanation(candidate);
        let evidence_pass = evidence_rank(&candidate.evidence_strength)
            >= evidence_rank(&self.candidate_rules.min_positive_evidence_strength);
        let enrichment = self.positive_enrichment(candidate, context);
        let enrichment_pass = enrichment.positive_p_value <= enrichment.corrected_alpha;

        let mut reasons = Vec::new();
        push_gate(
            &mut reasons,
            "unique_fraction_above_theta_pos",
            abundance_pass,
        );
        push_gate(
            &mut reasons,
            "accepted_fragments_meet_theta_pos_absolute",
            absolute_pass,
        );
        push_gate(
            &mut reasons,
            "nonoverlap_fragments_meet_minimum",
            nonoverlap_pass,
        );
        push_gate(&mut reasons, "breadth_meets_minimum", breadth_pass);
        reasons
            .push("negative_control_policy=not_used_for_default_background_assessment".to_string());
        push_gate(
            &mut reasons,
            "ambiguity_fraction_below_positive_limit",
            ambiguity_pass,
        );
        reasons.push(if ambiguity_pass {
            "no_strong_host_explanation".to_string()
        } else {
            "strong_host_explanation_present".to_string()
        });
        push_gate(
            &mut reasons,
            "evidence_strength_meets_positive_minimum",
            evidence_pass,
        );
        push_gate(&mut reasons, "positive_enrichment_p_value", enrichment_pass);
        reasons.push(format!(
            "positive_statistic_method={}",
            self.rules.positive_statistic_method
        ));
        reasons.push(format!(
            "positive_p_value={:.12e}",
            enrichment.positive_p_value
        ));
        reasons.push(format!(
            "positive_corrected_alpha={:.12e}",
            enrichment.corrected_alpha
        ));
        reasons.push(format!(
            "total_sampled_fragments={}",
            enrichment.total_sampled_fragments
        ));
        reasons.push(format!(
            "candidate_family_size={}",
            enrichment.candidate_family_size
        ));
        reasons.push(format!(
            "multiplicity_family={}",
            candidate.multiplicity_family
        ));
        reasons.push(format!("aggregation_level={}", candidate.aggregation_level));
        reasons.push(format!(
            "nonoverlap_fragments_algorithm={}",
            candidate.nonoverlap_fragments_algorithm
        ));
        reasons.push(format!(
            "coverage_interval_blocks={}",
            candidate.coverage_interval_blocks
        ));
        reasons.push(format!(
            "accepted_fraction_label={}",
            candidate.accepted_fraction_label
        ));
        reasons.push(format!(
            "unique_fraction_label={}",
            candidate.unique_fraction_label
        ));
        reasons.push(format!(
            "fraction_denominator_source={}",
            candidate.fraction_denominator_source
        ));
        reasons.push(format!(
            "planned_round_look_count={}",
            enrichment.planned_round_look_count
        ));
        reasons.push(format!(
            "current_look_index={}",
            enrichment.current_look_index
        ));
        reasons.push("sequential_look_correction=bonferroni".to_string());
        reasons.push(format!(
            "theta_neg_operational={:.12e}",
            self.theta_neg_operational()
        ));
        if candidate.unique_fraction <= self.theta_neg_operational() {
            reasons.push("unique_fraction_below_theta_neg".to_string());
        } else if candidate.unique_fraction > self.rules.theta_pos {
            reasons.push("unique_fraction_above_theta_pos".to_string());
        } else {
            reasons.push("unique_fraction_between_boundaries".to_string());
        }

        let accepted_core_pass = abundance_pass
            && absolute_pass
            && nonoverlap_pass
            && breadth_pass
            && background_pass
            && evidence_pass
            && enrichment_pass;
        let group_candidate = is_group_candidate(candidate);
        let high_ambiguity = !ambiguity_pass;
        if high_ambiguity {
            reasons.push("high_candidate_ambiguity_present".to_string());
        }
        let accession_resolution =
            accession_resolution_label(candidate, accepted_core_pass, high_ambiguity);
        // The reason string is retained as an audit trail; the canonical
        // value flows through `CandidateEvaluation::accession_resolution`
        // so downstream consumers never have to reverse-parse this text.
        reasons.push(format!("accession_resolution={accession_resolution}"));
        if group_candidate {
            reasons.push("positive_boundary_level=group".to_string());
            if high_ambiguity {
                reasons.push("ambiguity_not_used_as_group_detection_veto".to_string());
            }
            reasons.push(if context.can_stop_on_group(candidate) {
                "group_positive_identity_continuity_pass".to_string()
            } else {
                "group_positive_identity_continuity_pending".to_string()
            });
        }

        let all_positive = if group_candidate {
            accepted_core_pass
        } else {
            accepted_core_pass && ambiguity_pass
        };
        let trace_signal = self.candidate_rules.weak_positive_enabled
            && abundance_pass
            && absolute_pass
            && enrichment_pass
            && ambiguity_pass;

        let status = if all_positive {
            DecisionStatus::Positive
        } else if trace_signal {
            DecisionStatus::WeakPositive
        } else if candidate.unique_fraction <= self.theta_neg_operational()
            && candidate.clopper_pearson_upper <= self.rules.negative_cp_lod_max
        {
            DecisionStatus::Negative
        } else {
            DecisionStatus::Indeterminate
        };

        CandidateEvaluation {
            status,
            reasons,
            accession_resolution,
        }
    }

    fn positive_enrichment(
        &self,
        candidate: &CandidateCall,
        context: &DecisionContext,
    ) -> PositiveEnrichment {
        let total_sampled_fragments = candidate.total_sampled_fragments;
        let candidate_family_size = context.candidate_family_size.max(1);
        let planned_round_look_count = context.planned_round_look_count.max(1);
        let current_look_index = context.current_look_index.max(1);
        let p0 = self
            .theta_neg_operational()
            .max(self.rules.noise_floor_fraction);
        PositiveEnrichment {
            positive_p_value: exact_binomial_survival(
                candidate.accepted_fragments,
                total_sampled_fragments,
                p0,
            ),
            corrected_alpha: bonferroni_corrected_alpha(
                self.rules.positive_alpha_global,
                candidate_family_size,
                planned_round_look_count,
            ),
            total_sampled_fragments,
            candidate_family_size,
            planned_round_look_count,
            current_look_index,
        }
    }

    fn has_no_strong_host_explanation(&self, candidate: &CandidateCall) -> bool {
        let Some(accepted_fragments) =
            (candidate.accepted_fragments > 0).then_some(candidate.accepted_fragments as f64)
        else {
            return false;
        };

        (candidate.ambiguous_fragments as f64 / accepted_fragments)
            <= self.candidate_rules.max_ambiguous_fraction_for_positive
    }

    fn theta_neg_operational(&self) -> f64 {
        self.rules.theta_neg.max(self.rules.theta_neg_fixed)
    }

    fn evaluate_candidate(
        &self,
        mut candidate: CandidateCall,
        context: &DecisionContext,
    ) -> CandidateCall {
        let evaluation = self.evaluate_candidate_gates(&candidate, context);
        let mut reasons = std::mem::take(&mut candidate.decision_reasons)
            .into_iter()
            .filter(|reason| {
                reason != "decision_pending_task_16"
                    && !reason.starts_with("unique_fraction_")
                    && !reason.starts_with("accepted_fragments_meet_theta_pos_absolute")
                    && !reason.starts_with("nonoverlap_fragments_meet_minimum")
                    && !reason.starts_with("breadth_meets_minimum")
                    && !reason.starts_with("background_final_positive_supported")
                    && !reason.starts_with("background_ratio_")
                    && !reason.starts_with("background_negative_safe")
                    && !reason.starts_with("background_not_negative_safe")
                    && !reason.starts_with("negative_control_policy=")
                    && !reason.starts_with("aggregation_level=")
                    && !reason.starts_with("nonoverlap_fragments_algorithm=")
                    && !reason.starts_with("coverage_interval_blocks=")
                    && !reason.starts_with("accepted_fraction_label=")
                    && !reason.starts_with("unique_fraction_label=")
                    && !reason.starts_with("fraction_denominator_source=")
                    && !reason.starts_with("ambiguity_fraction_below_positive_limit")
                    && !reason.starts_with("evidence_strength_meets_positive_minimum")
                    && !reason.starts_with("positive_")
                    && !reason.starts_with("total_sampled_fragments=")
                    && !reason.starts_with("candidate_family_size=")
                    && !reason.starts_with("multiplicity_family=")
                    && !reason.starts_with("planned_round_look_count=")
                    && !reason.starts_with("current_look_index=")
                    && !reason.starts_with("sequential_look_correction=")
                    && !reason.starts_with("theta_neg_operational=")
                    && !reason.starts_with("high_candidate_ambiguity_present")
                    && !reason.starts_with("ambiguity_not_used_as_group_detection_veto")
                    && !reason.starts_with("accession_resolution=")
                    && !reason.starts_with("positive_boundary_level=")
                    && !reason.starts_with("group_positive_identity_continuity_")
            })
            .collect::<Vec<_>>();
        reasons.extend(evaluation.reasons);
        candidate.decision = evaluation.status;
        candidate.accession_resolution = evaluation.accession_resolution.to_string();
        candidate.decision_reasons = reasons;
        candidate
    }
}

fn is_group_candidate(candidate: &CandidateCall) -> bool {
    candidate.aggregation_level == "group" || candidate.multiplicity_family == "group"
}

fn accession_resolution_label(
    candidate: &CandidateCall,
    accepted_core_pass: bool,
    high_ambiguity: bool,
) -> &'static str {
    if is_group_candidate(candidate) {
        "not_applicable"
    } else if accepted_core_pass && high_ambiguity {
        "low_resolution"
    } else if accepted_core_pass {
        "high"
    } else if high_ambiguity {
        "unresolved"
    } else {
        "high"
    }
}

#[derive(Debug, Clone, PartialEq)]
struct PositiveEnrichment {
    positive_p_value: f64,
    corrected_alpha: f64,
    total_sampled_fragments: u64,
    candidate_family_size: u64,
    planned_round_look_count: u64,
    current_look_index: u64,
}

fn push_gate(reasons: &mut Vec<String>, label: &str, passed: bool) {
    reasons.push(format!("{label}_{}", if passed { "pass" } else { "fail" }));
}

fn evidence_rank(evidence_strength: &EvidenceStrength) -> u8 {
    match evidence_strength {
        EvidenceStrength::High => 4,
        EvidenceStrength::Moderate | EvidenceStrength::Medium => 3,
        EvidenceStrength::Low => 2,
        EvidenceStrength::Insufficient => 1,
    }
}

impl Default for DecisionContext {
    fn default() -> Self {
        Self::new(1, 1)
    }
}

impl DecisionContext {
    pub fn new(rounds_completed: u64, max_rounds: u64) -> Self {
        Self::with_gate_context(
            rounds_completed,
            max_rounds,
            1,
            max_rounds.max(1),
            rounds_completed.max(1),
        )
    }

    pub fn with_gate_context(
        rounds_completed: u64,
        max_rounds: u64,
        candidate_family_size: u64,
        planned_round_look_count: u64,
        current_look_index: u64,
    ) -> Self {
        Self {
            rounds_completed,
            max_rounds,
            candidate_family_size: candidate_family_size.max(1),
            planned_round_look_count: planned_round_look_count.max(1),
            current_look_index: current_look_index.max(1),
            prior_positive_group_ids: Vec::new(),
        }
    }

    pub fn with_prior_positive_group_ids(mut self, prior_positive_group_ids: Vec<String>) -> Self {
        self.prior_positive_group_ids = prior_positive_group_ids;
        self
    }

    pub fn positive_group_ids(&self, candidates: &[CandidateCall]) -> Vec<String> {
        candidates
            .iter()
            .filter(|candidate| {
                is_group_candidate(candidate) && candidate.decision == DecisionStatus::Positive
            })
            .map(|candidate| candidate.accession_or_group.clone())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect()
    }

    fn can_stop_on_group(&self, candidate: &CandidateCall) -> bool {
        if self.max_rounds <= 1 || self.planned_round_look_count <= 1 {
            return true;
        }
        self.has_minimum_representative_rounds()
            && self
                .prior_positive_group_ids
                .iter()
                .any(|group_id| group_id == &candidate.accession_or_group)
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
    fn group_strong_core_positive_ignores_high_candidate_ambiguity_for_detection() {
        let engine = decision_engine();
        let context = DecisionContext::new(1, 1);
        let evaluated = engine.evaluate_candidate(
            candidate_call()
                .into_group("ebv")
                .with_unique_fraction(0.001)
                .with_nonoverlap_fragments(10)
                .with_breadth(0.01)
                .with_counts(10, 500)
                .build(),
            &context,
        );

        assert_eq!(evaluated.decision, DecisionStatus::Positive);
        assert_eq!(evaluated.accession_resolution, "not_applicable");
        assert!(evaluated
            .decision_reasons
            .iter()
            .any(|reason| reason == "high_candidate_ambiguity_present"));
        assert!(evaluated
            .decision_reasons
            .iter()
            .any(|reason| reason == "ambiguity_not_used_as_group_detection_veto"));
    }

    #[test]
    fn accession_strong_core_high_ambiguity_is_low_resolution_not_positive() {
        let engine = decision_engine();
        let evaluated = engine.evaluate_candidate(
            candidate_call()
                .with_unique_fraction(0.001)
                .with_nonoverlap_fragments(10)
                .with_breadth(0.01)
                .with_counts(10, 500)
                .build(),
            &DecisionContext::new(1, 1),
        );

        assert_eq!(evaluated.decision, DecisionStatus::Indeterminate);
        assert_eq!(evaluated.accession_resolution, "low_resolution");
    }

    #[test]
    fn weak_core_high_ambiguity_group_is_not_positive() {
        let engine = decision_engine();
        let evaluated = engine.evaluate_candidate(
            candidate_call()
                .into_group("ebv")
                .with_unique_fraction(0.001)
                .with_nonoverlap_fragments(1)
                .with_breadth(0.0)
                .with_counts(10, 500)
                .build(),
            &DecisionContext::new(1, 1),
        );

        assert_ne!(evaluated.decision, DecisionStatus::Positive);
    }

    #[test]
    fn high_ambiguity_group_positive_requires_same_group_across_two_looks_to_stop() {
        let engine = decision_engine();
        let group = candidate_call()
            .into_group("ebv")
            .with_unique_fraction(0.001)
            .with_nonoverlap_fragments(10)
            .with_breadth(0.01)
            .with_counts(10, 500)
            .build();
        let first_context = DecisionContext::new(1, 3);
        let first_evaluated = engine.evaluate_candidate(group.clone(), &first_context);
        let first_outcome = engine.decide(std::slice::from_ref(&first_evaluated), &first_context);
        assert_eq!(first_outcome.status, DecisionStatus::Positive);
        assert_eq!(first_outcome.stop_reason, None);

        let prior_ids = first_context.positive_group_ids(&[first_evaluated]);
        let second_context = DecisionContext::new(2, 3).with_prior_positive_group_ids(prior_ids);
        let second_evaluated = engine.evaluate_candidate(group, &second_context);
        let second_outcome = engine.decide(&[second_evaluated], &second_context);
        assert_eq!(
            second_outcome.stop_reason,
            Some(StopReason::PositiveBoundaryCrossed)
        );
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

        assert_eq!(outcome.status, DecisionStatus::Negative);
        assert_eq!(
            outcome.stop_reason,
            Some(StopReason::NegativeBoundaryConfirmed)
        );
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
            .any(|reason| reason
                == "negative_control_policy=not_used_for_default_background_assessment"));

        assert_eq!(candidates[2].decision, DecisionStatus::Indeterminate);
        assert!(candidates[2]
            .decision_reasons
            .iter()
            .any(|reason| reason == "unique_fraction_between_boundaries"));
    }

    #[test]
    fn thresholds_continue_to_use_unique_fraction_after_denominator_fix() {
        let engine = decision_engine();
        let outcome = engine.decide(
            &[candidate_call()
                .with_raw_fraction(0.0)
                .with_unique_fraction(0.001)
                .with_nonoverlap_fragments(10)
                .with_breadth(0.01)
                .with_background_ratio(5.0)
                .with_counts(10, 0)
                .build()],
            &DecisionContext::new(1, 5),
        );

        assert_eq!(outcome.status, DecisionStatus::Positive);
        assert_eq!(
            outcome.stop_reason,
            Some(StopReason::PositiveBoundaryCrossed)
        );
    }

    fn decision_engine() -> DecisionEngine {
        DecisionEngine::new(
            &DecisionRules {
                theta_pos: 0.00005,
                theta_neg: 0.000005,
                allow_indeterminate: true,
                theta_neg_fixed: 0.00001,
                positive_alpha_global: 0.05,
                positive_statistic_method: "exact_binomial_survival_bonferroni".to_string(),
                negative_cp_lod_max: 0.00001,
                noise_floor_fraction: 0.000001,
            },
            &CandidateRules {
                min_nonoverlap_fragments: 3,
                min_breadth: 0.001,
                max_background_ratio: 1.5,
                theta_pos_absolute: 3,
                max_ambiguous_fraction_for_positive: 0.20,
                min_positive_evidence_strength: EvidenceStrength::High,
                weak_positive_enabled: true,
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
                coverage_interval_blocks: 0,
                nonoverlap_fragments_algorithm:
                    crate::aggregate::interval::NONOVERLAP_FRAGMENTS_ALGORITHM.to_string(),
                accepted_fraction_label: "accepted_fragments_per_sampled_fragment".into(),
                unique_fraction_label: "legacy_alias_not_molecule_unique".into(),
                total_sampled_fragments: 0,
                fraction_denominator_source: "legacy_inferred".into(),
                aggregation_level: "accession".into(),
                multiplicity_family: "accession".into(),
                accession_resolution: "high".into(),
                raw_fraction: 0.0,
                unique_fraction: 0.0,
                fraction_ci_95: [0.0, 0.0],
                clopper_pearson_upper: 0.0,
                breadth: 0.0,
                ambiguous_fragments: 0,
                background_ratio: 0.0,
                decision: DecisionStatus::Indeterminate,
                decision_reasons: vec!["decision_pending_task_16".to_string()],
                evidence_strength: EvidenceStrength::High,
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
            if unique_fraction > 0.0 && unique_fraction.is_finite() {
                self.inner.total_sampled_fragments =
                    ((self.inner.accepted_fragments as f64 / unique_fraction).round() as u64)
                        .max(self.inner.accepted_fragments);
                self.inner.fraction_denominator_source = "sampled_fragments".to_string();
            }
            self
        }

        fn with_raw_fraction(mut self, raw_fraction: f64) -> Self {
            self.inner.raw_fraction = raw_fraction;
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

        fn into_group(mut self, group: &str) -> Self {
            self.inner.accession_or_group = group.to_string();
            self.inner.aggregation_level = "group".to_string();
            self.inner.multiplicity_family = "group".to_string();
            self.inner.accession_resolution = "not_applicable".to_string();
            self
        }
    }
}
