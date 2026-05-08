use rvscreen::decision::{DecisionContext, DecisionEngine};
use rvscreen::types::{
    CandidateCall, CandidateRules, DecisionRules, DecisionStatus, EvidenceStrength, ReleaseStatus,
};
use serde::{Deserialize, Serialize};
use std::fs;
use std::path::Path;

#[derive(Debug, Deserialize)]
struct ReplaySample {
    sample_id: String,
    virus_name: String,
    truth: Truth,
    expected_status: DecisionStatus,
    expected_release_status: ReleaseStatus,
    accepted_fragments: u64,
    nonoverlap_fragments: u64,
    sampled_fragments: u64,
    breadth: f64,
    ambiguous_fragments: u64,
    background_ratio: f64,
    evidence_strength: EvidenceStrength,
    background: Background,
    candidate_family_size: u64,
    multiplicity_family: String,
    planned_round_look_count: u64,
    current_look_index: u64,
}

#[derive(Debug, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
enum Truth {
    Negative,
    Positive,
}

#[derive(Debug, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
enum Background {
    Observed,
    Unknown,
}

#[derive(Debug, Serialize)]
struct ReplayMetrics {
    total: u64,
    negatives: u64,
    positives: u64,
    true_positives: u64,
    false_positives: u64,
    true_negatives: u64,
    false_negatives: u64,
    weak_positive_count: u64,
    indeterminate_count: u64,
    final_count: u64,
    provisional_count: u64,
    blocked_count: u64,
    specificity: f64,
    sensitivity: f64,
    hpv18_like_detected_at_least_weak_positive: bool,
    per_virus_confusion: std::collections::BTreeMap<String, ConfusionCounts>,
}

#[derive(Debug, Default, Serialize)]
struct ConfusionCounts {
    tp: u64,
    fp: u64,
    tn: u64,
    fn_: u64,
}

#[test]
fn clinical_replay_enforces_decision_robustness_gates() {
    let samples = load_samples();
    let engine = DecisionEngine::new(&decision_rules(), &candidate_rules());
    let mut calls = Vec::new();

    for sample in &samples {
        let context = DecisionContext::with_gate_context(
            3,
            sample.planned_round_look_count,
            sample.candidate_family_size,
            sample.planned_round_look_count,
            sample.current_look_index,
        );
        let evaluated =
            engine.evaluate_candidates_with_context(&[candidate_from_sample(sample)], &context);
        let outcome = engine.decide(&evaluated, &context);
        let candidate = &evaluated[0];

        assert!(candidate
            .decision_reasons
            .iter()
            .any(|reason| reason.starts_with("positive_p_value=")));
        assert!(
            candidate
                .decision_reasons
                .iter()
                .any(|reason| reason
                    == &format!("multiplicity_family={}", sample.multiplicity_family))
        );
        assert!(candidate
            .decision_reasons
            .iter()
            .any(|reason| reason.starts_with("total_sampled_fragments=")));
        assert!(candidate
            .decision_reasons
            .iter()
            .any(|reason| reason
                == &format!("candidate_family_size={}", sample.candidate_family_size)));
        assert!(candidate.decision_reasons.iter().any(|reason| reason
            == &format!(
                "planned_round_look_count={}",
                sample.planned_round_look_count
            )));
        assert!(candidate
            .decision_reasons
            .iter()
            .any(|reason| reason == &format!("current_look_index={}", sample.current_look_index)));
        assert_eq!(candidate.decision, sample.expected_status);

        match sample.sample_id.as_str() {
            "single_fragment_5k" | "single_fragment_20k" => {
                assert_ne!(candidate.decision, DecisionStatus::Positive);
                assert_ne!(outcome.status, DecisionStatus::Positive);
            }
            "low_evidence_signal" => {
                assert_eq!(candidate.decision, DecisionStatus::WeakPositive);
            }
            "nonzero_low_noise_negative" => {
                assert_eq!(outcome.status, DecisionStatus::Negative);
            }
            "unknown_background_signal" => {
                assert_eq!(candidate.decision, DecisionStatus::Positive);
                assert!(candidate.decision_reasons.iter().any(|reason| {
                    reason == "negative_control_policy=not_used_for_default_background_assessment"
                }));
            }
            "hpv18_like_stable_low_signal" => {
                assert_eq!(candidate.decision, DecisionStatus::WeakPositive);
            }
            _ => {}
        }

        calls.push((sample, candidate.decision.clone()));
    }

    let metrics = metrics(&calls);
    fs::create_dir_all("target/clinical_replay").expect("metrics directory should be created");
    fs::write(
        "target/clinical_replay/metrics.json",
        serde_json::to_vec_pretty(&metrics).expect("metrics serialize"),
    )
    .expect("metrics should be written");

    assert_eq!(metrics.false_positives, 0);
    assert_eq!(metrics.false_negatives, 0);
    assert!(metrics.hpv18_like_detected_at_least_weak_positive);
}

fn load_samples() -> Vec<ReplaySample> {
    let fixture = Path::new("tests/fixtures/clinical_replay/samples.jsonl");
    fs::read_to_string(fixture)
        .expect("clinical replay fixture should be readable")
        .lines()
        .map(|line| serde_json::from_str(line).expect("clinical replay row should parse"))
        .collect()
}

fn candidate_from_sample(sample: &ReplaySample) -> CandidateCall {
    let unique_fraction = sample.accepted_fragments as f64 / sample.sampled_fragments as f64;
    let mut reasons = vec!["fraction_ci_95_label=sampling-only CI: accepted-fragment sampling uncertainty only; not diagnostic confidence".to_string()];
    if sample.background == Background::Unknown {
        reasons.push("negative_control_background_fraction_absent".to_string());
    }
    CandidateCall {
        virus_name: sample.virus_name.clone(),
        taxid: 18,
        accession_or_group: sample.sample_id.clone(),
        accepted_fragments: sample.accepted_fragments,
        nonoverlap_fragments: sample.nonoverlap_fragments,
        coverage_interval_blocks: 0,
        nonoverlap_fragments_algorithm:
            rvscreen::aggregate::interval::NONOVERLAP_FRAGMENTS_ALGORITHM.to_string(),
        accepted_fraction_label: "accepted_fragments_per_sampled_fragment".into(),
        unique_fraction_label: "legacy_alias_not_molecule_unique".into(),
        total_sampled_fragments: sample.sampled_fragments,
        fraction_denominator_source: "sampled_fragments".into(),
        aggregation_level: sample.multiplicity_family.clone(),
        multiplicity_family: sample.multiplicity_family.clone(),
        accession_resolution: "high".to_string(),
        raw_fraction: unique_fraction,
        unique_fraction,
        fraction_ci_95: [0.0, unique_fraction],
        clopper_pearson_upper: if sample.sample_id == "nonzero_low_noise_negative" {
            0.000002
        } else {
            unique_fraction
        },
        breadth: sample.breadth,
        ambiguous_fragments: sample.ambiguous_fragments,
        background_ratio: sample.background_ratio,
        decision: DecisionStatus::Indeterminate,
        decision_reasons: reasons,
        evidence_strength: sample.evidence_strength.clone(),
    }
}

fn decision_rules() -> DecisionRules {
    DecisionRules {
        theta_pos: 0.00005,
        theta_neg: 0.000005,
        allow_indeterminate: true,
        theta_neg_fixed: 0.00001,
        positive_alpha_global: 0.05,
        positive_statistic_method: "exact_binomial_survival_bonferroni".to_string(),
        negative_cp_lod_max: 0.00001,
        noise_floor_fraction: 0.000001,
    }
}

fn candidate_rules() -> CandidateRules {
    CandidateRules {
        min_nonoverlap_fragments: 3,
        min_breadth: 0.001,
        max_background_ratio: 1.5,
        theta_pos_absolute: 3,
        max_ambiguous_fraction_for_positive: 0.20,
        min_positive_evidence_strength: EvidenceStrength::High,
        weak_positive_enabled: true,
    }
}

fn metrics(calls: &[(&ReplaySample, DecisionStatus)]) -> ReplayMetrics {
    use std::collections::BTreeMap;

    let negatives = calls
        .iter()
        .filter(|(sample, _)| sample.truth == Truth::Negative)
        .count() as u64;
    let positives = calls
        .iter()
        .filter(|(sample, _)| sample.truth == Truth::Positive)
        .count() as u64;
    let true_positives = calls
        .iter()
        .filter(|(sample, status)| {
            sample.truth == Truth::Positive
                && matches!(
                    status,
                    DecisionStatus::Positive | DecisionStatus::WeakPositive
                )
        })
        .count() as u64;
    let false_positives = calls
        .iter()
        .filter(|(sample, status)| {
            sample.truth == Truth::Negative && *status == DecisionStatus::Positive
        })
        .count() as u64;
    let true_negatives = calls
        .iter()
        .filter(|(sample, status)| {
            sample.truth == Truth::Negative && !matches!(status, DecisionStatus::Positive)
        })
        .count() as u64;
    let false_negatives = calls
        .iter()
        .filter(|(sample, status)| {
            sample.truth == Truth::Positive
                && matches!(
                    status,
                    DecisionStatus::Negative | DecisionStatus::Indeterminate
                )
        })
        .count() as u64;
    let weak_positive_count = calls
        .iter()
        .filter(|(_, status)| *status == DecisionStatus::WeakPositive)
        .count() as u64;
    let indeterminate_count = calls
        .iter()
        .filter(|(_, status)| *status == DecisionStatus::Indeterminate)
        .count() as u64;
    let final_count = calls
        .iter()
        .filter(|(sample, _)| sample.expected_release_status == ReleaseStatus::Final)
        .count() as u64;
    let provisional_count = calls
        .iter()
        .filter(|(sample, _)| sample.expected_release_status == ReleaseStatus::Provisional)
        .count() as u64;
    let blocked_count = calls
        .iter()
        .filter(|(sample, _)| sample.expected_release_status == ReleaseStatus::Blocked)
        .count() as u64;
    let mut per_virus_confusion = BTreeMap::<String, ConfusionCounts>::new();
    for (sample, status) in calls {
        let entry = per_virus_confusion
            .entry(sample.virus_name.clone())
            .or_default();
        match (&sample.truth, status) {
            (Truth::Positive, DecisionStatus::Positive | DecisionStatus::WeakPositive) => {
                entry.tp += 1
            }
            (Truth::Positive, _) => entry.fn_ += 1,
            (Truth::Negative, DecisionStatus::Positive) => entry.fp += 1,
            (Truth::Negative, _) => entry.tn += 1,
        }
    }

    ReplayMetrics {
        total: calls.len() as u64,
        negatives,
        positives,
        true_positives,
        false_positives,
        true_negatives,
        false_negatives,
        weak_positive_count,
        indeterminate_count,
        final_count,
        provisional_count,
        blocked_count,
        specificity: 1.0 - false_positives as f64 / negatives as f64,
        sensitivity: 1.0 - false_negatives as f64 / positives as f64,
        hpv18_like_detected_at_least_weak_positive: calls.iter().any(|(sample, status)| {
            sample.sample_id == "hpv18_like_stable_low_signal"
                && matches!(
                    status,
                    DecisionStatus::WeakPositive | DecisionStatus::Positive
                )
        }),
        per_virus_confusion,
    }
}
