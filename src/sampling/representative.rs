use crate::error::{Result, RvScreenError};
use siphasher::sip::SipHasher13;
use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::hash::Hasher;

pub const DEFAULT_REPRESENTATIVE_ROUNDS: [u64; 5] = [50_000, 100_000, 200_000, 400_000, 800_000];
pub const DEFAULT_REPRESENTATIVE_ROUND_PROPORTIONS: [f64; 6] =
    [0.002, 0.005, 0.01, 0.02, 0.05, 0.10];
pub const PROPORTIONAL_MIN_EFFECTIVE_COUNT: u64 = 10_000;
pub const PROPORTIONAL_MAX_EFFECTIVE_COUNT: u64 = 10_000_000;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RepresentativeSamplingWarningKind {
    SmallDatasetRetainedAll,
    EffectiveCountClamped,
    SamplingThresholdOverride,
}

impl RepresentativeSamplingWarningKind {
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::SmallDatasetRetainedAll => "small_dataset_retained_all",
            Self::EffectiveCountClamped => "effective_count_clamped",
            Self::SamplingThresholdOverride => "sampling_threshold_override",
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RepresentativeSamplingWarning {
    pub kind: RepresentativeSamplingWarningKind,
    pub round: Option<usize>,
    pub message: String,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProportionalRoundResolution {
    pub raw_counts: Vec<u64>,
    pub effective_counts: Vec<u64>,
    pub warnings: Vec<RepresentativeSamplingWarning>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RepresentativeSampler {
    seed: u64,
    rounds: Vec<u64>,
    max_fragments: usize,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) struct RepresentativeRank {
    hash: u64,
    fragment_key: String,
}

#[derive(Debug)]
pub(crate) struct RankedRepresentative<T> {
    rank: RepresentativeRank,
    payload: T,
}

#[derive(Debug)]
pub(crate) struct RepresentativeRetention<T> {
    max_fragments: usize,
    retained: BinaryHeap<RankedRepresentative<T>>,
}

impl RepresentativeSampler {
    pub fn new(seed: u64, rounds: Vec<u64>) -> Result<Self> {
        validate_rounds(&rounds)?;
        let max_fragments = usize::try_from(
            *rounds
                .last()
                .expect("validated representative rounds must be non-empty"),
        )
        .map_err(|_| {
            RvScreenError::validation(
                "sampling.rounds",
                "representative sampling final round exceeds platform usize capacity",
            )
        })?;

        Ok(Self {
            seed,
            rounds,
            max_fragments,
        })
    }

    pub fn with_default_rounds(seed: u64) -> Self {
        Self::new(seed, DEFAULT_REPRESENTATIVE_ROUNDS.to_vec())
            .expect("default representative rounds must stay valid")
    }

    pub(crate) fn for_retention_capacity(seed: u64, max_fragments: usize) -> Self {
        Self {
            seed,
            rounds: vec![u64::try_from(max_fragments).unwrap_or(u64::MAX)],
            max_fragments,
        }
    }

    pub fn seed(&self) -> u64 {
        self.seed
    }

    pub fn rounds(&self) -> &[u64] {
        &self.rounds
    }

    pub fn round_count(&self) -> usize {
        self.rounds.len()
    }

    pub fn round_target(&self, round: usize) -> Option<u64> {
        self.rounds.get(round).copied()
    }

    pub fn max_fragments(&self) -> usize {
        self.max_fragments
    }

    pub(crate) fn rank_fragment(&self, fragment_key: &str) -> RepresentativeRank {
        RepresentativeRank {
            hash: self.hash_fragment_key(fragment_key),
            fragment_key: fragment_key.to_string(),
        }
    }

    pub(crate) fn prefix_sample_sizes(&self, retained_sample_len: usize) -> Vec<u64> {
        let retained_sample_len = u64::try_from(retained_sample_len).unwrap_or(u64::MAX);
        self.rounds
            .iter()
            .map(|target| (*target).min(retained_sample_len))
            .collect()
    }

    pub(crate) fn hash_fragment_key(&self, fragment_key: &str) -> u64 {
        let mut hasher = SipHasher13::new_with_keys(0, 0);
        hasher.write_u64(self.seed);
        hasher.write(fragment_key.as_bytes());
        hasher.finish()
    }
}

impl<T> RankedRepresentative<T> {
    pub(crate) fn new(rank: RepresentativeRank, payload: T) -> Self {
        Self { rank, payload }
    }

    pub(crate) fn rank(&self) -> &RepresentativeRank {
        &self.rank
    }

    pub(crate) fn into_payload(self) -> T {
        self.payload
    }
}

impl<T> PartialEq for RankedRepresentative<T> {
    fn eq(&self, other: &Self) -> bool {
        self.rank == other.rank
    }
}

impl<T> Eq for RankedRepresentative<T> {}

impl<T> PartialOrd for RankedRepresentative<T> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<T> Ord for RankedRepresentative<T> {
    fn cmp(&self, other: &Self) -> Ordering {
        self.rank.cmp(&other.rank)
    }
}

impl<T> RepresentativeRetention<T> {
    pub(crate) fn new(max_fragments: usize) -> Self {
        Self {
            max_fragments,
            retained: BinaryHeap::new(),
        }
    }

    pub(crate) fn would_retain_hash(&self, candidate_hash: u64, fragment_key: &str) -> bool {
        if self.max_fragments == 0 {
            return false;
        }

        self.retained.len() < self.max_fragments
            || self.retained.peek().is_some_and(|worst_retained| {
                candidate_hash < worst_retained.rank.hash
                    || (candidate_hash == worst_retained.rank.hash
                        && fragment_key < worst_retained.rank.fragment_key.as_str())
            })
    }

    pub(crate) fn consider(&mut self, candidate: RankedRepresentative<T>) {
        if self.max_fragments == 0 {
            return;
        }

        if self.retained.len() < self.max_fragments {
            self.retained.push(candidate);
            return;
        }

        let should_replace = self
            .retained
            .peek()
            .is_some_and(|worst_retained| candidate.rank() < worst_retained.rank());
        if should_replace {
            self.retained.pop();
            self.retained.push(candidate);
        }
    }

    pub(crate) fn into_sorted_vec(self) -> Vec<RankedRepresentative<T>> {
        let mut retained = self.retained.into_vec();
        retained.sort_by(|left, right| left.rank.cmp(&right.rank));
        retained
    }
}

pub fn resolve_proportional_round_targets(
    proportions: &[f64],
    qc_passing_fragments: u64,
    allow_sampling_threshold_override: bool,
) -> ProportionalRoundResolution {
    let raw_counts = proportions
        .iter()
        .copied()
        .map(|proportion| rounded_proportional_count(proportion, qc_passing_fragments))
        .collect::<Vec<_>>();

    if allow_sampling_threshold_override {
        let mut warnings = Vec::new();
        let effective_counts = raw_counts
            .iter()
            .copied()
            .enumerate()
            .map(|(index, raw_count)| {
                let effective_count = raw_count.min(qc_passing_fragments);
                if !(PROPORTIONAL_MIN_EFFECTIVE_COUNT..=PROPORTIONAL_MAX_EFFECTIVE_COUNT)
                    .contains(&raw_count)
                {
                    warnings.push(RepresentativeSamplingWarning {
                        kind: RepresentativeSamplingWarningKind::SamplingThresholdOverride,
                        round: Some(index + 1),
                        message: format!(
                            "representative proportional round {} raw count {} is outside default sampling thresholds [{}, {}]; --allow-sampling-threshold-override retained effective count {}",
                            index + 1,
                            raw_count,
                            PROPORTIONAL_MIN_EFFECTIVE_COUNT,
                            PROPORTIONAL_MAX_EFFECTIVE_COUNT,
                            effective_count
                        ),
                    });
                }
                effective_count
            })
            .collect();

        return ProportionalRoundResolution {
            raw_counts,
            effective_counts,
            warnings,
        };
    }

    if qc_passing_fragments < PROPORTIONAL_MIN_EFFECTIVE_COUNT {
        return ProportionalRoundResolution {
            effective_counts: vec![qc_passing_fragments; proportions.len()],
            raw_counts,
            warnings: vec![RepresentativeSamplingWarning {
                kind: RepresentativeSamplingWarningKind::SmallDatasetRetainedAll,
                round: None,
                message: format!(
                    "representative proportional sampling observed {} QC-passing fragments, below the {} minimum; retaining all QC-passing fragments for every round",
                    qc_passing_fragments, PROPORTIONAL_MIN_EFFECTIVE_COUNT
                ),
            }],
        };
    }

    let mut warnings = Vec::new();
    let effective_counts = raw_counts
        .iter()
        .copied()
        .enumerate()
        .map(|(index, raw_count)| {
            let clamped = raw_count.clamp(
                PROPORTIONAL_MIN_EFFECTIVE_COUNT,
                PROPORTIONAL_MAX_EFFECTIVE_COUNT,
            );
            let effective_count = clamped.min(qc_passing_fragments);
            if raw_count != clamped {
                warnings.push(RepresentativeSamplingWarning {
                    kind: RepresentativeSamplingWarningKind::EffectiveCountClamped,
                    round: Some(index + 1),
                    message: format!(
                        "representative proportional round {} raw count {} was clamped to effective count {} using sampling thresholds [{}, {}]",
                        index + 1,
                        raw_count,
                        effective_count,
                        PROPORTIONAL_MIN_EFFECTIVE_COUNT,
                        PROPORTIONAL_MAX_EFFECTIVE_COUNT
                    ),
                });
            }
            effective_count
        })
        .collect();

    ProportionalRoundResolution {
        raw_counts,
        effective_counts,
        warnings,
    }
}

fn rounded_proportional_count(proportion: f64, qc_passing_fragments: u64) -> u64 {
    (proportion * qc_passing_fragments as f64).round() as u64
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RoundManager {
    rounds: Vec<u64>,
    current_round: Option<usize>,
    current_round_completed: bool,
}

impl RoundManager {
    pub fn new(rounds: Vec<u64>) -> Result<Self> {
        validate_rounds(&rounds)?;

        Ok(Self {
            rounds,
            current_round: None,
            current_round_completed: false,
        })
    }

    pub fn with_default_rounds() -> Self {
        Self::new(DEFAULT_REPRESENTATIVE_ROUNDS.to_vec())
            .expect("default representative rounds must stay valid")
    }

    pub fn rounds(&self) -> &[u64] {
        &self.rounds
    }

    pub fn current_round(&self) -> Option<usize> {
        self.current_round
    }

    pub fn current_target(&self) -> Option<u64> {
        self.current_round
            .and_then(|round| self.round_target(round))
    }

    pub fn round_target(&self, round: usize) -> Option<u64> {
        self.rounds.get(round).copied()
    }

    pub fn advance_round(&mut self) -> Option<u64> {
        if self.current_round.is_some() && !self.current_round_completed {
            return None;
        }

        let next_round = match self.current_round {
            Some(round) => round + 1,
            None => 0,
        };

        let next_target = self.round_target(next_round)?;
        self.current_round = Some(next_round);
        self.current_round_completed = false;
        Some(next_target)
    }

    pub fn mark_current_round_complete(&mut self) -> bool {
        if self.current_round.is_none() || self.current_round_completed {
            return false;
        }

        self.current_round_completed = true;
        true
    }

    pub fn is_current_round_complete(&self) -> bool {
        self.current_round_completed
    }

    pub fn is_complete(&self) -> bool {
        matches!(
            self.current_round,
            Some(round)
                if self.current_round_completed && round + 1 == self.rounds.len()
        )
    }
}

fn validate_rounds(rounds: &[u64]) -> Result<()> {
    if rounds.is_empty() {
        return Err(RvScreenError::validation(
            "sampling.rounds",
            "representative sampling requires at least one round",
        ));
    }

    if rounds[0] == 0 {
        return Err(RvScreenError::validation(
            "sampling.rounds",
            "representative sampling rounds must be positive",
        ));
    }

    for window in rounds.windows(2) {
        if window[1] <= window[0] {
            return Err(RvScreenError::validation(
                "sampling.rounds",
                "representative sampling rounds must be strictly increasing",
            ));
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{
        resolve_proportional_round_targets, RankedRepresentative, RepresentativeRetention,
        RepresentativeSampler, RepresentativeSamplingWarningKind, RoundManager,
        DEFAULT_REPRESENTATIVE_ROUNDS, DEFAULT_REPRESENTATIVE_ROUND_PROPORTIONS,
    };
    use siphasher::sip::SipHasher13;
    use std::collections::BTreeSet;
    use std::hash::Hasher;

    #[test]
    fn test_representative_sampling_reproducible() {
        let sampler = RepresentativeSampler::with_default_rounds(42);
        let repeated = RepresentativeSampler::with_default_rounds(42);

        let first_pass: Vec<_> = (0..1_000)
            .map(fragment_key)
            .map(|key| sampler.rank_fragment(&key))
            .collect();
        let second_pass: Vec<_> = (0..1_000)
            .map(fragment_key)
            .map(|key| sampler.rank_fragment(&key))
            .collect();
        let recreated_pass: Vec<_> = (0..1_000)
            .map(fragment_key)
            .map(|key| repeated.rank_fragment(&key))
            .collect();

        assert_eq!(first_pass, second_pass);
        assert_eq!(first_pass, recreated_pass);
    }

    #[test]
    fn test_default_representative_round_proportions_are_task5_defaults() {
        assert_eq!(
            DEFAULT_REPRESENTATIVE_ROUND_PROPORTIONS,
            [0.002, 0.005, 0.01, 0.02, 0.05, 0.10]
        );
    }

    #[test]
    fn test_representative_sampling_rounds_are_strictly_nested() {
        let sampler = RepresentativeSampler::with_default_rounds(42);
        let retained_keys = retained_keys_for(&sampler, (0..1_000_000).map(fragment_key).collect());
        let prefix_sizes = sampler.prefix_sample_sizes(retained_keys.len());

        let accepted_by_round: Vec<BTreeSet<String>> = prefix_sizes
            .iter()
            .map(|target| {
                retained_keys
                    .iter()
                    .take(*target as usize)
                    .cloned()
                    .collect()
            })
            .collect();

        for pair in accepted_by_round.windows(2) {
            assert!(pair[0].is_subset(&pair[1]));
        }

        assert_eq!(prefix_sizes, sampler.rounds());
    }

    #[test]
    fn test_representative_exact_k_sampling_caps_final_round_at_kmax() {
        let sampler = RepresentativeSampler::new(42, vec![3, 5, 8])
            .expect("representative rounds should be valid");
        let retained_keys = retained_keys_for(&sampler, (0..128).map(fragment_key).collect());
        let prefix_sizes = sampler.prefix_sample_sizes(retained_keys.len());
        let expected = expected_top_k_keys(&sampler, 128);

        assert_eq!(retained_keys, expected);
        assert_eq!(retained_keys.len(), 8);
        assert_eq!(prefix_sizes, vec![3, 5, 8]);
    }

    #[test]
    fn test_representative_sampling_is_order_independent() {
        let sampler = RepresentativeSampler::with_default_rounds(42);
        let keys: Vec<String> = (0..1_000).map(fragment_key).collect();

        let forward = retained_keys_for(&sampler, keys.clone());
        let reverse = retained_keys_for(&sampler, keys.into_iter().rev().collect());
        assert_eq!(forward, reverse);
    }

    #[test]
    fn test_representative_exact_k_small_input_saturates_all_round_prefixes() {
        let sampler = RepresentativeSampler::new(42, vec![50, 100])
            .expect("representative rounds should be valid");
        let retained_keys = retained_keys_for(&sampler, (0..12).map(fragment_key).collect());

        assert_eq!(retained_keys.len(), 12);
        assert_eq!(
            sampler.prefix_sample_sizes(retained_keys.len()),
            vec![12, 12]
        );
    }

    #[test]
    fn test_round_manager_tracks_progression_and_completion() {
        let mut manager = RoundManager::with_default_rounds();

        assert_eq!(manager.rounds(), &DEFAULT_REPRESENTATIVE_ROUNDS);
        assert_eq!(manager.current_round(), None);
        assert!(!manager.is_complete());

        assert_eq!(manager.advance_round(), Some(50_000));
        assert_eq!(manager.current_round(), Some(0));
        assert_eq!(manager.current_target(), Some(50_000));
        assert!(!manager.is_current_round_complete());
        assert_eq!(manager.advance_round(), None);

        for expected_target in [50_000, 100_000, 200_000, 400_000, 800_000] {
            if manager.current_target() != Some(expected_target) {
                assert_eq!(manager.advance_round(), Some(expected_target));
            }

            assert!(manager.mark_current_round_complete());
            if expected_target != 800_000 {
                assert!(!manager.is_complete());
                assert_eq!(manager.advance_round(), Some(expected_target * 2));
            }
        }

        assert!(manager.is_complete());
        assert!(!manager.mark_current_round_complete());
        assert_eq!(manager.advance_round(), None);
    }

    #[test]
    fn test_invalid_rounds_are_rejected() {
        let err = RepresentativeSampler::new(7, vec![50_000, 50_000])
            .expect_err("duplicate rounds should fail validation");
        assert!(err.to_string().contains("strictly increasing"));

        let err = RoundManager::new(Vec::new()).expect_err("empty rounds should fail validation");
        assert!(err.to_string().contains("at least one round"));
    }

    #[test]
    fn test_proportional_resolution_retains_all_small_datasets() {
        let resolution = resolve_proportional_round_targets(&[0.25, 0.5, 1.0], 9_999, false);

        assert_eq!(resolution.raw_counts, vec![2_500, 5_000, 9_999]);
        assert_eq!(resolution.effective_counts, vec![9_999, 9_999, 9_999]);
        assert_eq!(resolution.warnings.len(), 1);
        assert_eq!(
            resolution.warnings[0].kind,
            RepresentativeSamplingWarningKind::SmallDatasetRetainedAll
        );
        assert!(resolution.warnings[0]
            .message
            .contains("retaining all QC-passing fragments"));
    }

    #[test]
    fn test_proportional_resolution_clamps_to_default_thresholds() {
        let resolution = resolve_proportional_round_targets(&[0.1, 0.75], 20_000_001, false);

        assert_eq!(resolution.raw_counts, vec![2_000_000, 15_000_001]);
        assert_eq!(resolution.effective_counts, vec![2_000_000, 10_000_000]);
        assert_eq!(resolution.warnings.len(), 1);
        assert_eq!(
            resolution.warnings[0].kind,
            RepresentativeSamplingWarningKind::EffectiveCountClamped
        );
        assert_eq!(resolution.warnings[0].round, Some(2));
    }

    #[test]
    fn test_proportional_resolution_clamps_low_raw_counts_upward() {
        let resolution = resolve_proportional_round_targets(&[0.01, 0.5], 20_000, false);

        assert_eq!(resolution.raw_counts, vec![200, 10_000]);
        assert_eq!(resolution.effective_counts, vec![10_000, 10_000]);
        assert_eq!(resolution.warnings.len(), 1);
        assert_eq!(
            resolution.warnings[0].kind,
            RepresentativeSamplingWarningKind::EffectiveCountClamped
        );
        assert_eq!(resolution.warnings[0].round, Some(1));
    }

    #[test]
    fn test_proportional_resolution_allow_override_uses_raw_counts() {
        let resolution = resolve_proportional_round_targets(&[0.0001, 1.0], 20_000_001, true);

        assert_eq!(resolution.raw_counts, vec![2_000, 20_000_001]);
        assert_eq!(resolution.effective_counts, vec![2_000, 20_000_001]);
        assert_eq!(resolution.warnings.len(), 2);
        assert!(resolution.warnings.iter().all(|warning| {
            warning.kind == RepresentativeSamplingWarningKind::SamplingThresholdOverride
                && warning
                    .message
                    .contains("--allow-sampling-threshold-override")
        }));
    }

    #[test]
    fn task8_threshold_resolution_covers_below_10k_and_above_10m_without_fixtures() {
        let below_minimum = resolve_proportional_round_targets(&[0.5, 1.0], 9_999, false);
        assert_eq!(below_minimum.raw_counts, vec![5_000, 9_999]);
        assert_eq!(below_minimum.effective_counts, vec![9_999, 9_999]);
        assert_eq!(below_minimum.warnings.len(), 1);
        assert_eq!(
            below_minimum.warnings[0].kind,
            RepresentativeSamplingWarningKind::SmallDatasetRetainedAll
        );

        let above_maximum = resolve_proportional_round_targets(&[0.5, 1.0], 20_000_000, false);
        assert_eq!(above_maximum.raw_counts, vec![10_000_000, 20_000_000]);
        assert_eq!(above_maximum.effective_counts, vec![10_000_000, 10_000_000]);
        assert_eq!(above_maximum.warnings.len(), 1);
        assert_eq!(
            above_maximum.warnings[0].kind,
            RepresentativeSamplingWarningKind::EffectiveCountClamped
        );
        assert_eq!(above_maximum.warnings[0].round, Some(2));

        let override_resolution =
            resolve_proportional_round_targets(&[0.0001, 1.0], 20_000_001, true);
        assert_eq!(override_resolution.raw_counts, vec![2_000, 20_000_001]);
        assert_eq!(
            override_resolution.effective_counts,
            vec![2_000, 20_000_001]
        );
        assert_eq!(override_resolution.warnings.len(), 2);
        assert!(override_resolution.warnings.iter().all(|warning| {
            warning.kind == RepresentativeSamplingWarningKind::SamplingThresholdOverride
        }));
    }

    #[test]
    fn test_proportional_retained_fragment_keys_are_repeatable_and_nested() {
        let seed = 20_260_420;
        let proportions = [0.25, 0.5, 1.0];
        let keys = (0..24).map(fragment_key).collect::<Vec<_>>();

        let run_a = proportional_retained_key_prefixes(seed, &keys, &proportions, true);
        let run_b = proportional_retained_key_prefixes(seed, &keys, &proportions, true);

        assert_eq!(run_a, run_b, "repeated runs must retain identical IDs");
        assert_eq!(
            run_a.iter().map(Vec::len).collect::<Vec<_>>(),
            vec![6, 12, 24]
        );
        for pair in run_a.windows(2) {
            let smaller = pair[0].iter().cloned().collect::<BTreeSet<_>>();
            let larger = pair[1].iter().cloned().collect::<BTreeSet<_>>();
            assert!(
                smaller.is_subset(&larger),
                "smaller proportional retained-ID set must be nested in larger set"
            );
        }

        print_retained_key_evidence("run_a", &run_a);
        print_retained_key_evidence("run_b", &run_b);
        println!("retained_id_equality run_a_vs_run_b=true");
        println!("retained_id_subset round1_in_round2=true round2_in_round3=true");
    }

    fn fragment_key(index: usize) -> String {
        format!("fragment-{index:06}")
    }

    fn retained_keys_for(sampler: &RepresentativeSampler, keys: Vec<String>) -> Vec<String> {
        let mut retention = RepresentativeRetention::new(sampler.max_fragments());
        for key in keys {
            retention.consider(RankedRepresentative::new(sampler.rank_fragment(&key), key));
        }

        retention
            .into_sorted_vec()
            .into_iter()
            .map(RankedRepresentative::into_payload)
            .collect()
    }

    fn proportional_retained_key_prefixes(
        seed: u64,
        keys: &[String],
        proportions: &[f64],
        allow_sampling_threshold_override: bool,
    ) -> Vec<Vec<String>> {
        let max_fragments = if allow_sampling_threshold_override {
            usize::MAX
        } else {
            super::PROPORTIONAL_MAX_EFFECTIVE_COUNT as usize
        };
        let sampler = RepresentativeSampler::for_retention_capacity(seed, max_fragments);
        let retained_keys = retained_keys_for(&sampler, keys.to_vec());
        let resolution = resolve_proportional_round_targets(
            proportions,
            keys.len() as u64,
            allow_sampling_threshold_override,
        );

        resolution
            .effective_counts
            .into_iter()
            .map(|prefix_size| {
                retained_keys
                    .iter()
                    .take(prefix_size as usize)
                    .cloned()
                    .collect::<Vec<_>>()
            })
            .collect()
    }

    fn print_retained_key_evidence(label: &str, rounds: &[Vec<String>]) {
        for (round_index, round_keys) in rounds.iter().enumerate() {
            println!(
                "retained_ids {label} round={} count={} checksum={:016x} ids={}",
                round_index + 1,
                round_keys.len(),
                retained_key_checksum(round_keys),
                round_keys.join(",")
            );
        }
    }

    fn retained_key_checksum(keys: &[String]) -> u64 {
        let mut hasher = SipHasher13::new_with_keys(0, 0);
        for key in keys {
            hasher.write_usize(key.len());
            hasher.write(key.as_bytes());
        }
        hasher.finish()
    }

    fn expected_top_k_keys(sampler: &RepresentativeSampler, total_keys: usize) -> Vec<String> {
        let mut ranked = (0..total_keys)
            .map(fragment_key)
            .map(|key| RankedRepresentative::new(sampler.rank_fragment(&key), key))
            .collect::<Vec<_>>();
        ranked.sort();
        ranked
            .into_iter()
            .take(sampler.max_fragments())
            .map(RankedRepresentative::into_payload)
            .collect()
    }
}
