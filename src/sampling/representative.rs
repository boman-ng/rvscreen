use crate::error::{Result, RvScreenError};
use siphasher::sip::SipHasher13;
use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::hash::Hasher;

pub const DEFAULT_REPRESENTATIVE_ROUNDS: [u64; 5] = [50_000, 100_000, 200_000, 400_000, 800_000];

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
            retained: BinaryHeap::with_capacity(max_fragments),
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
        RankedRepresentative, RepresentativeRetention, RepresentativeSampler, RoundManager,
        DEFAULT_REPRESENTATIVE_ROUNDS,
    };
    use std::collections::BTreeSet;

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
