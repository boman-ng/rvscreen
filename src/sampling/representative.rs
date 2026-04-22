use crate::error::{Result, RvScreenError};
use crate::sampling::Sampler;
use siphasher::sip::SipHasher13;
use std::hash::Hasher;

pub const DEFAULT_REPRESENTATIVE_ROUNDS: [u64; 5] = [50_000, 100_000, 200_000, 400_000, 800_000];

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RepresentativeSampler {
    seed: u64,
    rounds: Vec<u64>,
    bucket_universe: u64,
}

impl RepresentativeSampler {
    pub fn new(seed: u64, rounds: Vec<u64>) -> Result<Self> {
        validate_rounds(&rounds)?;

        Ok(Self {
            seed,
            bucket_universe: *rounds
                .last()
                .expect("validated representative rounds must be non-empty"),
            rounds,
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

    pub fn bucket_universe(&self) -> u64 {
        self.bucket_universe
    }

    pub fn bucket_for(&self, fragment_key: &str) -> u64 {
        self.hash_fragment_key(fragment_key) % self.bucket_universe
    }

    pub fn should_accept(&self, fragment_key: &str, round: usize) -> bool {
        self.round_target(round)
            .is_some_and(|target| self.bucket_for(fragment_key) < target)
    }

    fn hash_fragment_key(&self, fragment_key: &str) -> u64 {
        let mut hasher = SipHasher13::new_with_keys(0, 0);
        hasher.write_u64(self.seed);
        hasher.write(fragment_key.as_bytes());
        hasher.finish()
    }
}

impl Sampler for RepresentativeSampler {
    fn should_accept(&self, fragment_key: &str, _index: usize) -> bool {
        let last_round = self.round_count().saturating_sub(1);
        RepresentativeSampler::should_accept(self, fragment_key, last_round)
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
    use super::{RepresentativeSampler, RoundManager, DEFAULT_REPRESENTATIVE_ROUNDS};
    use std::collections::BTreeSet;

    #[test]
    fn test_representative_sampling_reproducible() {
        let sampler = RepresentativeSampler::with_default_rounds(42);
        let repeated = RepresentativeSampler::with_default_rounds(42);

        let first_pass: Vec<bool> = (0..1_000)
            .map(fragment_key)
            .map(|key| sampler.should_accept(&key, 0))
            .collect();
        let second_pass: Vec<bool> = (0..1_000)
            .map(fragment_key)
            .map(|key| sampler.should_accept(&key, 0))
            .collect();
        let recreated_pass: Vec<bool> = (0..1_000)
            .map(fragment_key)
            .map(|key| repeated.should_accept(&key, 0))
            .collect();

        assert_eq!(first_pass, second_pass);
        assert_eq!(first_pass, recreated_pass);
    }

    #[test]
    fn test_representative_sampling_rounds_are_strictly_nested() {
        let sampler = RepresentativeSampler::with_default_rounds(42);
        let keys: Vec<String> = (0..20_000).map(fragment_key).collect();

        let accepted_by_round: Vec<BTreeSet<String>> = (0..sampler.round_count())
            .map(|round| {
                keys.iter()
                    .filter(|key| sampler.should_accept(key, round))
                    .cloned()
                    .collect()
            })
            .collect();

        for pair in accepted_by_round.windows(2) {
            assert!(pair[0].is_subset(&pair[1]));
        }

        for round in 1..sampler.round_count() {
            let witness = (0..200_000).map(fragment_key).find(|key| {
                sampler.should_accept(key, round) && !sampler.should_accept(key, round - 1)
            });

            assert!(
                witness.is_some(),
                "expected round {round} to strictly expand round {}",
                round - 1
            );
        }
    }

    #[test]
    fn test_representative_sampling_is_order_independent() {
        let sampler = RepresentativeSampler::with_default_rounds(42);
        let keys: Vec<String> = (0..1_000).map(fragment_key).collect();

        let forward: BTreeSet<String> = keys
            .iter()
            .filter(|key| sampler.should_accept(key, 2))
            .cloned()
            .collect();
        let reverse: BTreeSet<String> = keys
            .iter()
            .rev()
            .filter(|key| sampler.should_accept(key, 2))
            .cloned()
            .collect();

        assert_eq!(forward, reverse);
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
}
