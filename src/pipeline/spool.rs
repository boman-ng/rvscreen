use crate::aggregate::CandidateAggregator;
use crate::align::FragmentAlignResult;
use crate::error::{Result, RvScreenError};
use crate::types::FragmentClass;

/// Bounded representative coordination state.
///
/// Each processed retained-sample fragment is routed exactly once into the
/// representative round where the sorted-sample prefix expands to include it.
/// Later representative rounds are reconstructed by cumulatively merging these
/// per-prefix-expansion deltas after the retained sample has been aligned in
/// stable `(hash, fragment_key)` order.
#[derive(Debug)]
pub(crate) struct RepresentativeRoundSpool {
    delta_sampled_fragments: Vec<u64>,
    delta_aggregators: Vec<CandidateAggregator>,
    next_batch_sequence: usize,
}

#[derive(Debug)]
pub(crate) struct RepresentativeBatchDelta {
    batch_sequence: usize,
    round_deltas: Vec<RepresentativeRoundContribution>,
}

#[derive(Debug)]
pub(crate) struct RepresentativeRoundContribution {
    pub sampled_fragments: u64,
    pub aggregator: CandidateAggregator,
}

impl RepresentativeRoundSpool {
    pub(crate) fn new(prefix_sample_sizes: Vec<u64>) -> Result<Self> {
        if prefix_sample_sizes.is_empty() {
            return Err(RvScreenError::validation(
                "sampling.rounds",
                "representative single-pass coordination requires at least one round",
            ));
        }

        let round_count = prefix_sample_sizes.len();
        let mut previous = 0u64;
        for (round_index, &size) in prefix_sample_sizes.iter().enumerate() {
            if size < previous {
                return Err(RvScreenError::validation(
                    "sampling.rounds",
                    format!(
                        "representative prefix sizes must be non-decreasing, but round {} size {} is smaller than prior size {}",
                        round_index + 1,
                        size,
                        previous
                    ),
                ));
            }
            previous = size;
        }

        Ok(Self {
            delta_sampled_fragments: vec![0; round_count],
            delta_aggregators: vec![CandidateAggregator::new(); round_count],
            next_batch_sequence: 0,
        })
    }

    pub(crate) fn round_count(&self) -> usize {
        self.delta_aggregators.len()
    }

    pub(crate) fn entry_round_for_prefix_index(
        prefix_sample_sizes: &[u64],
        prefix_index: u64,
    ) -> Option<usize> {
        let one_based_index = prefix_index.saturating_add(1);
        prefix_sample_sizes
            .iter()
            .position(|size| one_based_index <= *size)
    }

    pub(crate) fn merge_batch_delta(
        &mut self,
        batch_delta: RepresentativeBatchDelta,
    ) -> Result<()> {
        if batch_delta.batch_sequence != self.next_batch_sequence {
            return Err(RvScreenError::validation(
                "pipeline.parallel",
                format!(
                    "representative batch merge expected sequence {}, but received {}",
                    self.next_batch_sequence, batch_delta.batch_sequence
                ),
            ));
        }

        if batch_delta.round_deltas.len() != self.delta_aggregators.len() {
            return Err(RvScreenError::validation(
                "sampling.rounds",
                format!(
                    "representative batch sequence {} carried {} round delta(s) for {} configured round(s)",
                    batch_delta.batch_sequence,
                    batch_delta.round_deltas.len(),
                    self.delta_aggregators.len()
                ),
            ));
        }

        for (round_index, contribution) in batch_delta.round_deltas.into_iter().enumerate() {
            self.delta_sampled_fragments[round_index] = self.delta_sampled_fragments[round_index]
                .saturating_add(contribution.sampled_fragments);
            self.delta_aggregators[round_index].merge_from(contribution.aggregator)?;
        }

        self.next_batch_sequence = self.next_batch_sequence.saturating_add(1);
        Ok(())
    }

    pub(crate) fn into_round_contributions(self) -> Vec<RepresentativeRoundContribution> {
        self.delta_sampled_fragments
            .into_iter()
            .zip(self.delta_aggregators)
            .map(
                |(sampled_fragments, aggregator)| RepresentativeRoundContribution {
                    sampled_fragments,
                    aggregator,
                },
            )
            .collect()
    }
}

impl RepresentativeBatchDelta {
    pub(crate) fn new(batch_sequence: usize, round_count: usize) -> Self {
        Self {
            batch_sequence,
            round_deltas: (0..round_count)
                .map(|_| RepresentativeRoundContribution {
                    sampled_fragments: 0,
                    aggregator: CandidateAggregator::new(),
                })
                .collect(),
        }
    }

    pub(crate) fn record_processed_fragment(
        &mut self,
        entry_round: usize,
        class: FragmentClass,
        align_result: &FragmentAlignResult,
    ) -> Result<()> {
        let round_count = self.round_deltas.len();
        if entry_round >= round_count {
            return Err(RvScreenError::validation(
                "sampling.rounds",
                format!(
                    "representative batch sequence {} received invalid round index {} for {} configured round(s)",
                    self.batch_sequence,
                    entry_round,
                    round_count
                ),
            ));
        }

        let contribution = self
            .round_deltas
            .get_mut(entry_round)
            .expect("validated representative batch round index should exist");

        contribution.sampled_fragments = contribution.sampled_fragments.saturating_add(1);
        contribution.aggregator.add_fragment(class, align_result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn representative_round_spool_merges_batch_deltas_in_batch_sequence_order() {
        let mut spool = RepresentativeRoundSpool::new(vec![2, 4])
            .expect("representative spool should build for valid rounds");
        let mut batch_zero = RepresentativeBatchDelta::new(0, 2);
        batch_zero.round_deltas[0].sampled_fragments = 1;
        batch_zero.round_deltas[1].sampled_fragments = 2;
        let mut batch_one = RepresentativeBatchDelta::new(1, 2);
        batch_one.round_deltas[1].sampled_fragments = 1;

        spool
            .merge_batch_delta(batch_zero)
            .expect("first batch should merge");
        spool
            .merge_batch_delta(batch_one)
            .expect("later batch should merge after the first one");

        let contributions = spool.into_round_contributions();
        assert_eq!(contributions[0].sampled_fragments, 1);
        assert_eq!(contributions[1].sampled_fragments, 3);
    }

    #[test]
    fn representative_round_spool_rejects_out_of_order_batch_merge_sequence() {
        let mut spool = RepresentativeRoundSpool::new(vec![1, 2])
            .expect("representative spool should build for valid rounds");

        let error = spool
            .merge_batch_delta(RepresentativeBatchDelta::new(1, 2))
            .expect_err("batch sequence gaps should be rejected");

        assert!(error
            .to_string()
            .contains("expected sequence 0, but received 1"));
    }

    #[test]
    fn representative_round_spool_accepts_duplicate_and_zero_prefixes() {
        let spool = RepresentativeRoundSpool::new(vec![0, 10, 10, 25])
            .expect("proportional effective prefixes may be non-decreasing");

        assert_eq!(spool.round_count(), 4);
        assert_eq!(
            RepresentativeRoundSpool::entry_round_for_prefix_index(&[0, 10, 10, 25], 0),
            Some(1)
        );
        assert_eq!(
            RepresentativeRoundSpool::entry_round_for_prefix_index(&[0, 10, 10, 25], 9),
            Some(1)
        );
        assert_eq!(
            RepresentativeRoundSpool::entry_round_for_prefix_index(&[0, 10, 10, 25], 10),
            Some(3)
        );
    }
}
