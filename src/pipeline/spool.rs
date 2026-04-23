use crate::aggregate::CandidateAggregator;
use crate::align::FragmentAlignResult;
use crate::error::{Result, RvScreenError};
use crate::types::FragmentClass;

/// Bounded representative coordination state.
///
/// Each processed fragment is routed exactly once into the first round that
/// accepts its stable hash bucket. Later representative rounds are reconstructed
/// by cumulatively merging these per-entry-round deltas after the single reader
/// pass completes.
#[derive(Debug)]
pub(crate) struct RepresentativeRoundSpool {
    delta_sampled_fragments: Vec<u64>,
    delta_aggregators: Vec<CandidateAggregator>,
}

#[derive(Debug)]
pub(crate) struct RepresentativeRoundContribution {
    pub sampled_fragments: u64,
    pub aggregator: CandidateAggregator,
}

impl RepresentativeRoundSpool {
    pub(crate) fn new(round_count: usize) -> Result<Self> {
        if round_count == 0 {
            return Err(RvScreenError::validation(
                "sampling.rounds",
                "representative single-pass coordination requires at least one round",
            ));
        }

        Ok(Self {
            delta_sampled_fragments: vec![0; round_count],
            delta_aggregators: vec![CandidateAggregator::new(); round_count],
        })
    }

    pub(crate) fn record_processed_fragment(
        &mut self,
        entry_round: usize,
        class: FragmentClass,
        align_result: &FragmentAlignResult,
    ) -> Result<()> {
        if entry_round >= self.delta_aggregators.len() {
            return Err(RvScreenError::validation(
                "sampling.rounds",
                format!(
                    "representative single-pass coordination received invalid round index {} for {} configured round(s)",
                    entry_round,
                    self.delta_aggregators.len()
                ),
            ));
        }

        self.delta_sampled_fragments[entry_round] = self.delta_sampled_fragments[entry_round]
            .saturating_add(1);
        self.delta_aggregators[entry_round].add_fragment(class, align_result)
    }

    pub(crate) fn into_round_contributions(self) -> Vec<RepresentativeRoundContribution> {
        self.delta_sampled_fragments
            .into_iter()
            .zip(self.delta_aggregators)
            .map(|(sampled_fragments, aggregator)| RepresentativeRoundContribution {
                sampled_fragments,
                aggregator,
            })
            .collect()
    }
}
