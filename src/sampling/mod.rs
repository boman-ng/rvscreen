pub mod representative;
pub mod streaming;

pub use representative::{
    resolve_proportional_round_targets, ProportionalRoundResolution, RepresentativeSampler,
    RepresentativeSamplingWarning, RepresentativeSamplingWarningKind, RoundManager,
    DEFAULT_REPRESENTATIVE_ROUNDS, DEFAULT_REPRESENTATIVE_ROUND_PROPORTIONS,
    PROPORTIONAL_MAX_EFFECTIVE_COUNT,
};
pub(crate) use representative::{RankedRepresentative, RepresentativeRetention};
pub use streaming::StreamingSampler;

/// Unified interface for all sampling strategies.
pub trait Sampler {
    /// Returns `true` if the fragment should be accepted.
    ///
    /// `fragment_key`: normalised QNAME.
    /// `index`: 0-based stream position (used by streaming mode; ignored by representative).
    fn should_accept(&self, fragment_key: &str, index: usize) -> bool;
}
