pub mod representative;
pub mod streaming;

pub use representative::{RepresentativeSampler, RoundManager, DEFAULT_REPRESENTATIVE_ROUNDS};
pub use streaming::StreamingSampler;

/// Unified interface for all sampling strategies.
pub trait Sampler {
    /// Returns `true` if the fragment should be accepted.
    ///
    /// `fragment_key`: normalised QNAME.
    /// `index`: 0-based stream position (used by streaming mode; ignored by representative).
    fn should_accept(&self, fragment_key: &str, index: usize) -> bool;
}
