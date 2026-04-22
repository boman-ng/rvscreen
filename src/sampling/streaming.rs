use crate::sampling::Sampler;

/// Accepts only the first `max_fragments` inputs in stream order.
///
/// No reproducibility guarantee across re-ordered inputs.
/// Callers must treat results as `Provisional` (never `Final`).
pub struct StreamingSampler {
    max_fragments: usize,
}

impl StreamingSampler {
    pub fn new(max_fragments: usize) -> Self {
        Self { max_fragments }
    }

    pub fn max_fragments(&self) -> usize {
        self.max_fragments
    }
}

impl Sampler for StreamingSampler {
    fn should_accept(&self, _fragment_key: &str, index: usize) -> bool {
        index < self.max_fragments
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sampling::Sampler as SamplerTrait;

    #[test]
    fn accepts_exactly_first_n() {
        let sampler = StreamingSampler::new(10);
        let accepted: Vec<usize> = (0..100)
            .filter(|&i| sampler.should_accept("key", i))
            .collect();
        assert_eq!(accepted.len(), 10);
        assert_eq!(accepted, (0..10).collect::<Vec<_>>());
    }

    #[test]
    fn rejects_all_when_max_is_zero() {
        let sampler = StreamingSampler::new(0);
        for i in 0..10 {
            assert!(!sampler.should_accept("key", i));
        }
    }

    #[test]
    fn accepts_all_when_max_exceeds_stream() {
        let sampler = StreamingSampler::new(1000);
        for i in 0..100 {
            assert!(sampler.should_accept("key", i));
        }
    }

    #[test]
    fn fragment_key_is_irrelevant_for_cutoff() {
        let sampler = StreamingSampler::new(3);
        assert!(sampler.should_accept("key-a", 0));
        assert!(sampler.should_accept("key-b", 1));
        assert!(sampler.should_accept("key-c", 2));
        assert!(!sampler.should_accept("key-d", 3));
        assert!(!sampler.should_accept("key-e", 100));
    }

    #[test]
    fn boundary_index_exactly_at_max_is_rejected() {
        let sampler = StreamingSampler::new(5);
        assert!(sampler.should_accept("k", 4));
        assert!(!sampler.should_accept("k", 5));
    }

    #[test]
    fn trait_object_works() {
        let sampler: Box<dyn SamplerTrait> = Box::new(StreamingSampler::new(2));
        assert!(sampler.should_accept("a", 0));
        assert!(sampler.should_accept("b", 1));
        assert!(!sampler.should_accept("c", 2));
    }

    #[test]
    fn max_fragments_accessor() {
        let sampler = StreamingSampler::new(42);
        assert_eq!(sampler.max_fragments(), 42);
    }
}
