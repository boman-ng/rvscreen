mod complexity;

use crate::error::{Result, RvScreenError};

pub use complexity::shannon_entropy_2mer;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Mate {
    Read1,
    Read2,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct QcConfig {
    pub min_read_length: usize,
    pub max_read_length: usize,
    pub max_n_fraction: f64,
    pub low_complexity_entropy_bits: f64,
    pub require_pair: bool,
}

impl Default for QcConfig {
    fn default() -> Self {
        Self {
            min_read_length: 100,
            max_read_length: 250,
            max_n_fraction: 0.10,
            low_complexity_entropy_bits: 1.0,
            require_pair: true,
        }
    }
}

impl QcConfig {
    pub fn validate(&self) -> Result<()> {
        if self.min_read_length == 0 {
            return Err(RvScreenError::validation(
                "qc.min_read_length",
                "must be greater than zero",
            ));
        }

        if self.min_read_length > self.max_read_length {
            return Err(RvScreenError::validation(
                "qc.read_length_range",
                format!(
                    "min_read_length {} cannot exceed max_read_length {}",
                    self.min_read_length, self.max_read_length
                ),
            ));
        }

        if !self.max_n_fraction.is_finite() || !(0.0..=1.0).contains(&self.max_n_fraction) {
            return Err(RvScreenError::validation(
                "qc.max_n_fraction",
                format!(
                    "must be finite and within 0.0..=1.0, got {}",
                    self.max_n_fraction
                ),
            ));
        }

        if !self.low_complexity_entropy_bits.is_finite() || self.low_complexity_entropy_bits < 0.0 {
            return Err(RvScreenError::validation(
                "qc.low_complexity_entropy_bits",
                format!(
                    "must be finite and non-negative, got {}",
                    self.low_complexity_entropy_bits
                ),
            ));
        }

        Ok(())
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FragmentRecord<'a> {
    pub read1: &'a [u8],
    pub read2: Option<&'a [u8]>,
}

impl<'a> FragmentRecord<'a> {
    pub fn paired(read1: &'a [u8], read2: &'a [u8]) -> Self {
        Self {
            read1,
            read2: Some(read2),
        }
    }

    pub fn single(read1: &'a [u8]) -> Self {
        Self { read1, read2: None }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum QcResult {
    Pass(QcPass),
    Filtered(QcFailure),
}

impl QcResult {
    pub fn is_pass(&self) -> bool {
        matches!(self, Self::Pass(_))
    }

    pub fn reason(&self) -> Option<&QcFilterReason> {
        match self {
            Self::Pass(_) => None,
            Self::Filtered(failure) => Some(&failure.reason),
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct QcPass {
    pub n_fraction: f64,
    pub read1_entropy_bits: f64,
    pub read2_entropy_bits: Option<f64>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct QcFailure {
    pub reason: QcFilterReason,
}

#[derive(Debug, Clone, PartialEq)]
pub enum QcFilterReason {
    PairIncomplete {
        missing_or_empty: Mate,
    },
    LengthOutOfBounds {
        mate: Mate,
        observed: usize,
        min: usize,
        max: usize,
    },
    NContent {
        observed_fraction: f64,
        max_fraction: f64,
    },
    LowComplexity {
        mate: Mate,
        observed_entropy_bits: f64,
        min_entropy_bits: f64,
    },
}

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct QcStats {
    pub total_fragments: u64,
    pub passed_fragments: u64,
    pub low_complexity_filtered: u64,
    pub n_content_filtered: u64,
    pub length_filtered: u64,
    pub pair_incomplete: u64,
}

impl QcStats {
    pub fn record(&mut self, result: &QcResult) {
        self.total_fragments += 1;

        match result {
            QcResult::Pass(_) => self.passed_fragments += 1,
            QcResult::Filtered(failure) => match &failure.reason {
                QcFilterReason::PairIncomplete { .. } => self.pair_incomplete += 1,
                QcFilterReason::LengthOutOfBounds { .. } => self.length_filtered += 1,
                QcFilterReason::NContent { .. } => self.n_content_filtered += 1,
                QcFilterReason::LowComplexity { .. } => self.low_complexity_filtered += 1,
            },
        }
    }
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct QcFilter {
    config: QcConfig,
}

impl QcFilter {
    pub fn new(config: QcConfig) -> Result<Self> {
        config.validate()?;
        Ok(Self { config })
    }

    pub fn config(&self) -> &QcConfig {
        &self.config
    }

    pub fn filter(&self, fragment: &FragmentRecord<'_>) -> QcResult {
        if self.config.require_pair && fragment.read1.is_empty() {
            return QcResult::Filtered(QcFailure {
                reason: QcFilterReason::PairIncomplete {
                    missing_or_empty: Mate::Read1,
                },
            });
        }

        if let Some(reason) = self.length_failure(Mate::Read1, fragment.read1) {
            return QcResult::Filtered(QcFailure { reason });
        }

        let read2 = match fragment.read2 {
            Some(read2) if self.config.require_pair && read2.is_empty() => {
                return QcResult::Filtered(QcFailure {
                    reason: QcFilterReason::PairIncomplete {
                        missing_or_empty: Mate::Read2,
                    },
                });
            }
            None if self.config.require_pair => {
                return QcResult::Filtered(QcFailure {
                    reason: QcFilterReason::PairIncomplete {
                        missing_or_empty: Mate::Read2,
                    },
                });
            }
            read2 => read2,
        };

        if let Some(read2) = read2 {
            if let Some(reason) = self.length_failure(Mate::Read2, read2) {
                return QcResult::Filtered(QcFailure { reason });
            }
        }

        let n_fraction = fragment_n_fraction(fragment.read1, read2);
        if n_fraction > self.config.max_n_fraction {
            return QcResult::Filtered(QcFailure {
                reason: QcFilterReason::NContent {
                    observed_fraction: n_fraction,
                    max_fraction: self.config.max_n_fraction,
                },
            });
        }

        let read1_entropy_bits = shannon_entropy_2mer(fragment.read1);
        if read1_entropy_bits < self.config.low_complexity_entropy_bits {
            return QcResult::Filtered(QcFailure {
                reason: QcFilterReason::LowComplexity {
                    mate: Mate::Read1,
                    observed_entropy_bits: read1_entropy_bits,
                    min_entropy_bits: self.config.low_complexity_entropy_bits,
                },
            });
        }

        let read2_entropy_bits = read2.map(shannon_entropy_2mer);
        if let Some(read2_entropy_bits) = read2_entropy_bits {
            if read2_entropy_bits < self.config.low_complexity_entropy_bits {
                return QcResult::Filtered(QcFailure {
                    reason: QcFilterReason::LowComplexity {
                        mate: Mate::Read2,
                        observed_entropy_bits: read2_entropy_bits,
                        min_entropy_bits: self.config.low_complexity_entropy_bits,
                    },
                });
            }
        }

        QcResult::Pass(QcPass {
            n_fraction,
            read1_entropy_bits,
            read2_entropy_bits,
        })
    }

    pub fn filter_and_record(
        &self,
        fragment: &FragmentRecord<'_>,
        stats: &mut QcStats,
    ) -> QcResult {
        let result = self.filter(fragment);
        stats.record(&result);
        result
    }

    fn length_failure(&self, mate: Mate, read: &[u8]) -> Option<QcFilterReason> {
        let observed = read.len();
        if observed < self.config.min_read_length || observed > self.config.max_read_length {
            return Some(QcFilterReason::LengthOutOfBounds {
                mate,
                observed,
                min: self.config.min_read_length,
                max: self.config.max_read_length,
            });
        }

        None
    }
}

fn fragment_n_fraction(read1: &[u8], read2: Option<&[u8]>) -> f64 {
    let total_bases = read1.len() + read2.map_or(0, |read| read.len());
    let n_bases = count_ns(read1) + read2.map_or(0, count_ns);

    n_bases as f64 / total_bases as f64
}

fn count_ns(read: &[u8]) -> usize {
    read.iter()
        .filter(|base| matches!(base.to_ascii_uppercase(), b'N'))
        .count()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn normal_fragment_passes_qc() {
        let filter = QcFilter::default();
        let read1 = diverse_sequence(150, 7);
        let read2 = diverse_sequence(150, 19);

        let result = filter.filter(&FragmentRecord::paired(read1.as_bytes(), read2.as_bytes()));

        match result {
            QcResult::Pass(metrics) => {
                assert!(
                    metrics.n_fraction < 0.10,
                    "expected low N fraction, got {metrics:?}"
                );
                assert!(
                    metrics.read1_entropy_bits > filter.config().low_complexity_entropy_bits,
                    "expected read1 to remain above low-complexity threshold, got {metrics:?}"
                );
                assert!(
                    metrics
                        .read2_entropy_bits
                        .expect("paired fragment should report read2 entropy")
                        > filter.config().low_complexity_entropy_bits,
                    "expected read2 to remain above low-complexity threshold, got {metrics:?}"
                );
            }
            other => panic!("expected pass for diverse fragment, got {other:?}"),
        }
    }

    #[test]
    fn all_n_fragment_is_filtered_by_n_content() {
        let filter = QcFilter::default();
        let read1 = "N".repeat(150);
        let read2 = "N".repeat(150);

        let result = filter.filter(&FragmentRecord::paired(read1.as_bytes(), read2.as_bytes()));

        assert_eq!(
            result.reason(),
            Some(&QcFilterReason::NContent {
                observed_fraction: 1.0,
                max_fraction: 0.10,
            })
        );
    }

    #[test]
    fn low_complexity_fragment_is_filtered_distinctly() {
        let filter = QcFilter::default();
        let read1 = "AT".repeat(75);
        let read2 = diverse_sequence(150, 41);

        let result = filter.filter(&FragmentRecord::paired(read1.as_bytes(), read2.as_bytes()));

        match result.reason() {
            Some(QcFilterReason::LowComplexity {
                mate,
                observed_entropy_bits,
                min_entropy_bits,
            }) => {
                assert_eq!(*mate, Mate::Read1);
                assert!(
                    *observed_entropy_bits < *min_entropy_bits,
                    "expected entropy below threshold, got {observed_entropy_bits} >= {min_entropy_bits}"
                );
            }
            other => panic!("expected low-complexity failure, got {other:?}"),
        }
    }

    #[test]
    fn read_length_bounds_are_inclusive() {
        let filter = QcFilter::default();
        let min_r1 = diverse_sequence(100, 43);
        let min_r2 = diverse_sequence(100, 47);
        let max_r1 = diverse_sequence(250, 51);
        let max_r2 = diverse_sequence(250, 57);

        assert!(filter
            .filter(&FragmentRecord::paired(
                min_r1.as_bytes(),
                min_r2.as_bytes()
            ))
            .is_pass());
        assert!(filter
            .filter(&FragmentRecord::paired(
                max_r1.as_bytes(),
                max_r2.as_bytes()
            ))
            .is_pass());
    }

    #[test]
    fn short_fragment_is_filtered_by_length() {
        let filter = QcFilter::default();
        let read1 = diverse_sequence(50, 13);
        let read2 = diverse_sequence(50, 29);

        let result = filter.filter(&FragmentRecord::paired(read1.as_bytes(), read2.as_bytes()));

        assert_eq!(
            result.reason(),
            Some(&QcFilterReason::LengthOutOfBounds {
                mate: Mate::Read1,
                observed: 50,
                min: 100,
                max: 250,
            })
        );
    }

    #[test]
    fn overly_long_fragment_is_filtered_by_length() {
        let filter = QcFilter::default();
        let read1 = diverse_sequence(251, 17);
        let read2 = diverse_sequence(251, 23);

        let result = filter.filter(&FragmentRecord::paired(read1.as_bytes(), read2.as_bytes()));

        assert_eq!(
            result.reason(),
            Some(&QcFilterReason::LengthOutOfBounds {
                mate: Mate::Read1,
                observed: 251,
                min: 100,
                max: 250,
            })
        );
    }

    #[test]
    fn missing_second_mate_is_filtered_as_pair_incomplete() {
        let filter = QcFilter::default();
        let read1 = diverse_sequence(150, 53);

        let result = filter.filter(&FragmentRecord::single(read1.as_bytes()));

        assert_eq!(
            result.reason(),
            Some(&QcFilterReason::PairIncomplete {
                missing_or_empty: Mate::Read2,
            })
        );
    }

    #[test]
    fn empty_second_mate_is_rejected_as_pair_incomplete() {
        let filter = QcFilter::default();
        let read1 = diverse_sequence(150, 59);

        let result = filter.filter(&FragmentRecord::paired(read1.as_bytes(), b""));

        assert_eq!(
            result.reason(),
            Some(&QcFilterReason::PairIncomplete {
                missing_or_empty: Mate::Read2,
            })
        );
    }

    #[test]
    fn qc_stats_count_passes_and_each_failure_bucket() {
        let filter = QcFilter::default();
        let mut stats = QcStats::default();

        let pass_r1 = diverse_sequence(150, 61);
        let pass_r2 = diverse_sequence(150, 71);
        assert!(filter
            .filter_and_record(
                &FragmentRecord::paired(pass_r1.as_bytes(), pass_r2.as_bytes()),
                &mut stats,
            )
            .is_pass());

        let all_n = "N".repeat(150);
        filter.filter_and_record(
            &FragmentRecord::paired(all_n.as_bytes(), all_n.as_bytes()),
            &mut stats,
        );

        let low_complexity = "AT".repeat(75);
        let diverse = diverse_sequence(150, 83);
        filter.filter_and_record(
            &FragmentRecord::paired(low_complexity.as_bytes(), diverse.as_bytes()),
            &mut stats,
        );

        let short = diverse_sequence(50, 97);
        filter.filter_and_record(
            &FragmentRecord::paired(short.as_bytes(), short.as_bytes()),
            &mut stats,
        );

        let missing_mate = diverse_sequence(150, 109);
        filter.filter_and_record(&FragmentRecord::single(missing_mate.as_bytes()), &mut stats);

        assert_eq!(
            stats,
            QcStats {
                total_fragments: 5,
                passed_fragments: 1,
                low_complexity_filtered: 1,
                n_content_filtered: 1,
                length_filtered: 1,
                pair_incomplete: 1,
            }
        );
    }

    #[test]
    fn qc_config_validation_rejects_invalid_thresholds() {
        let invalid_n_fraction = QcConfig {
            max_n_fraction: 1.1,
            ..QcConfig::default()
        };
        assert!(QcFilter::new(invalid_n_fraction).is_err());

        let invalid_lengths = QcConfig {
            min_read_length: 251,
            max_read_length: 250,
            ..QcConfig::default()
        };
        assert!(QcFilter::new(invalid_lengths).is_err());
    }

    fn diverse_sequence(length: usize, seed: u64) -> String {
        let mut state = seed;
        let mut sequence = String::with_capacity(length);

        for _ in 0..length {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);

            let base = match ((state >> 62) & 0b11) as u8 {
                0 => 'A',
                1 => 'C',
                2 => 'G',
                _ => 'T',
            };
            sequence.push(base);
        }

        sequence
    }
}
