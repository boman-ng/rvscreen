use statrs::distribution::{Beta, ContinuousCDF};

pub const SAMPLING_ONLY_CI_LABEL: &str =
    "sampling-only CI: accepted-fragment sampling uncertainty only; not diagnostic confidence";

const WILSON_Z_95: f64 = 1.96;
const ONE_SIDED_CONFIDENCE_95: f64 = 0.95;

#[derive(Debug, Clone, PartialEq)]
pub struct ProportionEstimate {
    pub raw_fraction: f64,
    pub unique_fraction: f64,
    pub wilson_ci_lower: f64,
    pub wilson_ci_upper: f64,
    pub clopper_pearson_upper: f64,
    pub ci_label: &'static str,
}

impl ProportionEstimate {
    pub fn from_counts(
        accepted_fragments: u64,
        total_sampled_fragments: u64,
        qc_passing_fragments: u64,
    ) -> Self {
        let raw_fraction = ratio(accepted_fragments, total_sampled_fragments);
        let unique_fraction = ratio(accepted_fragments, qc_passing_fragments);
        let [wilson_ci_lower, wilson_ci_upper] =
            wilson_ci_95(accepted_fragments, qc_passing_fragments);
        let clopper_pearson_upper =
            clopper_pearson_upper_95_one_sided(accepted_fragments, qc_passing_fragments);

        Self {
            raw_fraction,
            unique_fraction,
            wilson_ci_lower,
            wilson_ci_upper,
            clopper_pearson_upper,
            ci_label: SAMPLING_ONLY_CI_LABEL,
        }
    }

    pub fn fraction_ci_95(&self) -> [f64; 2] {
        [self.wilson_ci_lower, self.wilson_ci_upper]
    }
}

pub fn wilson_ci_95(successes: u64, trials: u64) -> [f64; 2] {
    if trials == 0 {
        return [0.0, 0.0];
    }

    let successes = successes.min(trials);
    let n = trials as f64;
    let p = successes as f64 / n;
    let z2 = WILSON_Z_95 * WILSON_Z_95;
    let denominator = 1.0 + z2 / n;
    let center = (p + z2 / (2.0 * n)) / denominator;
    let half_width =
        (WILSON_Z_95 * ((p * (1.0 - p) / n) + z2 / (4.0 * n * n)).sqrt()) / denominator;

    [
        clamp_probability(center - half_width),
        clamp_probability(center + half_width),
    ]
}

pub fn clopper_pearson_upper_95_one_sided(successes: u64, trials: u64) -> f64 {
    if trials == 0 {
        return 0.0;
    }

    let successes = successes.min(trials);
    if successes == 0 {
        return clamp_probability(1.0 - (1.0 - ONE_SIDED_CONFIDENCE_95).powf(1.0 / trials as f64));
    }
    if successes == trials {
        return 1.0;
    }

    let distribution = Beta::new((successes + 1) as f64, (trials - successes) as f64)
        .expect("valid beta parameters for one-sided Clopper-Pearson upper bound");

    let mut lower = 0.0;
    let mut upper = 1.0;
    for _ in 0..128 {
        let midpoint = (lower + upper) / 2.0;
        if distribution.cdf(midpoint) < ONE_SIDED_CONFIDENCE_95 {
            lower = midpoint;
        } else {
            upper = midpoint;
        }
    }

    clamp_probability((lower + upper) / 2.0)
}

fn ratio(numerator: u64, denominator: u64) -> f64 {
    if denominator == 0 {
        0.0
    } else {
        numerator as f64 / denominator as f64
    }
}

fn clamp_probability(value: f64) -> f64 {
    value.clamp(0.0, 1.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn wilson_ci_matches_known_value() {
        let estimate = ProportionEstimate::from_counts(50, 10_000, 10_000);

        assert!((estimate.raw_fraction - 0.005).abs() < 1e-12);
        assert!((estimate.unique_fraction - 0.005).abs() < 1e-12);
        assert!((estimate.wilson_ci_lower - 0.003_794_881_950_804_586).abs() < 1e-12);
        assert!((estimate.wilson_ci_upper - 0.006_585_290_402_184_288).abs() < 1e-12);
        assert!(estimate.wilson_ci_lower < estimate.unique_fraction);
        assert!(estimate.unique_fraction < estimate.wilson_ci_upper);
    }

    #[test]
    fn zero_successes_keep_unique_fraction_and_lower_bound_zero() {
        let estimate = ProportionEstimate::from_counts(0, 100_000, 100_000);

        assert_eq!(estimate.raw_fraction, 0.0);
        assert_eq!(estimate.unique_fraction, 0.0);
        assert_eq!(estimate.wilson_ci_lower, 0.0);
        assert!(estimate.wilson_ci_upper > 0.0);
        assert!(estimate.clopper_pearson_upper > 0.0);
        assert_eq!(estimate.ci_label, SAMPLING_ONLY_CI_LABEL);
    }

    #[test]
    fn clopper_pearson_upper_for_zero_successes_matches_reference() {
        let estimate = ProportionEstimate::from_counts(0, 400_000, 400_000);

        assert!((estimate.clopper_pearson_upper - 7.489_302_638_941_098e-6).abs() < 1e-10);
    }

    #[test]
    fn clopper_pearson_upper_for_very_low_count_matches_reference() {
        let estimate = ProportionEstimate::from_counts(1, 10_000, 10_000);

        assert!((estimate.clopper_pearson_upper - 4.742_976_591_654_568_5e-4).abs() < 1e-10);
    }
}
