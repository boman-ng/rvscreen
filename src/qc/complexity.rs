pub fn shannon_entropy_2mer(sequence: &[u8]) -> f64 {
    let mut counts = [0usize; 16];
    let mut total = 0usize;

    for window in sequence.windows(2) {
        if let (Some(left), Some(right)) = (base_index(window[0]), base_index(window[1])) {
            counts[left * 4 + right] += 1;
            total += 1;
        }
    }

    if total == 0 {
        return 0.0;
    }

    counts
        .into_iter()
        .filter(|count| *count > 0)
        .map(|count| {
            let probability = count as f64 / total as f64;
            -probability * probability.log2()
        })
        .sum()
}

fn base_index(base: u8) -> Option<usize> {
    match base.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::shannon_entropy_2mer;

    #[test]
    fn homopolymer_entropy_is_zero() {
        assert_eq!(shannon_entropy_2mer(b"AAAAAAAAAAAA"), 0.0);
    }

    #[test]
    fn alternating_repeat_stays_below_one_bit() {
        assert!(shannon_entropy_2mer(&b"AT".repeat(75)) < 1.0);
    }

    #[test]
    fn diverse_sequence_exceeds_low_complexity_threshold() {
        let sequence = b"ACGTTGCATGTCGCATGATGCATGAGAGCTACGTAGCTAGTCGATCGTACGATCGTAGCTAGCAT";
        assert!(shannon_entropy_2mer(sequence) > 1.0);
    }
}
