#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Interval {
    pub start: u64,
    pub end: u64,
}

impl Interval {
    pub fn new(start: u64, end: u64) -> Self {
        Self { start, end }
    }

    pub fn len(self) -> u64 {
        self.end.saturating_sub(self.start)
    }

    pub fn is_empty(self) -> bool {
        self.len() == 0
    }
}

/// Algorithm label emitted in v3 reports for corrected non-overlap support.
pub const NONOVERLAP_FRAGMENTS_ALGORITHM: &str = "max_nonoverlapping_intervals_v1";

/// Count the maximum-cardinality set of non-overlapping half-open intervals.
///
/// Empty intervals are ignored. Intervals are sorted by `(end asc, start asc)`
/// and greedily accepted when `start >= last_accepted_end`.
pub fn max_nonoverlapping_intervals(intervals: impl IntoIterator<Item = Interval>) -> u64 {
    let mut sorted = intervals
        .into_iter()
        .filter(|interval| !interval.is_empty())
        .collect::<Vec<_>>();
    sorted.sort_by(|left, right| left.end.cmp(&right.end).then(left.start.cmp(&right.start)));

    let mut accepted = 0u64;
    let mut last_accepted_end: Option<u64> = None;
    for interval in sorted {
        if last_accepted_end.is_none_or(|end| interval.start >= end) {
            accepted = accepted.saturating_add(1);
            last_accepted_end = Some(interval.end);
        }
    }
    accepted
}

pub fn merge_intervals(intervals: impl IntoIterator<Item = Interval>) -> Vec<Interval> {
    let mut sorted = intervals.into_iter().collect::<Vec<_>>();
    sorted.sort_by(|left, right| left.start.cmp(&right.start).then(left.end.cmp(&right.end)));

    let mut merged: Vec<Interval> = Vec::with_capacity(sorted.len());
    for interval in sorted {
        match merged.last_mut() {
            Some(current) if interval.start < current.end => {
                current.end = current.end.max(interval.end);
            }
            _ => merged.push(interval),
        }
    }

    merged
}

pub fn covered_bases(intervals: &[Interval]) -> u64 {
    intervals.iter().copied().map(Interval::len).sum()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn merge_intervals_coalesces_only_overlaps() {
        let merged = merge_intervals([
            Interval::new(0, 100),
            Interval::new(50, 150),
            Interval::new(150, 200),
            Interval::new(180, 240),
        ]);

        assert_eq!(
            merged,
            vec![Interval::new(0, 150), Interval::new(150, 240),]
        );
        assert_eq!(covered_bases(&merged), 240);
    }

    #[test]
    fn max_nonoverlapping_intervals_excludes_empty_and_prefers_earliest_end() {
        let count = max_nonoverlapping_intervals([
            Interval::new(0, 10),
            Interval::new(1, 2),
            Interval::new(2, 3),
            Interval::new(3, 3),
            Interval::new(3, 4),
            Interval::new(3, 9),
            Interval::new(9, 10),
        ]);
        assert_eq!(count, 4);
    }
}
