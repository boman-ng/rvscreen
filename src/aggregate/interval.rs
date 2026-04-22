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
}
