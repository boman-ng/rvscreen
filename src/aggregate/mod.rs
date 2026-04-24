pub mod interval;

use crate::align::{AlignHit, FragmentAlignResult};
use crate::decision::ProportionEstimate;
use crate::error::{Result, RvScreenError};
use crate::qc::QcStats;
use crate::types::{CandidateCall, DecisionStatus, EvidenceStrength, FragmentClass, HotPathId};
use interval::{covered_bases, merge_intervals, Interval};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum AggregationLevel {
    #[default]
    Accession,
    VirusGroup,
}

#[derive(Debug, Clone, Default)]
pub struct CandidateAggregator {
    level: AggregationLevel,
    total_sampled_fragments: Option<u64>,
    fragments_seen: u64,
    candidates: BTreeMap<HotPathId, CandidateAccumulator>,
}

#[derive(Debug, Clone)]
struct CandidateAccumulator {
    candidate_id: HotPathId,
    level: AggregationLevel,
    accession_or_group: String,
    representative_virus_name: String,
    mixed_virus_names: bool,
    taxid: u64,
    accepted_fragments: u64,
    ambiguous_fragments: u64,
    contig_lengths: BTreeMap<HotPathId, u64>,
    contig_labels: BTreeMap<HotPathId, String>,
    intervals_by_contig: BTreeMap<HotPathId, Vec<Interval>>,
}

impl CandidateAggregator {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn with_level(level: AggregationLevel) -> Self {
        Self {
            level,
            ..Self::default()
        }
    }

    pub fn level(&self) -> AggregationLevel {
        self.level
    }

    pub fn with_total_sampled_fragments(mut self, total_sampled_fragments: u64) -> Self {
        self.total_sampled_fragments = Some(total_sampled_fragments);
        self
    }

    pub fn set_total_sampled_fragments(&mut self, total_sampled_fragments: u64) {
        self.total_sampled_fragments = Some(total_sampled_fragments);
    }

    pub(crate) fn merge_from(&mut self, other: Self) -> Result<()> {
        if self.level != other.level {
            return Err(RvScreenError::validation(
                "aggregate.level",
                format!(
                    "cannot merge aggregation levels `{}` and `{}`",
                    self.level.label(),
                    other.level.label()
                ),
            ));
        }

        self.fragments_seen = self.fragments_seen.saturating_add(other.fragments_seen);
        self.total_sampled_fragments = None;

        for (candidate_key, candidate) in other.candidates {
            match self.candidates.entry(candidate_key) {
                std::collections::btree_map::Entry::Vacant(entry) => {
                    entry.insert(candidate);
                }
                std::collections::btree_map::Entry::Occupied(mut entry) => {
                    entry.get_mut().merge_from(candidate)?;
                }
            }
        }

        Ok(())
    }

    pub fn add_fragment(
        &mut self,
        class: FragmentClass,
        align_result: &FragmentAlignResult,
    ) -> Result<()> {
        self.fragments_seen = self.fragments_seen.saturating_add(1);

        match class {
            FragmentClass::Virus => {
                let hit = align_result.best_virus_hit.as_ref().ok_or_else(|| {
                    RvScreenError::validation(
                        "aggregate.fragment_class",
                        "Virus fragment is missing best_virus_hit alignment evidence",
                    )
                })?;
                self.record_hit(hit, false)
            }
            FragmentClass::Ambiguous => {
                let mut seen_candidates = BTreeSet::new();
                for hit in ambiguous_candidate_hits(align_result) {
                    let candidate_id = candidate_id(self.level, hit)?;
                    if seen_candidates.insert(candidate_id) {
                        self.record_hit(hit, true)?;
                    }
                }
                Ok(())
            }
            FragmentClass::Host | FragmentClass::Unmapped => Ok(()),
        }
    }

    pub fn finalize(&self, _qc_stats: &QcStats) -> Vec<CandidateCall> {
        let total_sampled_fragments = self.total_sampled_fragments.unwrap_or(self.fragments_seen);

        let mut calls = self
            .candidates
            .values()
            .map(|candidate| candidate.to_call(total_sampled_fragments))
            .collect::<Vec<_>>();

        calls.sort_by(|left, right| {
            right
                .accepted_fragments
                .cmp(&left.accepted_fragments)
                .then_with(|| right.nonoverlap_fragments.cmp(&left.nonoverlap_fragments))
                .then_with(|| right.breadth.total_cmp(&left.breadth))
                .then_with(|| left.accession_or_group.cmp(&right.accession_or_group))
        });
        calls
    }

    fn record_hit(&mut self, hit: &AlignHit, ambiguous_only: bool) -> Result<()> {
        let candidate_id = candidate_id(self.level, hit)?;
        let accumulator = self
            .candidates
            .entry(candidate_id)
            .or_insert_with(|| CandidateAccumulator::new(self.level, hit));
        accumulator.observe_hit(hit, ambiguous_only)
    }
}

impl CandidateAccumulator {
    fn new(level: AggregationLevel, hit: &AlignHit) -> Self {
        let target = hit
            .virus_target
            .as_ref()
            .expect("validated virus-target metadata before accumulator initialization");

        Self {
            candidate_id: candidate_id_from_target(level, target),
            level,
            accession_or_group: candidate_label(level, target),
            representative_virus_name: target.virus_name.clone(),
            mixed_virus_names: false,
            taxid: target.taxid,
            accepted_fragments: 0,
            ambiguous_fragments: 0,
            contig_lengths: BTreeMap::new(),
            contig_labels: BTreeMap::new(),
            intervals_by_contig: BTreeMap::new(),
        }
    }

    fn observe_hit(&mut self, hit: &AlignHit, ambiguous_only: bool) -> Result<()> {
        let target = hit.virus_target.as_ref().ok_or_else(|| {
            RvScreenError::validation(
                "aggregate.virus_target",
                format!(
                    "candidate aggregation requires manifest-linked virus metadata for contig `{}`",
                    hit.contig
                ),
            )
        })?;

        let expected_candidate_id = candidate_id_from_target(self.level, target);
        let expected_label = candidate_label(self.level, target);
        if self.candidate_id != expected_candidate_id {
            return Err(RvScreenError::validation(
                "aggregate.candidate_key",
                format!(
                    "candidate key drift for contig `{}`: expected `{}` ({expected_candidate_id}), observed `{}` ({})",
                    hit.contig, expected_label, self.accession_or_group, self.candidate_id
                ),
            ));
        }
        debug_assert_eq!(self.accession_or_group, expected_label);

        if self.taxid != target.taxid {
            return Err(RvScreenError::validation(
                "aggregate.taxid",
                format!(
                    "conflicting taxid values within candidate `{}`: {} vs {}",
                    self.accession_or_group, self.taxid, target.taxid
                ),
            ));
        }

        if self.level == AggregationLevel::Accession
            && self.representative_virus_name != target.virus_name
        {
            return Err(RvScreenError::validation(
                "aggregate.virus_name",
                format!(
                    "accession candidate `{}` saw conflicting virus_name values: `{}` vs `{}`",
                    self.accession_or_group, self.representative_virus_name, target.virus_name
                ),
            ));
        }
        if self.level == AggregationLevel::VirusGroup
            && self.representative_virus_name != target.virus_name
        {
            self.mixed_virus_names = true;
        }

        match self.contig_lengths.get(&hit.contig_id) {
            Some(existing) if *existing != target.genome_length => {
                return Err(RvScreenError::validation(
                    "aggregate.genome_length",
                    format!(
                        "contig `{}` saw conflicting genome_length values: {} vs {}",
                        hit.contig, existing, target.genome_length
                    ),
                ));
            }
            Some(_) => {}
            None => {
                self.contig_lengths
                    .insert(hit.contig_id, target.genome_length);
                self.contig_labels.insert(hit.contig_id, hit.contig.clone());
            }
        }

        if ambiguous_only {
            self.ambiguous_fragments = self.ambiguous_fragments.saturating_add(1);
            return Ok(());
        }

        if hit.reference_end < hit.reference_start {
            return Err(RvScreenError::validation(
                "aggregate.interval",
                format!(
                    "contig `{}` emitted inverted interval [{}, {})",
                    hit.contig, hit.reference_start, hit.reference_end
                ),
            ));
        }

        self.accepted_fragments = self.accepted_fragments.saturating_add(1);
        self.intervals_by_contig
            .entry(hit.contig_id)
            .or_default()
            .push(Interval::new(hit.reference_start, hit.reference_end));
        Ok(())
    }

    fn merge_from(&mut self, other: Self) -> Result<()> {
        if self.level != other.level {
            return Err(RvScreenError::validation(
                "aggregate.level",
                format!(
                    "cannot merge candidate `{}` across aggregation levels `{}` and `{}`",
                    self.accession_or_group,
                    self.level.label(),
                    other.level.label()
                ),
            ));
        }

        if self.candidate_id != other.candidate_id {
            return Err(RvScreenError::validation(
                "aggregate.candidate_key",
                format!(
                    "cannot merge candidate accumulators `{}` ({}) and `{}` ({})",
                    self.accession_or_group,
                    self.candidate_id,
                    other.accession_or_group,
                    other.candidate_id
                ),
            ));
        }

        if self.taxid != other.taxid {
            return Err(RvScreenError::validation(
                "aggregate.taxid",
                format!(
                    "conflicting taxid values within candidate `{}`: {} vs {}",
                    self.accession_or_group, self.taxid, other.taxid
                ),
            ));
        }

        if self.level == AggregationLevel::Accession
            && self.representative_virus_name != other.representative_virus_name
        {
            return Err(RvScreenError::validation(
                "aggregate.virus_name",
                format!(
                    "accession candidate `{}` saw conflicting virus_name values: `{}` vs `{}`",
                    self.accession_or_group,
                    self.representative_virus_name,
                    other.representative_virus_name
                ),
            ));
        }

        if self.level == AggregationLevel::VirusGroup
            && self.representative_virus_name != other.representative_virus_name
        {
            self.mixed_virus_names = true;
        }
        self.mixed_virus_names |= other.mixed_virus_names;

        self.accepted_fragments = self
            .accepted_fragments
            .saturating_add(other.accepted_fragments);
        self.ambiguous_fragments = self
            .ambiguous_fragments
            .saturating_add(other.ambiguous_fragments);

        for (contig, genome_length) in other.contig_lengths {
            match self.contig_lengths.get(&contig) {
                Some(existing) if *existing != genome_length => {
                    let contig_label = other
                        .contig_labels
                        .get(&contig)
                        .or_else(|| self.contig_labels.get(&contig))
                        .cloned()
                        .unwrap_or_else(|| format!("contig#{contig}"));
                    return Err(RvScreenError::validation(
                        "aggregate.genome_length",
                        format!(
                            "contig `{}` saw conflicting genome_length values: {} vs {}",
                            contig_label, existing, genome_length
                        ),
                    ));
                }
                Some(_) => {}
                None => {
                    self.contig_lengths.insert(contig, genome_length);
                }
            }
        }

        for (contig, label) in other.contig_labels {
            self.contig_labels.entry(contig).or_insert(label);
        }

        for (contig, mut intervals) in other.intervals_by_contig {
            self.intervals_by_contig
                .entry(contig)
                .or_default()
                .append(&mut intervals);
        }

        Ok(())
    }

    fn to_call(&self, total_sampled_fragments: u64) -> CandidateCall {
        let mut nonoverlap_fragments = 0u64;
        let mut covered_bases_total = 0u64;
        for (contig, intervals) in &self.intervals_by_contig {
            let merged = merge_intervals(intervals.iter().copied());
            nonoverlap_fragments = nonoverlap_fragments.saturating_add(merged.len() as u64);
            covered_bases_total = covered_bases_total.saturating_add(covered_bases(&merged));
            debug_assert!(self.contig_lengths.contains_key(contig));
        }

        let genome_length = self.contig_lengths.values().copied().sum::<u64>();
        let breadth = ratio(covered_bases_total, genome_length);
        let stats =
            ProportionEstimate::from_counts(self.accepted_fragments, total_sampled_fragments);

        CandidateCall {
            virus_name: self.output_virus_name(),
            taxid: self.taxid,
            accession_or_group: self.accession_or_group.clone(),
            accepted_fragments: self.accepted_fragments,
            nonoverlap_fragments,
            raw_fraction: stats.raw_fraction,
            unique_fraction: stats.unique_fraction,
            fraction_ci_95: stats.fraction_ci_95(),
            clopper_pearson_upper: stats.clopper_pearson_upper,
            breadth,
            ambiguous_fragments: self.ambiguous_fragments,
            background_ratio: 0.0,
            decision: DecisionStatus::Indeterminate,
            decision_reasons: vec![
                format!("aggregation_level={}", self.level.label()),
                format!("fraction_ci_95_label={}", stats.ci_label),
                "background_ratio_pending_task_19".to_string(),
                "decision_pending_task_16".to_string(),
            ],
            evidence_strength: evidence_strength(nonoverlap_fragments, breadth),
        }
    }

    fn output_virus_name(&self) -> String {
        if self.level == AggregationLevel::VirusGroup && self.mixed_virus_names {
            self.accession_or_group.clone()
        } else {
            self.representative_virus_name.clone()
        }
    }
}

impl AggregationLevel {
    fn label(self) -> &'static str {
        match self {
            Self::Accession => "accession",
            Self::VirusGroup => "virus_group",
        }
    }
}

fn ambiguous_candidate_hits(align_result: &FragmentAlignResult) -> Vec<&AlignHit> {
    align_result
        .best_virus_hit
        .iter()
        .chain(
            align_result
                .secondary_hits
                .iter()
                .filter(|hit| hit.virus_target.is_some()),
        )
        .filter(|hit| hit.virus_target.is_some())
        .collect()
}

fn candidate_id(level: AggregationLevel, hit: &AlignHit) -> Result<HotPathId> {
    let target = hit.virus_target.as_ref().ok_or_else(|| {
        RvScreenError::validation(
            "aggregate.virus_target",
            format!(
                "candidate aggregation requires manifest-linked virus metadata for contig `{}`",
                hit.contig
            ),
        )
    })?;
    let label = candidate_label(level, target);
    if label.trim().is_empty() {
        return Err(RvScreenError::validation(
            "aggregate.candidate_key",
            format!(
                "candidate aggregation produced an empty {} key for contig `{}`",
                level.label(),
                hit.contig
            ),
        ));
    }
    Ok(candidate_id_from_target(level, target))
}

fn candidate_id_from_target(
    level: AggregationLevel,
    target: &crate::align::VirusTarget,
) -> HotPathId {
    match level {
        AggregationLevel::Accession => target.target_id,
        AggregationLevel::VirusGroup => target.group_id,
    }
}

fn candidate_label(level: AggregationLevel, target: &crate::align::VirusTarget) -> String {
    match level {
        AggregationLevel::Accession => target.accession.clone(),
        AggregationLevel::VirusGroup => target.group.clone(),
    }
}

fn evidence_strength(nonoverlap_fragments: u64, breadth: f64) -> EvidenceStrength {
    if nonoverlap_fragments >= 3 && breadth > 0.0 {
        EvidenceStrength::High
    } else if nonoverlap_fragments >= 2 || breadth > 0.0 {
        EvidenceStrength::Medium
    } else {
        EvidenceStrength::Low
    }
}

fn ratio(numerator: u64, denominator: u64) -> f64 {
    if denominator == 0 {
        0.0
    } else {
        numerator as f64 / denominator as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::VirusTarget;

    #[test]
    fn test_nonoverlap_fragments_merge_overlapping_intervals() {
        let mut aggregator = CandidateAggregator::new().with_total_sampled_fragments(5);

        for (start, end) in [(0, 100), (50, 150), (200, 300), (250, 350), (500, 600)] {
            aggregator
                .add_fragment(
                    FragmentClass::Virus,
                    &virus_alignment("virus-a", "ACC-A", "ebv", 10_376, 1000, start, end),
                )
                .expect("virus fragment should aggregate");
        }

        let calls = aggregator.finalize(&qc_stats(5));
        assert_eq!(calls.len(), 1);
        assert_eq!(calls[0].accepted_fragments, 5);
        assert_eq!(calls[0].nonoverlap_fragments, 3);
        assert!((calls[0].breadth - 0.4).abs() < 1e-9);
        assert_eq!(calls[0].evidence_strength, EvidenceStrength::High);
    }

    #[test]
    fn test_breadth_uses_covered_bases_over_genome_length() {
        let mut aggregator = CandidateAggregator::new().with_total_sampled_fragments(4);

        aggregator
            .add_fragment(
                FragmentClass::Virus,
                &virus_alignment("virus-b", "ACC-B", "influenza-a", 11_320, 1000, 0, 100),
            )
            .expect("first breadth fragment should aggregate");
        aggregator
            .add_fragment(
                FragmentClass::Virus,
                &virus_alignment("virus-b", "ACC-B", "influenza-a", 11_320, 1000, 200, 300),
            )
            .expect("second breadth fragment should aggregate");

        let calls = aggregator.finalize(&qc_stats(4));
        assert_eq!(calls.len(), 1);
        assert_eq!(calls[0].nonoverlap_fragments, 2);
        assert!((calls[0].breadth - 0.2).abs() < 1e-9);
        assert!((calls[0].raw_fraction - 0.5).abs() < 1e-9);
        assert!((calls[0].unique_fraction - 0.5).abs() < 1e-9);
        assert!((calls[0].fraction_ci_95[0] - 0.150_035_708_820_171_48).abs() < 1e-12);
        assert!((calls[0].fraction_ci_95[1] - 0.849_964_291_179_828_5).abs() < 1e-12);
        assert!(calls[0]
            .decision_reasons
            .iter()
            .any(|reason| reason.contains("sampling-only CI")));
    }

    #[test]
    fn test_group_aggregation_merges_related_accessions() {
        let mut aggregator = CandidateAggregator::with_level(AggregationLevel::VirusGroup)
            .with_total_sampled_fragments(6);

        aggregator
            .add_fragment(
                FragmentClass::Virus,
                &virus_alignment("virus-ebv-a", "A1", "ebv", 10_376, 400, 0, 90),
            )
            .expect("first group fragment should aggregate");
        aggregator
            .add_fragment(
                FragmentClass::Virus,
                &virus_alignment("virus-ebv-b", "A2", "ebv", 10_376, 400, 100, 180),
            )
            .expect("second group fragment should aggregate");
        aggregator
            .add_fragment(
                FragmentClass::Virus,
                &virus_alignment("virus-ebv-a", "A1", "ebv", 10_376, 400, 250, 320),
            )
            .expect("third group fragment should aggregate");

        let calls = aggregator.finalize(&qc_stats(6));
        assert_eq!(calls.len(), 1);
        assert_eq!(calls[0].accession_or_group, "ebv");
        assert_eq!(calls[0].accepted_fragments, 3);
        assert_eq!(calls[0].taxid, 10_376);
    }

    #[test]
    fn test_ambiguous_fragments_count_once_per_candidate() {
        let mut aggregator = CandidateAggregator::new().with_total_sampled_fragments(3);
        aggregator
            .add_fragment(
                FragmentClass::Ambiguous,
                &ambiguous_alignment_same_candidate("virus-c", "ACC-C", "hsv", 10_294, 700),
            )
            .expect("ambiguous fragment should aggregate once");

        let calls = aggregator.finalize(&qc_stats(3));
        assert_eq!(calls.len(), 1);
        assert_eq!(calls[0].accepted_fragments, 0);
        assert_eq!(calls[0].ambiguous_fragments, 1);
        assert_eq!(calls[0].nonoverlap_fragments, 0);
        assert_eq!(calls[0].evidence_strength, EvidenceStrength::Low);
    }

    #[test]
    fn test_merge_from_combines_round_deltas_into_cumulative_candidate_state() {
        let mut early_round = CandidateAggregator::new();
        early_round
            .add_fragment(
                FragmentClass::Virus,
                &virus_alignment("virus-d", "ACC-D", "adeno", 10_580, 500, 0, 90),
            )
            .expect("early representative delta should aggregate");

        let mut later_round = CandidateAggregator::new();
        later_round
            .add_fragment(
                FragmentClass::Virus,
                &virus_alignment("virus-d", "ACC-D", "adeno", 10_580, 500, 200, 290),
            )
            .expect("later representative delta should aggregate");
        later_round
            .add_fragment(
                FragmentClass::Ambiguous,
                &ambiguous_alignment_same_candidate("virus-d", "ACC-D", "adeno", 10_580, 500),
            )
            .expect("later ambiguous representative delta should aggregate");

        early_round
            .merge_from(later_round)
            .expect("round deltas should merge cleanly");
        early_round.set_total_sampled_fragments(4);

        let calls = early_round.finalize(&qc_stats(4));
        assert_eq!(calls.len(), 1);
        assert_eq!(calls[0].accepted_fragments, 2);
        assert_eq!(calls[0].ambiguous_fragments, 1);
        assert_eq!(calls[0].nonoverlap_fragments, 2);
        assert!((calls[0].breadth - 0.36).abs() < 1e-9);
        assert!((calls[0].raw_fraction - 0.5).abs() < 1e-9);
    }

    #[test]
    fn denominator_uses_sampled_fragments_instead_of_qc_pass_volume() {
        let mut aggregator = CandidateAggregator::new().with_total_sampled_fragments(5);

        aggregator
            .add_fragment(
                FragmentClass::Virus,
                &virus_alignment("virus-e", "ACC-E", "noro", 12_753, 500, 0, 90),
            )
            .expect("first sampled fragment should aggregate");
        aggregator
            .add_fragment(
                FragmentClass::Virus,
                &virus_alignment("virus-e", "ACC-E", "noro", 12_753, 500, 120, 210),
            )
            .expect("second sampled fragment should aggregate");

        let calls = aggregator.finalize(&qc_stats(10));
        assert_eq!(calls.len(), 1);
        assert!((calls[0].raw_fraction - 0.4).abs() < 1e-12);
        assert!((calls[0].unique_fraction - 0.4).abs() < 1e-12);
        assert_eq!(
            calls[0].fraction_ci_95,
            crate::decision::stats::wilson_ci_95(2, 5)
        );
        assert_eq!(
            calls[0].clopper_pearson_upper,
            crate::decision::stats::clopper_pearson_upper_95_one_sided(2, 5)
        );
    }

    fn virus_alignment(
        contig: &str,
        accession: &str,
        group: &str,
        taxid: u64,
        genome_length: u64,
        reference_start: u64,
        reference_end: u64,
    ) -> FragmentAlignResult {
        FragmentAlignResult {
            best_host_hit: None,
            best_virus_hit: Some(virus_hit(
                contig,
                accession,
                group,
                taxid,
                genome_length,
                reference_start,
                reference_end,
            )),
            secondary_hits: Vec::new(),
        }
    }

    fn ambiguous_alignment_same_candidate(
        contig: &str,
        accession: &str,
        group: &str,
        taxid: u64,
        genome_length: u64,
    ) -> FragmentAlignResult {
        let best = virus_hit(contig, accession, group, taxid, genome_length, 50, 120);
        let secondary = virus_hit(contig, accession, group, taxid, genome_length, 75, 140);

        FragmentAlignResult {
            best_host_hit: None,
            best_virus_hit: Some(best),
            secondary_hits: vec![AlignHit {
                is_primary: false,
                ..secondary
            }],
        }
    }

    fn virus_hit(
        contig: &str,
        accession: &str,
        group: &str,
        taxid: u64,
        genome_length: u64,
        reference_start: u64,
        reference_end: u64,
    ) -> AlignHit {
        AlignHit {
            contig_id: 1,
            contig: contig.to_string(),
            mapq: 60,
            as_score: Some(100),
            nm: Some(1),
            is_primary: true,
            is_supplementary: false,
            r1_mapped: true,
            r2_mapped: true,
            pair_consistent: true,
            reference_start,
            reference_end,
            virus_target: Some(VirusTarget {
                target_id: 11,
                contig_id: 1,
                group_id: 21,
                accession: accession.to_string(),
                group: group.to_string(),
                taxid,
                virus_name: format!("{group}-virus"),
                genome_length,
            }),
        }
    }

    fn qc_stats(passed_fragments: u64) -> QcStats {
        QcStats {
            total_fragments: passed_fragments,
            passed_fragments,
            ..QcStats::default()
        }
    }
}
