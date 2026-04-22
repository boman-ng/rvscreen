use crate::align::{AlignHit, FragmentAlignResult};
use crate::types::{FragmentClass, FragmentRules};

#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct AdjudicationInputs {
    pub low_complexity: bool,
}

#[derive(Debug, Clone, PartialEq)]
pub struct AdjudicationResult {
    pub class: FragmentClass,
    pub virus_contig: Option<String>,
    pub confidence_flags: Vec<String>,
    pub as_diff: Option<i32>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct FragmentAdjudicator {
    rules: FragmentRules,
}

impl FragmentAdjudicator {
    pub fn new(rules: &FragmentRules) -> Self {
        Self {
            rules: rules.clone(),
        }
    }

    pub fn rules(&self) -> &FragmentRules {
        &self.rules
    }

    pub fn classify(
        &self,
        result: &FragmentAlignResult,
        inputs: AdjudicationInputs,
    ) -> FragmentClass {
        self.adjudicate(result, inputs).class
    }

    pub fn adjudicate(
        &self,
        result: &FragmentAlignResult,
        inputs: AdjudicationInputs,
    ) -> AdjudicationResult {
        let host_hit = result.best_host_hit.as_ref();
        let virus_hit = result.best_virus_hit.as_ref();
        let host_primary = host_hit.is_some_and(is_default_counting_hit);
        let virus_primary = virus_hit.is_some_and(is_default_counting_hit);
        let has_secondary_or_supplementary = !result.secondary_hits.is_empty();
        let stronger_host =
            virus_hit.is_some_and(|virus_hit| stronger_host_explains(host_hit, virus_hit));
        let as_diff = calculate_as_diff(virus_hit, host_hit);
        let mapq_ok = virus_hit.is_some_and(|hit| hit.mapq >= self.rules.min_mapq);
        let as_diff_ok = as_diff.is_some_and(|diff| diff >= self.rules.min_as_diff);
        let nm_ok = virus_hit.is_some_and(|hit| nm_within_limit(hit, self.rules.max_nm));
        let pair_ok = virus_hit
            .is_some_and(|hit| !self.rules.require_pair_consistency || hit.pair_consistent);
        let low_complexity_ok = !inputs.low_complexity;
        let host_not_stronger = virus_hit.is_some() && !stronger_host;
        let virus_high_confidence = virus_hit.is_some()
            && virus_primary
            && mapq_ok
            && as_diff_ok
            && nm_ok
            && low_complexity_ok
            && host_not_stronger
            && pair_ok;

        let mut confidence_flags = Vec::new();
        if has_secondary_or_supplementary {
            confidence_flags.push("secondary_or_supplementary_evidence_present".to_string());
        }

        match host_hit {
            Some(hit) if is_default_counting_hit(hit) => {
                confidence_flags.push("host_primary".to_string())
            }
            Some(_) => confidence_flags.push("host_not_primary".to_string()),
            None => confidence_flags.push("no_host_hit".to_string()),
        }

        match virus_hit {
            Some(hit) if is_default_counting_hit(hit) => {
                confidence_flags.push("virus_primary".to_string())
            }
            Some(_) => confidence_flags.push("virus_not_primary".to_string()),
            None => confidence_flags.push("no_virus_hit".to_string()),
        }

        if virus_hit.is_some() {
            confidence_flags.push(if mapq_ok {
                "mapq_ok".to_string()
            } else {
                "mapq_below_min".to_string()
            });
            confidence_flags.push(match as_diff {
                Some(diff) if diff >= self.rules.min_as_diff => "as_diff_ok".to_string(),
                Some(_) => "as_diff_below_min".to_string(),
                None => "as_diff_unavailable".to_string(),
            });
            confidence_flags.push(if nm_ok {
                "nm_ok".to_string()
            } else {
                "nm_above_max_or_unavailable".to_string()
            });
            confidence_flags.push(if low_complexity_ok {
                "not_low_complexity".to_string()
            } else {
                "low_complexity".to_string()
            });
            confidence_flags.push(if stronger_host {
                "host_stronger_or_equal".to_string()
            } else {
                "host_not_stronger".to_string()
            });
            confidence_flags.push(if pair_ok {
                "pair_consistent".to_string()
            } else {
                "pair_inconsistent".to_string()
            });
        }

        let class = if virus_high_confidence {
            FragmentClass::Virus
        } else if host_primary && (!virus_primary || stronger_host) {
            FragmentClass::Host
        } else if virus_hit.is_some() {
            if stronger_host {
                FragmentClass::Host
            } else {
                FragmentClass::Ambiguous
            }
        } else if host_primary {
            FragmentClass::Host
        } else if host_hit.is_some() || has_secondary_or_supplementary {
            FragmentClass::Ambiguous
        } else {
            confidence_flags.push("unmapped".to_string());
            FragmentClass::Unmapped
        };

        AdjudicationResult {
            class,
            virus_contig: virus_hit.map(|hit| hit.contig.clone()),
            confidence_flags,
            as_diff,
        }
    }
}

fn is_default_counting_hit(hit: &AlignHit) -> bool {
    hit.is_primary && !hit.is_supplementary
}

fn stronger_host_explains(host_hit: Option<&AlignHit>, virus_hit: &AlignHit) -> bool {
    let Some(host_hit) = host_hit.filter(|hit| is_default_counting_hit(hit)) else {
        return false;
    };
    let (Some(host_as), Some(virus_as)) = (host_hit.as_score, virus_hit.as_score) else {
        return false;
    };

    host_as >= virus_as
}

fn calculate_as_diff(virus_hit: Option<&AlignHit>, host_hit: Option<&AlignHit>) -> Option<i32> {
    let virus_as = virus_hit?.as_score?;

    match host_hit {
        Some(host_hit) => host_hit.as_score.map(|host_as| virus_as - host_as),
        None => Some(virus_as),
    }
}

fn nm_within_limit(hit: &AlignHit, max_nm: u32) -> bool {
    let Some(nm) = hit.nm else {
        return false;
    };
    let Ok(nm) = u32::try_from(nm) else {
        return false;
    };

    nm <= max_nm
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_high_confidence_virus_all_conditions_met() {
        let adjudicator = FragmentAdjudicator::new(&fragment_rules());
        let result =
            adjudicator.adjudicate(&high_confidence_result(), AdjudicationInputs::default());

        assert_eq!(result.class, FragmentClass::Virus);
        assert_eq!(result.virus_contig.as_deref(), Some("virus-alpha"));
        assert_eq!(result.as_diff, Some(30));
        assert!(result
            .confidence_flags
            .contains(&"virus_primary".to_string()));
        assert!(result.confidence_flags.contains(&"mapq_ok".to_string()));
        assert!(result.confidence_flags.contains(&"as_diff_ok".to_string()));
        assert!(result.confidence_flags.contains(&"nm_ok".to_string()));
        assert!(result
            .confidence_flags
            .contains(&"not_low_complexity".to_string()));
        assert!(result
            .confidence_flags
            .contains(&"host_not_stronger".to_string()));
        assert!(result
            .confidence_flags
            .contains(&"pair_consistent".to_string()));
        assert_eq!(
            adjudicator.classify(&high_confidence_result(), AdjudicationInputs::default()),
            FragmentClass::Virus
        );
    }

    #[test]
    fn test_non_primary_virus_hit_is_ambiguous() {
        let adjudicator = FragmentAdjudicator::new(&fragment_rules());
        let mut alignment = high_confidence_result_without_host();
        alignment.best_virus_hit.as_mut().unwrap().is_primary = false;

        let result = adjudicator.adjudicate(&alignment, AdjudicationInputs::default());

        assert_eq!(result.class, FragmentClass::Ambiguous);
        assert!(result
            .confidence_flags
            .contains(&"virus_not_primary".to_string()));
    }

    #[test]
    fn test_low_mapq_downgrades_to_ambiguous() {
        let adjudicator = FragmentAdjudicator::new(&fragment_rules());
        let mut alignment = high_confidence_result();
        alignment.best_virus_hit.as_mut().unwrap().mapq = 10;

        let result = adjudicator.adjudicate(&alignment, AdjudicationInputs::default());

        assert_eq!(result.class, FragmentClass::Ambiguous);
        assert!(result
            .confidence_flags
            .contains(&"mapq_below_min".to_string()));
    }

    #[test]
    fn test_insufficient_as_diff_downgrades_to_ambiguous() {
        let adjudicator = FragmentAdjudicator::new(&fragment_rules());
        let mut alignment = high_confidence_result();
        alignment.best_host_hit.as_mut().unwrap().as_score = Some(82);

        let result = adjudicator.adjudicate(&alignment, AdjudicationInputs::default());

        assert_eq!(result.class, FragmentClass::Ambiguous);
        assert_eq!(result.as_diff, Some(8));
        assert!(result
            .confidence_flags
            .contains(&"as_diff_below_min".to_string()));
    }

    #[test]
    fn test_high_nm_downgrades_to_ambiguous() {
        let adjudicator = FragmentAdjudicator::new(&fragment_rules());
        let mut alignment = high_confidence_result();
        alignment.best_virus_hit.as_mut().unwrap().nm = Some(9);

        let result = adjudicator.adjudicate(&alignment, AdjudicationInputs::default());

        assert_eq!(result.class, FragmentClass::Ambiguous);
        assert!(result
            .confidence_flags
            .contains(&"nm_above_max_or_unavailable".to_string()));
    }

    #[test]
    fn test_low_complexity_downgrades_to_ambiguous() {
        let adjudicator = FragmentAdjudicator::new(&fragment_rules());

        let result = adjudicator.adjudicate(
            &high_confidence_result(),
            AdjudicationInputs {
                low_complexity: true,
            },
        );

        assert_eq!(result.class, FragmentClass::Ambiguous);
        assert!(result
            .confidence_flags
            .contains(&"low_complexity".to_string()));
    }

    #[test]
    fn test_stronger_host_hit_wins() {
        let adjudicator = FragmentAdjudicator::new(&fragment_rules());
        let mut alignment = high_confidence_result();
        alignment.best_host_hit.as_mut().unwrap().as_score = Some(95);

        let result = adjudicator.adjudicate(&alignment, AdjudicationInputs::default());

        assert_eq!(result.class, FragmentClass::Host);
        assert_eq!(result.as_diff, Some(-5));
        assert!(result
            .confidence_flags
            .contains(&"host_stronger_or_equal".to_string()));
    }

    #[test]
    fn test_pair_inconsistency_downgrades_to_ambiguous() {
        let adjudicator = FragmentAdjudicator::new(&fragment_rules());
        let mut alignment = high_confidence_result();
        alignment.best_virus_hit.as_mut().unwrap().pair_consistent = false;

        let result = adjudicator.adjudicate(&alignment, AdjudicationInputs::default());

        assert_eq!(result.class, FragmentClass::Ambiguous);
        assert!(result
            .confidence_flags
            .contains(&"pair_inconsistent".to_string()));
    }

    #[test]
    fn test_secondary_or_supplementary_only_hits_do_not_count_as_virus() {
        let adjudicator = FragmentAdjudicator::new(&fragment_rules());
        let alignment = FragmentAlignResult {
            best_host_hit: None,
            best_virus_hit: None,
            secondary_hits: vec![AlignHit {
                contig: "virus-secondary".to_string(),
                mapq: 42,
                as_score: Some(88),
                nm: Some(2),
                is_primary: false,
                is_supplementary: false,
                r1_mapped: true,
                r2_mapped: false,
                pair_consistent: false,
                reference_start: 0,
                reference_end: 100,
                virus_target: None,
            }],
        };

        let result = adjudicator.adjudicate(&alignment, AdjudicationInputs::default());

        assert_eq!(result.class, FragmentClass::Ambiguous);
        assert!(result
            .confidence_flags
            .contains(&"secondary_or_supplementary_evidence_present".to_string()));
        assert_eq!(result.virus_contig, None);
    }

    #[test]
    fn test_host_only_primary_hit_classifies_as_host() {
        let adjudicator = FragmentAdjudicator::new(&fragment_rules());
        let alignment = FragmentAlignResult {
            best_host_hit: Some(host_hit(70)),
            best_virus_hit: None,
            secondary_hits: Vec::new(),
        };

        let result = adjudicator.adjudicate(&alignment, AdjudicationInputs::default());

        assert_eq!(result.class, FragmentClass::Host);
        assert!(result
            .confidence_flags
            .contains(&"no_virus_hit".to_string()));
    }

    #[test]
    fn test_without_host_or_virus_hits_is_unmapped() {
        let adjudicator = FragmentAdjudicator::new(&fragment_rules());
        let result = adjudicator.adjudicate(
            &FragmentAlignResult::default(),
            AdjudicationInputs::default(),
        );

        assert_eq!(result.class, FragmentClass::Unmapped);
        assert!(result.confidence_flags.contains(&"unmapped".to_string()));
    }

    fn fragment_rules() -> FragmentRules {
        FragmentRules {
            min_mapq: 20,
            min_as_diff: 12,
            max_nm: 8,
            require_pair_consistency: true,
        }
    }

    fn high_confidence_result() -> FragmentAlignResult {
        FragmentAlignResult {
            best_host_hit: Some(host_hit(60)),
            best_virus_hit: Some(virus_hit(90)),
            secondary_hits: Vec::new(),
        }
    }

    fn high_confidence_result_without_host() -> FragmentAlignResult {
        FragmentAlignResult {
            best_host_hit: None,
            best_virus_hit: Some(virus_hit(90)),
            secondary_hits: Vec::new(),
        }
    }

    fn host_hit(as_score: i32) -> AlignHit {
        AlignHit {
            contig: "host-chr1".to_string(),
            mapq: 45,
            as_score: Some(as_score),
            nm: Some(2),
            is_primary: true,
            is_supplementary: false,
            r1_mapped: true,
            r2_mapped: true,
            pair_consistent: true,
            reference_start: 0,
            reference_end: 100,
            virus_target: None,
        }
    }

    fn virus_hit(as_score: i32) -> AlignHit {
        AlignHit {
            contig: "virus-alpha".to_string(),
            mapq: 30,
            as_score: Some(as_score),
            nm: Some(3),
            is_primary: true,
            is_supplementary: false,
            r1_mapped: true,
            r2_mapped: true,
            pair_consistent: true,
            reference_start: 0,
            reference_end: 100,
            virus_target: None,
        }
    }
}
