use crate::error::{Result, RvScreenError};
use crate::io::fastq::FragmentRecord;
use crate::types::{BundleManifest, ContigEntry};
use minimap2::{Aligner, Built, Mapping};
use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::Arc;

const DEFAULT_BEST_N: i32 = 5;
const DEFAULT_CAP_KALLOC: i64 = 512 * 1024 * 1024;
const TOTAL_KALLOC_BUDGET: u64 = 1_000_000_000;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum CompetitivePreset {
    #[default]
    SrConservative,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum ContigClass {
    Host,
    Virus,
    Decoy,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VirusTarget {
    pub accession: String,
    pub group: String,
    pub taxid: u64,
    pub virus_name: String,
    pub genome_length: u64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AlignHit {
    pub contig: String,
    pub mapq: u32,
    pub as_score: Option<i32>,
    pub nm: Option<i32>,
    pub is_primary: bool,
    pub is_supplementary: bool,
    pub r1_mapped: bool,
    pub r2_mapped: bool,
    pub pair_consistent: bool,
    pub reference_start: u64,
    pub reference_end: u64,
    pub virus_target: Option<VirusTarget>,
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct FragmentAlignResult {
    pub best_host_hit: Option<AlignHit>,
    pub best_virus_hit: Option<AlignHit>,
    pub secondary_hits: Vec<AlignHit>,
}

#[derive(Clone)]
pub struct CompetitiveAligner {
    aligner: Arc<Aligner<Built>>,
    contig_classes: BTreeMap<String, ContigClass>,
    virus_targets: BTreeMap<String, VirusTarget>,
    reference_input: PathBuf,
    manifest_path: PathBuf,
    preset: CompetitivePreset,
}

#[derive(Debug, Clone)]
struct ContigAggregate {
    class: ContigClass,
    contig: String,
    representative: Mapping,
    r1_mapped: bool,
    r2_mapped: bool,
    r1_best_mapq: u32,
    r2_best_mapq: u32,
    r1_best_as: Option<i32>,
    r2_best_as: Option<i32>,
    min_target_start: Option<u64>,
    max_target_end: Option<u64>,
    primary_supports: u8,
    supplementary_supports: u8,
}

impl CompetitiveAligner {
    pub fn new<P>(reference_bundle_path: P, preset: CompetitivePreset) -> Result<Self>
    where
        P: AsRef<Path>,
    {
        Self::new_with_threads(reference_bundle_path, preset, 1)
    }

    pub fn new_with_threads<P>(
        reference_bundle_path: P,
        preset: CompetitivePreset,
        threads: usize,
    ) -> Result<Self>
    where
        P: AsRef<Path>,
    {
        if threads == 0 {
            return Err(RvScreenError::validation(
                "threads",
                "--threads must be greater than zero",
            ));
        }

        let (reference_input, manifest_path) =
            resolve_bundle_inputs(reference_bundle_path.as_ref())?;
        let manifest = load_manifest(&manifest_path)?;
        let (contig_classes, virus_targets) = build_contig_indexes(&manifest)?;
        let aligner = build_competitive_aligner(&reference_input, preset, threads)?;

        Ok(Self {
            aligner: Arc::new(aligner),
            contig_classes,
            virus_targets,
            reference_input,
            manifest_path,
            preset,
        })
    }

    pub fn align_fragment(&self, fragment: &FragmentRecord) -> Result<FragmentAlignResult> {
        let (read1_mappings, read2_mappings) = self
            .aligner
            .map_pair(
                &fragment.r1_seq,
                &fragment.r2_seq,
                false,
                false,
                None,
                None,
                Some(fragment.fragment_key.as_bytes()),
            )
            .map_err(|reason| RvScreenError::alignment(&fragment.fragment_key, reason))?;

        let aggregates = self.aggregate_hits(&read1_mappings, &read2_mappings)?;

        Ok(FragmentAlignResult {
            best_host_hit: select_best_hit(&aggregates, ContigClass::Host, &self.virus_targets),
            best_virus_hit: select_best_hit(&aggregates, ContigClass::Virus, &self.virus_targets),
            secondary_hits: collect_secondary_hits(
                &read1_mappings,
                &read2_mappings,
                &aggregates,
                &self.virus_targets,
            ),
        })
    }

    pub fn reference_input(&self) -> &Path {
        &self.reference_input
    }

    pub fn manifest_path(&self) -> &Path {
        &self.manifest_path
    }

    pub fn preset(&self) -> CompetitivePreset {
        self.preset
    }

    fn aggregate_hits(
        &self,
        read1_mappings: &[Mapping],
        read2_mappings: &[Mapping],
    ) -> Result<BTreeMap<String, ContigAggregate>> {
        let mut aggregates = BTreeMap::new();

        for mapping in read1_mappings {
            self.observe_mapping(&mut aggregates, mapping, Mate::Read1)?;
        }
        for mapping in read2_mappings {
            self.observe_mapping(&mut aggregates, mapping, Mate::Read2)?;
        }

        Ok(aggregates)
    }

    fn observe_mapping(
        &self,
        aggregates: &mut BTreeMap<String, ContigAggregate>,
        mapping: &Mapping,
        mate: Mate,
    ) -> Result<()> {
        let Some(contig) = mapping.target_name.as_deref() else {
            return Ok(());
        };
        let class = *self.contig_classes.get(contig).ok_or_else(|| {
            RvScreenError::validation(
                "manifest",
                format!(
                    "mapping emitted contig `{contig}` that was not classified from the manifest"
                ),
            )
        })?;

        aggregates
            .entry(contig.to_string())
            .and_modify(|aggregate| aggregate.observe(mapping, mate))
            .or_insert_with(|| ContigAggregate::new(contig.to_string(), class, mapping, mate));

        Ok(())
    }
}

impl ContigAggregate {
    fn new(contig: String, class: ContigClass, mapping: &Mapping, mate: Mate) -> Self {
        let mut aggregate = Self {
            class,
            contig,
            representative: mapping.clone(),
            r1_mapped: false,
            r2_mapped: false,
            r1_best_mapq: 0,
            r2_best_mapq: 0,
            r1_best_as: None,
            r2_best_as: None,
            min_target_start: None,
            max_target_end: None,
            primary_supports: 0,
            supplementary_supports: 0,
        };
        aggregate.observe(mapping, mate);
        aggregate
    }

    fn observe(&mut self, mapping: &Mapping, mate: Mate) {
        match mate {
            Mate::Read1 => {
                self.r1_mapped = true;
                self.r1_best_mapq = self.r1_best_mapq.max(mapping.mapq);
                self.r1_best_as = max_optional_score(self.r1_best_as, alignment_score(mapping));
            }
            Mate::Read2 => {
                self.r2_mapped = true;
                self.r2_best_mapq = self.r2_best_mapq.max(mapping.mapq);
                self.r2_best_as = max_optional_score(self.r2_best_as, alignment_score(mapping));
            }
        }

        if mapping.is_primary {
            self.primary_supports = self.primary_supports.saturating_add(1);
        }
        if mapping.is_supplementary {
            self.supplementary_supports = self.supplementary_supports.saturating_add(1);
        }

        if let Some((target_start, target_end)) = mapping_interval(mapping) {
            self.min_target_start = Some(
                self.min_target_start
                    .map_or(target_start, |current| current.min(target_start)),
            );
            self.max_target_end = Some(
                self.max_target_end
                    .map_or(target_end, |current| current.max(target_end)),
            );
        }

        if compare_mapping_priority(mapping, &self.representative).is_gt() {
            self.representative = mapping.clone();
        }
    }

    fn mapped_mates(&self) -> u8 {
        u8::from(self.r1_mapped) + u8::from(self.r2_mapped)
    }

    fn pair_consistent(&self) -> bool {
        self.r1_mapped && self.r2_mapped
    }

    fn combined_mapq(&self) -> u64 {
        u64::from(self.r1_best_mapq) + u64::from(self.r2_best_mapq)
    }

    fn combined_as(&self) -> i32 {
        self.r1_best_as.unwrap_or(i32::MIN / 4) + self.r2_best_as.unwrap_or(i32::MIN / 4)
    }

    fn to_hit(&self, virus_target: Option<&VirusTarget>) -> AlignHit {
        let (reference_start, reference_end) = self.reference_interval();

        AlignHit {
            contig: self.contig.clone(),
            mapq: self.representative.mapq,
            as_score: alignment_score(&self.representative),
            nm: edit_distance(&self.representative),
            is_primary: self.representative.is_primary,
            is_supplementary: self.representative.is_supplementary,
            r1_mapped: self.r1_mapped,
            r2_mapped: self.r2_mapped,
            pair_consistent: self.pair_consistent(),
            reference_start,
            reference_end,
            virus_target: virus_target.cloned(),
        }
    }

    fn reference_interval(&self) -> (u64, u64) {
        match (self.min_target_start, self.max_target_end) {
            (Some(start), Some(end)) if end >= start => (start, end),
            _ => mapping_interval(&self.representative).unwrap_or((0, 0)),
        }
    }
}

#[derive(Debug, Clone, Copy)]
enum Mate {
    Read1,
    Read2,
}

fn resolve_bundle_inputs(bundle_path: &Path) -> Result<(PathBuf, PathBuf)> {
    let metadata = fs::metadata(bundle_path).map_err(|err| RvScreenError::io(bundle_path, err))?;

    if metadata.is_file() {
        let manifest_path = bundle_path
            .parent()
            .map(|parent| parent.join("manifest.json"))
            .ok_or_else(|| {
                RvScreenError::validation(
                    "reference_bundle",
                    format!(
                        "cannot infer manifest.json next to reference input `{}`",
                        bundle_path.display()
                    ),
                )
            })?;
        validate_manifest_path(&manifest_path)?;
        return Ok((bundle_path.to_path_buf(), manifest_path));
    }

    if !metadata.is_dir() {
        return Err(RvScreenError::validation(
            "reference_bundle",
            format!(
                "`{}` is neither a file nor a directory",
                bundle_path.display()
            ),
        ));
    }

    let manifest_path = bundle_path.join("manifest.json");
    validate_manifest_path(&manifest_path)?;
    let reference_input = discover_reference_input(bundle_path)?;
    Ok((reference_input, manifest_path))
}

fn validate_manifest_path(path: &Path) -> Result<()> {
    let metadata = fs::metadata(path).map_err(|err| RvScreenError::io(path, err))?;
    if !metadata.is_file() {
        return Err(RvScreenError::validation(
            "manifest",
            format!("`{}` is not a regular file", path.display()),
        ));
    }
    Ok(())
}

fn discover_reference_input(bundle_dir: &Path) -> Result<PathBuf> {
    let preferred_candidates = [
        bundle_dir.join("index/minimap2/composite.mmi"),
        bundle_dir.join("composite.fa"),
        bundle_dir.join("composite.mmi"),
        bundle_dir.join("mini_reference.fa"),
    ];

    for candidate in preferred_candidates {
        if candidate.is_file() {
            return Ok(candidate);
        }
    }

    let mut discovered = fs::read_dir(bundle_dir)
        .map_err(|err| RvScreenError::io(bundle_dir, err))?
        .filter_map(|entry| entry.ok().map(|entry| entry.path()))
        .filter(|path| path.is_file())
        .filter(|path| {
            path.extension()
                .and_then(|extension| extension.to_str())
                .is_some_and(|extension| matches!(extension, "fa" | "fasta" | "fna" | "mmi"))
        })
        .collect::<Vec<_>>();
    discovered.sort();
    discovered.dedup();

    match discovered.as_slice() {
        [single] => Ok(single.clone()),
        [] => Err(RvScreenError::validation(
            "reference_bundle",
            format!(
                "`{}` did not contain a composite reference input (.mmi/.fa/.fasta/.fna)",
                bundle_dir.display()
            ),
        )),
        many => Err(RvScreenError::validation(
            "reference_bundle",
            format!(
                "multiple candidate reference inputs found in `{}`: {}",
                bundle_dir.display(),
                many.iter()
                    .map(|path| path
                        .file_name()
                        .unwrap_or_default()
                        .to_string_lossy()
                        .into_owned())
                    .collect::<Vec<_>>()
                    .join(", ")
            ),
        )),
    }
}

fn load_manifest(path: &Path) -> Result<BundleManifest> {
    let bytes = fs::read(path).map_err(|err| RvScreenError::io(path, err))?;
    serde_json::from_slice(&bytes).map_err(|err| RvScreenError::parse(path, 1, err.to_string()))
}

fn build_contig_indexes(
    manifest: &BundleManifest,
) -> Result<(BTreeMap<String, ContigClass>, BTreeMap<String, VirusTarget>)> {
    let mut contig_classes = BTreeMap::new();
    let mut virus_targets = BTreeMap::new();
    for entry in &manifest.0 {
        let class = classify_contig(entry)?;
        if contig_classes.insert(entry.contig.clone(), class).is_some() {
            return Err(RvScreenError::validation(
                "manifest",
                format!(
                    "duplicate contig classification entry for `{}`",
                    entry.contig
                ),
            ));
        }

        if class == ContigClass::Virus {
            virus_targets.insert(
                entry.contig.clone(),
                VirusTarget {
                    accession: entry.accession.clone(),
                    group: entry.group.clone(),
                    taxid: entry.taxid,
                    virus_name: entry.virus_name.clone(),
                    genome_length: entry.genome_length,
                },
            );
        }
    }
    Ok((contig_classes, virus_targets))
}

fn classify_contig(entry: &ContigEntry) -> Result<ContigClass> {
    let group = normalize_metadata_value(&entry.group);
    let source_type = normalize_metadata_value(&entry.source_type);
    let from_group = classify_group(&group);
    let from_source_type = classify_source_type(&source_type);

    match (from_group, from_source_type) {
        (Some(left), Some(right)) if left != right => Err(RvScreenError::validation(
            "manifest",
            format!(
                "contig `{}` has conflicting classification metadata: group=`{}` source_type=`{}`",
                entry.contig, entry.group, entry.source_type
            ),
        )),
        (Some(class), _) | (_, Some(class)) => Ok(class),
        (None, None) => Err(RvScreenError::validation(
            "manifest",
            format!(
                "contig `{}` has unhandled source classification metadata: group=`{}` source_type=`{}`",
                entry.contig, entry.group, entry.source_type
            ),
        )),
    }
}

fn classify_group(group: &str) -> Option<ContigClass> {
    match group {
        "human" | "host" | "host-backbone" => Some(ContigClass::Host),
        "virus" | "viral" | "viral-panel" => Some(ContigClass::Virus),
        "decoy" | "background" | "background-decoy" => Some(ContigClass::Decoy),
        _ => None,
    }
}

fn classify_source_type(source_type: &str) -> Option<ContigClass> {
    match source_type {
        "host" | "host-genome" | "host-backbone" => Some(ContigClass::Host),
        "virus" | "virus-panel" | "viral-panel" | "refseq-virus" | "spike-in-virus" => {
            Some(ContigClass::Virus)
        }
        "decoy" | "decoy-panel" | "adapter-decoy" | "phix-decoy" | "vector-decoy" => {
            Some(ContigClass::Decoy)
        }
        _ => None,
    }
}

fn normalize_metadata_value(value: &str) -> String {
    value.trim().to_ascii_lowercase().replace(['_', ' '], "-")
}

fn build_competitive_aligner(
    reference_input: &Path,
    preset: CompetitivePreset,
    threads: usize,
) -> Result<Aligner<Built>> {
    let mut builder = match preset {
        CompetitivePreset::SrConservative => Aligner::builder().sr().with_cigar(),
    };
    builder.mapopt.best_n = DEFAULT_BEST_N;
    builder.mapopt.cap_kalloc = per_thread_kalloc_cap(threads);

    builder.with_index(reference_input, None).map_err(|reason| {
        RvScreenError::validation(
            "reference_bundle",
            format!(
                "failed to load minimap2 reference input `{}`: {reason}",
                reference_input.display()
            ),
        )
    })
}

fn per_thread_kalloc_cap(threads: usize) -> i64 {
    let threads = threads.max(1);
    let cap = TOTAL_KALLOC_BUDGET / u64::try_from(threads).unwrap_or(u64::MAX).max(1);
    let cap = cap.max(1);
    i64::try_from(cap).unwrap_or(DEFAULT_CAP_KALLOC)
}

fn select_best_hit(
    aggregates: &BTreeMap<String, ContigAggregate>,
    class: ContigClass,
    virus_targets: &BTreeMap<String, VirusTarget>,
) -> Option<AlignHit> {
    aggregates
        .values()
        .filter(|aggregate| aggregate.class == class)
        .max_by(|left, right| compare_aggregate_priority(left, right))
        .map(|aggregate| aggregate.to_hit(virus_targets.get(&aggregate.contig)))
}

fn collect_secondary_hits(
    read1_mappings: &[Mapping],
    read2_mappings: &[Mapping],
    aggregates: &BTreeMap<String, ContigAggregate>,
    virus_targets: &BTreeMap<String, VirusTarget>,
) -> Vec<AlignHit> {
    let mut secondary_hits = read1_mappings
        .iter()
        .chain(read2_mappings.iter())
        .filter(|mapping| !mapping.is_primary || mapping.is_supplementary)
        .filter_map(|mapping| {
            let contig = mapping.target_name.as_deref()?;
            let aggregate = aggregates.get(contig)?;
            Some(AlignHit {
                contig: contig.to_string(),
                mapq: mapping.mapq,
                as_score: alignment_score(mapping),
                nm: edit_distance(mapping),
                is_primary: mapping.is_primary,
                is_supplementary: mapping.is_supplementary,
                r1_mapped: aggregate.r1_mapped,
                r2_mapped: aggregate.r2_mapped,
                pair_consistent: aggregate.pair_consistent(),
                reference_start: mapping_interval(mapping).map_or(0, |(start, _)| start),
                reference_end: mapping_interval(mapping).map_or(0, |(_, end)| end),
                virus_target: virus_targets.get(contig).cloned(),
            })
        })
        .collect::<Vec<_>>();

    secondary_hits.sort_by(compare_align_hit_priority);
    secondary_hits
}

fn compare_aggregate_priority(left: &ContigAggregate, right: &ContigAggregate) -> Ordering {
    left.pair_consistent()
        .cmp(&right.pair_consistent())
        .then_with(|| left.mapped_mates().cmp(&right.mapped_mates()))
        .then_with(|| left.primary_supports.cmp(&right.primary_supports))
        .then_with(|| {
            right
                .supplementary_supports
                .cmp(&left.supplementary_supports)
        })
        .then_with(|| left.combined_as().cmp(&right.combined_as()))
        .then_with(|| left.combined_mapq().cmp(&right.combined_mapq()))
        .then_with(|| compare_mapping_priority(&left.representative, &right.representative))
        .then_with(|| right.contig.cmp(&left.contig))
}

fn compare_mapping_priority(left: &Mapping, right: &Mapping) -> Ordering {
    left.is_primary
        .cmp(&right.is_primary)
        .then_with(|| right.is_supplementary.cmp(&left.is_supplementary))
        .then_with(|| left.mapq.cmp(&right.mapq))
        .then_with(|| compare_optional_high(alignment_score(left), alignment_score(right)))
        .then_with(|| compare_optional_low(edit_distance(left), edit_distance(right)))
        .then_with(|| right.segment_id.cmp(&left.segment_id))
        .then_with(|| right.target_start.cmp(&left.target_start))
        .then_with(|| right.target_end.cmp(&left.target_end))
}

fn compare_align_hit_priority(left: &AlignHit, right: &AlignHit) -> Ordering {
    left.contig
        .cmp(&right.contig)
        .then_with(|| right.is_primary.cmp(&left.is_primary))
        .then_with(|| left.is_supplementary.cmp(&right.is_supplementary))
        .then_with(|| right.mapq.cmp(&left.mapq))
        .then_with(|| compare_optional_high(right.as_score, left.as_score))
        .then_with(|| {
            left.nm
                .unwrap_or(i32::MAX)
                .cmp(&right.nm.unwrap_or(i32::MAX))
        })
}

fn alignment_score(mapping: &Mapping) -> Option<i32> {
    mapping
        .alignment
        .as_ref()
        .and_then(|alignment| alignment.alignment_score)
}

fn edit_distance(mapping: &Mapping) -> Option<i32> {
    mapping.alignment.as_ref().map(|alignment| alignment.nm)
}

fn mapping_interval(mapping: &Mapping) -> Option<(u64, u64)> {
    let start = u64::try_from(mapping.target_start).ok()?;
    let end = u64::try_from(mapping.target_end).ok()?;
    (end >= start).then_some((start, end))
}

fn max_optional_score(left: Option<i32>, right: Option<i32>) -> Option<i32> {
    match (left, right) {
        (Some(left), Some(right)) => Some(left.max(right)),
        (Some(left), None) => Some(left),
        (None, Some(right)) => Some(right),
        (None, None) => None,
    }
}

fn compare_optional_high(left: Option<i32>, right: Option<i32>) -> Ordering {
    left.unwrap_or(i32::MIN).cmp(&right.unwrap_or(i32::MIN))
}

fn compare_optional_low(left: Option<i32>, right: Option<i32>) -> Ordering {
    right.unwrap_or(i32::MAX).cmp(&left.unwrap_or(i32::MAX))
}

#[cfg(test)]
#[path = "../../tests/testutil/mod.rs"]
mod testutil;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::fastq::FastqPairReader;
    use crate::types::ContigEntry;
    use std::fs;
    use tempfile::{tempdir, TempDir};
    use testutil::{
        generate_fastq_pair, generate_mini_reference, FastqPairConfig, ReadComponent,
        SyntheticSource, VirusSelector,
    };

    #[test]
    fn test_competitive_alignment_separates_host_and_virus_hits() {
        let reference_dir =
            generate_mini_reference().expect("mini reference generation should succeed");
        let aligner = CompetitiveAligner::new(&reference_dir, CompetitivePreset::SrConservative)
            .expect("competitive aligner should load mini reference bundle");

        assert!(aligner.reference_input().ends_with("mini_reference.fa"));
        assert!(aligner.manifest_path().ends_with("manifest.json"));

        let human_fragment = first_fragment_from_config(
            FastqPairConfig::new(1, 75, 101)
                .with_components(vec![ReadComponent::new(SyntheticSource::Human, 1.0)]),
        );
        let virus_fragment =
            first_fragment_from_config(FastqPairConfig::new(1, 75, 202).with_components(vec![
                ReadComponent::new(
                    SyntheticSource::Virus(VirusSelector::Accession("NC_SYNTHV1.1".to_string())),
                    1.0,
                ),
            ]));

        let human_result = aligner
            .align_fragment(&human_fragment)
            .expect("human fragment should align competitively");
        let virus_result = aligner
            .align_fragment(&virus_fragment)
            .expect("virus fragment should align competitively");

        let best_host = human_result
            .best_host_hit
            .as_ref()
            .expect("human fragment should retain a host hit");
        assert_eq!(best_host.contig, "mini_human_chr1");
        assert!(best_host.mapq > 0);
        assert_eq!(best_host.r1_mapped, true);
        assert_eq!(best_host.r2_mapped, true);
        assert!(best_host.pair_consistent);
        if let Some(virus_hit) = human_result.best_virus_hit.as_ref() {
            assert!(
                (best_host.mapq, best_host.as_score.unwrap_or(i32::MIN))
                    >= (virus_hit.mapq, virus_hit.as_score.unwrap_or(i32::MIN)),
                "human host hit should outrank any competing virus hit: host={best_host:?} virus={virus_hit:?}"
            );
        }

        let best_virus = virus_result
            .best_virus_hit
            .as_ref()
            .expect("virus fragment should retain a virus hit");
        assert_eq!(best_virus.contig, "mini_virus_alpha");
        assert!(best_virus.mapq > 0);
        assert_eq!(best_virus.r1_mapped, true);
        assert_eq!(best_virus.r2_mapped, true);
        assert!(best_virus.pair_consistent);
        if let Some(host_hit) = virus_result.best_host_hit.as_ref() {
            assert!(
                (best_virus.mapq, best_virus.as_score.unwrap_or(i32::MIN))
                    >= (host_hit.mapq, host_hit.as_score.unwrap_or(i32::MIN)),
                "virus hit should outrank any competing host hit: virus={best_virus:?} host={host_hit:?}"
            );
        }
    }

    #[test]
    fn test_competitive_alignment_retains_secondary_hits() {
        let bundle_dir =
            duplicate_virus_bundle().expect("duplicate virus bundle should be created");
        let aligner = CompetitiveAligner::new(bundle_dir.path(), CompetitivePreset::SrConservative)
            .expect("competitive aligner should load duplicate-virus bundle");
        let fragment = fragment_from_template(
            "secondary-fragment",
            "ACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCA",
            75,
        );

        let result = aligner
            .align_fragment(&fragment)
            .expect("duplicate-virus fragment should align competitively");

        let best_virus = result
            .best_virus_hit
            .as_ref()
            .expect("duplicate-virus fragment should keep a best virus hit");
        assert!(
            matches!(best_virus.contig.as_str(), "virus_dup_a" | "virus_dup_b"),
            "best virus hit should be one of the duplicated virus contigs: {best_virus:?}"
        );
        assert!(!result.secondary_hits.is_empty());
        let alternate = if best_virus.contig == "virus_dup_a" {
            "virus_dup_b"
        } else {
            "virus_dup_a"
        };
        assert!(result
            .secondary_hits
            .iter()
            .any(|hit| hit.contig == alternate));
        assert!(result
            .secondary_hits
            .iter()
            .all(|hit| !hit.is_primary || hit.is_supplementary));
    }

    fn first_fragment_from_config(config: FastqPairConfig) -> FragmentRecord {
        let tempdir = tempdir().expect("tempdir should be created");
        let (r1_path, r2_path) =
            generate_fastq_pair(&config.with_output_dir(tempdir.path().join("generated")))
                .expect("FASTQ generation should succeed");
        let mut reader =
            FastqPairReader::open(&r1_path, &r2_path).expect("FASTQ pair should open cleanly");
        reader
            .next()
            .expect("FASTQ pair should contain a first fragment")
            .expect("first fragment should parse cleanly")
    }

    fn duplicate_virus_bundle() -> std::io::Result<TempDir> {
        let tempdir = tempdir()?;
        let fasta_path = tempdir.path().join("composite.fa");
        let manifest_path = tempdir.path().join("manifest.json");

        let host_seq = "TTGGAACCTTGGAACCTTGGAACCTTGGAACCTTGGAACCTTGGAACCTTGGAACCTTGGAACCTTGGAACCTTGGAACCTTGGAACCTTGGAACCTTGGAACCTTGGAACCTTGGAA";
        let virus_seq = "ACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCA";

        fs::write(
            &fasta_path,
            format!(
                ">host_unique\n{host_seq}\n>virus_dup_a\n{virus_seq}\n>virus_dup_b\n{virus_seq}\n"
            ),
        )?;

        let manifest = BundleManifest(vec![
            duplicate_manifest_entry("host_unique", "HOST-001", 9606, "human", "host-genome"),
            duplicate_manifest_entry(
                "virus_dup_a",
                "NC_DUPVIRUS_A.1",
                12_345,
                "virus",
                "refseq-virus",
            ),
            duplicate_manifest_entry(
                "virus_dup_b",
                "NC_DUPVIRUS_B.1",
                12_345,
                "virus",
                "refseq-virus",
            ),
        ]);
        let manifest_json = serde_json::to_vec_pretty(&manifest)
            .map_err(|err| std::io::Error::other(err.to_string()))?;
        fs::write(&manifest_path, manifest_json)?;

        Ok(tempdir)
    }

    fn duplicate_manifest_entry(
        contig: &str,
        accession: &str,
        taxid: u64,
        group: &str,
        source_type: &str,
    ) -> ContigEntry {
        ContigEntry {
            contig: contig.to_string(),
            accession: accession.to_string(),
            taxid,
            virus_name: contig.to_string(),
            segment: None,
            group: group.to_string(),
            genome_length: 152,
            source_release: "task-11-test".to_string(),
            source_type: source_type.to_string(),
            masked_regions: Vec::new(),
        }
    }

    fn fragment_from_template(
        fragment_key: &str,
        template: &str,
        read_length: usize,
    ) -> FragmentRecord {
        let r1 = template[..read_length].as_bytes().to_vec();
        let r2 = reverse_complement(&template[template.len() - read_length..]).into_bytes();
        let qual = vec![b'I'; read_length];

        FragmentRecord {
            fragment_key: fragment_key.to_string(),
            r1_seq: r1,
            r1_qual: qual.clone(),
            r2_seq: r2,
            r2_qual: qual,
        }
    }

    fn reverse_complement(sequence: &str) -> String {
        sequence
            .chars()
            .rev()
            .map(|base| match base {
                'A' => 'T',
                'C' => 'G',
                'G' => 'C',
                'T' => 'A',
                other => other,
            })
            .collect()
    }

    #[test]
    fn competitive_aligner_rejects_zero_threads() {
        let reference_dir =
            generate_mini_reference().expect("mini reference generation should succeed");

        let error = CompetitiveAligner::new_with_threads(
            &reference_dir,
            CompetitivePreset::SrConservative,
            0,
        )
        .err()
        .expect("zero threads should be rejected");

        assert!(error
            .to_string()
            .contains("--threads must be greater than zero"));
    }

    #[test]
    fn per_thread_kalloc_cap_scales_with_threads() {
        assert_eq!(per_thread_kalloc_cap(1), 1_000_000_000);
        assert_eq!(per_thread_kalloc_cap(4), 250_000_000);
    }
}
