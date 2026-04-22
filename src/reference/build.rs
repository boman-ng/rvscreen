use crate::error::{Result, RvScreenError};
use crate::types::{BundleManifest, BundleToml};
use minimap2::Aligner;
use sha2::{Digest, Sha256};
use std::collections::BTreeSet;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

const FASTA_LINE_BASES: usize = 80;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BuildReferenceBundleRequest {
    pub host_fasta: PathBuf,
    pub virus_fasta: PathBuf,
    pub decoy_fasta: Option<PathBuf>,
    pub manifest: PathBuf,
    pub taxonomy: PathBuf,
    pub out_dir: PathBuf,
}

#[derive(Debug, Clone, PartialEq)]
pub struct BuildReferenceBundleOutcome {
    pub bundle_dir: PathBuf,
    pub bundle: BundleToml,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct LayerSpec {
    label: &'static str,
    bundle_name: &'static str,
    path: PathBuf,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct FastaRecord {
    header: String,
    sequence: String,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct FaiRecord {
    header: String,
    sequence_len: usize,
    sequence_offset: u64,
    line_bases: usize,
    line_width: usize,
}

pub fn build_reference_bundle(
    request: &BuildReferenceBundleRequest,
) -> Result<BuildReferenceBundleOutcome> {
    let layers = collect_layers(request);
    for layer in &layers {
        validate_readable_file(&layer.path, layer.label)?;
    }
    validate_readable_file(&request.manifest, "manifest")?;
    validate_readable_file(&request.taxonomy, "taxonomy")?;
    validate_output_directory(&request.out_dir)?;

    let manifest = read_manifest(&request.manifest)?;
    let merged_records = read_layers(&layers)?;
    validate_manifest_headers(&manifest, &merged_records)?;

    fs::create_dir_all(&request.out_dir).map_err(|err| RvScreenError::io(&request.out_dir, err))?;
    let composite_path = request.out_dir.join("composite.fa");
    let fai_path = request.out_dir.join("composite.fa.fai");
    let manifest_out = request.out_dir.join("manifest.json");
    let taxonomy_out = request.out_dir.join("taxonomy.tsv");
    let index_dir = request.out_dir.join("index/minimap2");
    let mmi_path = index_dir.join("composite.mmi");
    let bundle_toml_path = request.out_dir.join("bundle.toml");
    let checksum_path = request.out_dir.join("checksum.sha256");

    let fai_records = write_composite_fasta(&composite_path, &merged_records)?;
    write_fai(&fai_path, &fai_records)?;
    fs::copy(&request.manifest, &manifest_out)
        .map_err(|err| RvScreenError::io(&manifest_out, err))?;
    fs::copy(&request.taxonomy, &taxonomy_out)
        .map_err(|err| RvScreenError::io(&taxonomy_out, err))?;

    fs::create_dir_all(&index_dir).map_err(|err| RvScreenError::io(&index_dir, err))?;
    build_minimap2_index(&composite_path, &mmi_path)?;

    let (created_at, date_stamp) = current_utc_timestamp()?;
    let bundle = BundleToml {
        version: format!("rvscreen_ref_{date_stamp}-r1"),
        created_at,
        included_layers: layers
            .iter()
            .map(|layer| layer.bundle_name.to_string())
            .collect(),
    };
    write_bundle_toml(&bundle_toml_path, &bundle)?;

    write_checksum_file(
        &checksum_path,
        &request.out_dir,
        &[
            bundle_toml_path,
            composite_path,
            fai_path,
            manifest_out,
            taxonomy_out,
            mmi_path,
        ],
    )?;

    Ok(BuildReferenceBundleOutcome {
        bundle_dir: request.out_dir.clone(),
        bundle,
    })
}

fn collect_layers(request: &BuildReferenceBundleRequest) -> Vec<LayerSpec> {
    let mut layers = vec![
        LayerSpec {
            label: "host_fasta",
            bundle_name: "host_backbone",
            path: request.host_fasta.clone(),
        },
        LayerSpec {
            label: "virus_fasta",
            bundle_name: "viral_panel",
            path: request.virus_fasta.clone(),
        },
    ];

    if let Some(decoy_fasta) = &request.decoy_fasta {
        layers.push(LayerSpec {
            label: "decoy_fasta",
            bundle_name: "decoy_panel",
            path: decoy_fasta.clone(),
        });
    }

    layers
}

fn validate_readable_file(path: &Path, field: &str) -> Result<()> {
    let metadata = fs::metadata(path).map_err(|err| RvScreenError::io(path, err))?;
    if !metadata.is_file() {
        return Err(RvScreenError::validation(
            field,
            format!("`{}` is not a regular file", path.display()),
        ));
    }
    if metadata.len() == 0 {
        return Err(RvScreenError::validation(
            field,
            format!("`{}` is empty", path.display()),
        ));
    }

    File::open(path).map_err(|err| RvScreenError::io(path, err))?;
    Ok(())
}

fn validate_output_directory(path: &Path) -> Result<()> {
    if path.exists() {
        if !path.is_dir() {
            return Err(RvScreenError::validation(
                "out",
                format!("`{}` exists and is not a directory", path.display()),
            ));
        }

        let mut entries = fs::read_dir(path).map_err(|err| RvScreenError::io(path, err))?;
        if entries
            .next()
            .transpose()
            .map_err(|err| RvScreenError::io(path, err))?
            .is_some()
        {
            return Err(RvScreenError::validation(
                "out",
                format!("output directory `{}` must be empty", path.display()),
            ));
        }
    }

    Ok(())
}

fn read_manifest(path: &Path) -> Result<BundleManifest> {
    let bytes = fs::read(path).map_err(|err| RvScreenError::io(path, err))?;
    serde_json::from_slice(&bytes).map_err(|err| RvScreenError::parse(path, 1, err.to_string()))
}

fn read_layers(layers: &[LayerSpec]) -> Result<Vec<FastaRecord>> {
    let mut merged = Vec::new();
    let mut seen = BTreeSet::new();

    for layer in layers {
        for record in read_fasta_records(&layer.path, layer.label)? {
            if !seen.insert(record.header.clone()) {
                return Err(RvScreenError::validation(
                    layer.label,
                    format!(
                        "duplicate FASTA header `{}` detected across input layers",
                        record.header
                    ),
                ));
            }
            merged.push(record);
        }
    }

    if merged.is_empty() {
        return Err(RvScreenError::validation(
            "reference_fasta",
            "no FASTA records were loaded from the supplied layers",
        ));
    }

    Ok(merged)
}

fn read_fasta_records(path: &Path, field: &str) -> Result<Vec<FastaRecord>> {
    let file = File::open(path).map_err(|err| RvScreenError::io(path, err))?;
    let reader = BufReader::new(file);
    let mut records = Vec::new();
    let mut current_header: Option<String> = None;
    let mut current_sequence = String::new();

    for (index, line) in reader.lines().enumerate() {
        let line_number = (index + 1) as u64;
        let line = line.map_err(|err| RvScreenError::io(path, err))?;
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }

        if let Some(rest) = trimmed.strip_prefix('>') {
            if let Some(header) = current_header.take() {
                if current_sequence.is_empty() {
                    return Err(RvScreenError::parse(
                        path,
                        line_number,
                        format!("FASTA record `{header}` has no sequence lines"),
                    ));
                }
                records.push(FastaRecord {
                    header,
                    sequence: std::mem::take(&mut current_sequence),
                });
            }

            let header = rest.trim();
            if header.is_empty() {
                return Err(RvScreenError::parse(
                    path,
                    line_number,
                    "empty FASTA header",
                ));
            }
            current_header = Some(header.to_string());
            continue;
        }

        if current_header.is_none() {
            return Err(RvScreenError::parse(
                path,
                line_number,
                format!("encountered sequence content before any FASTA header in field `{field}`"),
            ));
        }

        current_sequence.push_str(trimmed);
    }

    if let Some(header) = current_header {
        if current_sequence.is_empty() {
            return Err(RvScreenError::parse(
                path,
                1,
                format!("FASTA record `{header}` has no sequence lines"),
            ));
        }
        records.push(FastaRecord {
            header,
            sequence: current_sequence,
        });
    }

    if records.is_empty() {
        return Err(RvScreenError::validation(
            field,
            format!("`{}` did not contain any FASTA records", path.display()),
        ));
    }

    Ok(records)
}

fn validate_manifest_headers(
    manifest: &BundleManifest,
    fasta_records: &[FastaRecord],
) -> Result<()> {
    let fasta_headers: BTreeSet<_> = fasta_records
        .iter()
        .map(|record| record.header.clone())
        .collect();
    let manifest_headers: BTreeSet<_> = manifest
        .0
        .iter()
        .map(|entry| entry.contig.clone())
        .collect();

    if fasta_headers.len() != fasta_records.len() {
        return Err(RvScreenError::validation(
            "composite_fasta",
            "duplicate FASTA headers detected in merged composite",
        ));
    }
    if manifest_headers.len() != manifest.0.len() {
        return Err(RvScreenError::validation(
            "manifest",
            "manifest contains duplicate contig names",
        ));
    }

    if manifest_headers != fasta_headers {
        let missing_in_fasta: Vec<_> = manifest_headers
            .difference(&fasta_headers)
            .cloned()
            .collect();
        let missing_in_manifest: Vec<_> = fasta_headers
            .difference(&manifest_headers)
            .cloned()
            .collect();
        let mut details = Vec::new();
        if !missing_in_fasta.is_empty() {
            details.push(format!(
                "manifest-only contigs: {}",
                missing_in_fasta.join(", ")
            ));
        }
        if !missing_in_manifest.is_empty() {
            details.push(format!(
                "fasta-only contigs: {}",
                missing_in_manifest.join(", ")
            ));
        }
        return Err(RvScreenError::validation(
            "manifest",
            format!(
                "manifest contig names must exactly match FASTA headers; {}",
                details.join("; ")
            ),
        ));
    }

    Ok(())
}

fn write_composite_fasta(path: &Path, records: &[FastaRecord]) -> Result<Vec<FaiRecord>> {
    let file = File::create(path).map_err(|err| RvScreenError::io(path, err))?;
    let mut writer = BufWriter::new(file);
    let mut position = 0u64;
    let mut fai_records = Vec::with_capacity(records.len());

    for record in records {
        writeln!(writer, ">{header}", header = record.header)
            .map_err(|err| RvScreenError::io(path, err))?;
        position += 1 + record.header.len() as u64 + 1;

        let sequence_offset = position;
        let line_bases = record.sequence.len().min(FASTA_LINE_BASES);
        let line_width = line_bases + 1;
        fai_records.push(FaiRecord {
            header: record.header.clone(),
            sequence_len: record.sequence.len(),
            sequence_offset,
            line_bases,
            line_width,
        });

        for chunk in record.sequence.as_bytes().chunks(FASTA_LINE_BASES) {
            writer
                .write_all(chunk)
                .and_then(|_| writer.write_all(b"\n"))
                .map_err(|err| RvScreenError::io(path, err))?;
            position += chunk.len() as u64 + 1;
        }
    }

    writer.flush().map_err(|err| RvScreenError::io(path, err))?;
    Ok(fai_records)
}

fn write_fai(path: &Path, records: &[FaiRecord]) -> Result<()> {
    let file = File::create(path).map_err(|err| RvScreenError::io(path, err))?;
    let mut writer = BufWriter::new(file);

    for record in records {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}",
            record.header,
            record.sequence_len,
            record.sequence_offset,
            record.line_bases,
            record.line_width
        )
        .map_err(|err| RvScreenError::io(path, err))?;
    }

    writer.flush().map_err(|err| RvScreenError::io(path, err))
}

fn build_minimap2_index(composite_fasta: &Path, mmi_path: &Path) -> Result<()> {
    let output = mmi_path.to_str().ok_or_else(|| {
        RvScreenError::validation(
            "reference_fasta",
            format!("index path `{}` is not valid UTF-8", mmi_path.display()),
        )
    })?;

    Aligner::builder()
        .sr()
        .with_index(composite_fasta, Some(output))
        .map(|_| ())
        .map_err(|reason| {
            RvScreenError::validation(
                "reference_fasta",
                format!(
                    "failed to build minimap2 index `{}` from `{}`: {reason}",
                    mmi_path.display(),
                    composite_fasta.display()
                ),
            )
        })
}

fn write_bundle_toml(path: &Path, bundle: &BundleToml) -> Result<()> {
    let toml = toml::to_string_pretty(bundle)
        .map_err(|err| RvScreenError::bundle(&bundle.version, err.to_string()))?;
    fs::write(path, toml).map_err(|err| RvScreenError::io(path, err))
}

fn write_checksum_file(checksum_path: &Path, bundle_dir: &Path, files: &[PathBuf]) -> Result<()> {
    let file = File::create(checksum_path).map_err(|err| RvScreenError::io(checksum_path, err))?;
    let mut writer = BufWriter::new(file);

    for path in files {
        let digest = sha256_file(path)?;
        let relative = path.strip_prefix(bundle_dir).map_err(|_| {
            RvScreenError::validation(
                "checksum",
                format!(
                    "`{}` is not inside `{}`",
                    path.display(),
                    bundle_dir.display()
                ),
            )
        })?;
        writeln!(writer, "{digest}  {}", relative.display())
            .map_err(|err| RvScreenError::io(checksum_path, err))?;
    }

    writer
        .flush()
        .map_err(|err| RvScreenError::io(checksum_path, err))
}

fn sha256_file(path: &Path) -> Result<String> {
    let file = File::open(path).map_err(|err| RvScreenError::io(path, err))?;
    let mut reader = BufReader::new(file);
    let mut hasher = Sha256::new();
    let mut buffer = [0u8; 8192];

    loop {
        let read = reader
            .read(&mut buffer)
            .map_err(|err| RvScreenError::io(path, err))?;
        if read == 0 {
            break;
        }
        hasher.update(&buffer[..read]);
    }

    Ok(format!("{:x}", hasher.finalize()))
}

fn current_utc_timestamp() -> Result<(String, String)> {
    let seconds = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_err(|err| RvScreenError::validation("time", err.to_string()))?
        .as_secs() as i64;
    let (year, month, day, hour, minute, second) = unix_to_utc(seconds);
    Ok((
        format!("{year:04}-{month:02}-{day:02}T{hour:02}:{minute:02}:{second:02}Z"),
        format!("{year:04}.{month:02}.{day:02}"),
    ))
}

fn unix_to_utc(seconds: i64) -> (i64, i64, i64, i64, i64, i64) {
    let days = seconds.div_euclid(86_400);
    let seconds_of_day = seconds.rem_euclid(86_400);
    let (year, month, day) = civil_from_days(days);
    let hour = seconds_of_day / 3_600;
    let minute = (seconds_of_day % 3_600) / 60;
    let second = seconds_of_day % 60;
    (year, month, day, hour, minute, second)
}

fn civil_from_days(days_since_epoch: i64) -> (i64, i64, i64) {
    let z = days_since_epoch + 719_468;
    let era = if z >= 0 { z } else { z - 146_096 } / 146_097;
    let doe = z - era * 146_097;
    let yoe = (doe - doe / 1_460 + doe / 36_524 - doe / 146_096) / 365;
    let y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let day = doy - (153 * mp + 2) / 5 + 1;
    let month = mp + if mp < 10 { 3 } else { -9 };
    let year = y + if month <= 2 { 1 } else { 0 };
    (year, month, day)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::{BundleManifest, BundleToml, ContigEntry};
    use std::fs;
    use tempfile::tempdir;

    #[test]
    fn test_build_reference_bundle_success() {
        let tmp = tempdir().expect("temp dir should be created");
        let fixtures = write_reference_fixtures(tmp.path(), false);
        let out_dir = tmp.path().join("bundle");

        let outcome = build_reference_bundle(&BuildReferenceBundleRequest {
            host_fasta: fixtures.host_fasta.clone(),
            virus_fasta: fixtures.virus_fasta.clone(),
            decoy_fasta: Some(fixtures.decoy_fasta.clone()),
            manifest: fixtures.manifest.clone(),
            taxonomy: fixtures.taxonomy.clone(),
            out_dir: out_dir.clone(),
        })
        .expect("bundle build should succeed");

        assert_eq!(outcome.bundle_dir, out_dir);
        assert!(outcome.bundle.version.starts_with("rvscreen_ref_"));
        assert_eq!(
            outcome.bundle.included_layers,
            vec!["host_backbone", "viral_panel", "decoy_panel"]
        );

        let bundle_toml: BundleToml = toml::from_str(
            &fs::read_to_string(out_dir.join("bundle.toml")).expect("bundle TOML should exist"),
        )
        .expect("bundle TOML should deserialize");
        assert_eq!(bundle_toml, outcome.bundle);

        let manifest: BundleManifest = serde_json::from_str(
            &fs::read_to_string(out_dir.join("manifest.json"))
                .expect("manifest copy should be readable"),
        )
        .expect("manifest copy should deserialize");
        assert_eq!(manifest.0.len(), 6);

        let taxonomy = fs::read_to_string(out_dir.join("taxonomy.tsv"))
            .expect("taxonomy copy should be readable");
        assert!(taxonomy.contains("9606\thuman_background"));
        assert!(taxonomy.contains("32630\tadapter_decoy"));

        let headers = fasta_headers(&out_dir.join("composite.fa"));
        assert_eq!(
            headers,
            vec![
                "host_chr1",
                "host_chr2",
                "host_chr3",
                "virus_alpha",
                "virus_beta",
                "decoy_adapter",
            ]
        );

        let fai =
            fs::read_to_string(out_dir.join("composite.fa.fai")).expect("fai should be readable");
        assert_eq!(fai.lines().count(), 6);

        let index_metadata =
            fs::metadata(out_dir.join("index/minimap2/composite.mmi")).expect("mmi should exist");
        assert!(index_metadata.len() > 0, "mmi should not be empty");

        let checksum = fs::read_to_string(out_dir.join("checksum.sha256"))
            .expect("checksum file should be readable");
        assert_eq!(checksum.lines().count(), 6);
        for relative in [
            "bundle.toml",
            "composite.fa",
            "composite.fa.fai",
            "manifest.json",
            "taxonomy.tsv",
            "index/minimap2/composite.mmi",
        ] {
            assert!(
                checksum.contains(relative),
                "checksum file should include {relative}: {checksum}"
            );
        }
    }

    #[test]
    fn test_build_reference_bundle_rejects_manifest_mismatch() {
        let tmp = tempdir().expect("temp dir should be created");
        let fixtures = write_reference_fixtures(tmp.path(), true);
        let err = build_reference_bundle(&BuildReferenceBundleRequest {
            host_fasta: fixtures.host_fasta,
            virus_fasta: fixtures.virus_fasta,
            decoy_fasta: Some(fixtures.decoy_fasta),
            manifest: fixtures.manifest,
            taxonomy: fixtures.taxonomy,
            out_dir: tmp.path().join("bundle"),
        })
        .expect_err("manifest mismatch should fail");

        let message = err.to_string();
        assert!(message.contains("manifest contig names must exactly match FASTA headers"));
        assert!(message.contains("ghost_contig"));
    }

    struct FixturePaths {
        host_fasta: PathBuf,
        virus_fasta: PathBuf,
        decoy_fasta: PathBuf,
        manifest: PathBuf,
        taxonomy: PathBuf,
    }

    fn write_reference_fixtures(base: &Path, mismatch: bool) -> FixturePaths {
        let host_fasta = base.join("host.fa");
        let virus_fasta = base.join("virus.fa");
        let decoy_fasta = base.join("decoy.fa");
        let manifest = base.join("manifest.json");
        let taxonomy = base.join("taxonomy.tsv");

        fs::write(
            &host_fasta,
            concat!(
                ">host_chr1\n",
                "ACGTACGTACGTACGTACGTACGTACGTACGT\n",
                ">host_chr2\n",
                "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG\n",
                ">host_chr3\n",
                "GATCGATCGATCGATCGATCGATCGATCGATC\n",
            ),
        )
        .expect("host FASTA should be written");
        fs::write(
            &virus_fasta,
            concat!(
                ">virus_alpha\n",
                "AACCGGTTAACCGGTTAACCGGTTAACCGGTT\n",
                ">virus_beta\n",
                "GGCCAATTGGCCAATTGGCCAATTGGCCAATT\n",
            ),
        )
        .expect("virus FASTA should be written");
        fs::write(
            &decoy_fasta,
            concat!(">decoy_adapter\n", "AGATCGGAAGAGCACACGTCTGAACTCCAGTC\n",),
        )
        .expect("decoy FASTA should be written");

        let manifest_entries = vec![
            contig(
                "host_chr1",
                "HOST-001",
                9606,
                "human_background",
                "human",
                32,
                "host",
            ),
            contig(
                "host_chr2",
                "HOST-002",
                9606,
                "human_background",
                "human",
                32,
                "host",
            ),
            contig(
                "host_chr3",
                "HOST-003",
                9606,
                "human_background",
                "human",
                32,
                "host",
            ),
            contig(
                "virus_alpha",
                "NC_ALPHA.1",
                11103,
                "Virus alpha",
                "virus",
                32,
                "virus",
            ),
            contig(
                "virus_beta",
                "NC_BETA.1",
                11104,
                "Virus beta",
                "virus",
                32,
                "virus",
            ),
            contig(
                if mismatch {
                    "ghost_contig"
                } else {
                    "decoy_adapter"
                },
                "DECOY-001",
                32630,
                "adapter_decoy",
                "decoy",
                32,
                "decoy",
            ),
        ];
        let manifest_json = serde_json::to_string_pretty(&BundleManifest(manifest_entries))
            .expect("manifest should serialize");
        fs::write(&manifest, manifest_json).expect("manifest should be written");

        fs::write(
            &taxonomy,
            concat!(
                "taxid\tname\n",
                "9606\thuman_background\n",
                "11103\tVirus alpha\n",
                "11104\tVirus beta\n",
                "32630\tadapter_decoy\n",
            ),
        )
        .expect("taxonomy should be written");

        FixturePaths {
            host_fasta,
            virus_fasta,
            decoy_fasta,
            manifest,
            taxonomy,
        }
    }

    fn contig(
        contig: &str,
        accession: &str,
        taxid: u64,
        virus_name: &str,
        group: &str,
        genome_length: u64,
        source_type: &str,
    ) -> ContigEntry {
        ContigEntry {
            contig: contig.to_string(),
            accession: accession.to_string(),
            taxid,
            virus_name: virus_name.to_string(),
            segment: None,
            group: group.to_string(),
            genome_length,
            source_release: "test-2026-04".to_string(),
            source_type: source_type.to_string(),
            masked_regions: Vec::new(),
        }
    }

    fn fasta_headers(path: &Path) -> Vec<String> {
        fs::read_to_string(path)
            .expect("composite FASTA should be readable")
            .lines()
            .filter_map(|line| line.strip_prefix('>').map(str::to_string))
            .collect()
    }
}
