#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 1 ]; then
  echo "usage: $0 <fixture-root>" >&2
  exit 64
fi

fixture_root="$1"

# Safety guard: refuse to clobber a target that does not look like a
# previously-generated fixture. We only allow `rm -rf` when the path either
# does not exist yet or already contains the marker file we drop below.
if [ -e "${fixture_root}" ]; then
  if [ ! -f "${fixture_root}/.rvscreen_audit_fixture" ]; then
    echo "refusing to overwrite ${fixture_root}: missing .rvscreen_audit_fixture marker" >&2
    echo "(remove the directory manually if you really meant to point at it)" >&2
    exit 1
  fi
fi
rm -rf "${fixture_root}"
mkdir -p "${fixture_root}"
touch "${fixture_root}/.rvscreen_audit_fixture"

python3 - "${fixture_root}" <<'PY'
import json
import pathlib
import sys

root = pathlib.Path(sys.argv[1])
inputs = root / "inputs"
fastq = root / "fastq"
calibration = root / "calibration-profile"
inputs.mkdir(parents=True, exist_ok=True)
fastq.mkdir(parents=True, exist_ok=True)
calibration.mkdir(parents=True, exist_ok=True)

alphabet = "ACGT"
def dna(length: int, seed: int) -> str:
    state = seed & 0xFFFFFFFFFFFFFFFF
    out = []
    for _ in range(length):
        state = (state * 6364136223846793005 + 1442695040888963407) & 0xFFFFFFFFFFFFFFFF
        out.append(alphabet[(state >> 32) & 3])
    return "".join(out)

def wrap(seq: str, width: int = 80) -> str:
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))

records = {
    "mini_human_chr1": dna(1024, 0xA1B2_C3D4),
    "mini_virus_alpha": dna(512, 0x1100_2200),
    "mini_decoy_adapter": ("AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGTTAACCGGTTCCAAGGTT" * 4)[:224],
}

(inputs / "host.fa").write_text(f">mini_human_chr1\n{wrap(records['mini_human_chr1'])}\n")
(inputs / "virus.fa").write_text(f">mini_virus_alpha\n{wrap(records['mini_virus_alpha'])}\n")
(inputs / "decoy.fa").write_text(f">mini_decoy_adapter\n{wrap(records['mini_decoy_adapter'])}\n")

manifest = [
    {
        "contig": "mini_human_chr1",
        "accession": "CHM13-MINI-001",
        "taxid": 9606,
        "virus_name": "human_background",
        "segment": None,
        "group": "human",
        "genome_length": 1024,
        "source_release": "ci-2026-05",
        "source_type": "host-genome",
        "masked_regions": [],
    },
    {
        "contig": "mini_virus_alpha",
        "accession": "NC_SYNTHV1.1",
        "taxid": 2697049,
        "virus_name": "Synthetic virus alpha",
        "segment": "segment-A",
        "group": "virus",
        "genome_length": 512,
        "source_release": "ci-2026-05",
        "source_type": "refseq-virus",
        "masked_regions": [],
    },
    {
        "contig": "mini_decoy_adapter",
        "accession": "DECOY-MINI-001",
        "taxid": 32630,
        "virus_name": "adapter_decoy",
        "segment": None,
        "group": "decoy",
        "genome_length": 224,
        "source_release": "ci-2026-05",
        "source_type": "adapter-decoy",
        "masked_regions": [],
    },
]
(inputs / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
(inputs / "taxonomy.tsv").write_text(
    "taxid\tname\n"
    "9606\thuman_background\n"
    "2697049\tSynthetic virus alpha\n"
    "32630\tadapter_decoy\n"
)

def write_fastq(path: pathlib.Path, reads):
    with path.open("w") as handle:
        for idx, seq in enumerate(reads, 1):
            handle.write(f"@ci-read-{idx}\n{seq}\n+\n{'I' * len(seq)}\n")

read_len = 100
r1_reads = []
r2_reads = []
for idx in range(200):
    source = records["mini_virus_alpha"] if idx % 20 == 0 else records["mini_human_chr1"]
    offset = (idx * 7) % (len(source) - read_len - 20)
    r1_reads.append(source[offset:offset + read_len])
    r2_reads.append(source[offset + 20:offset + 20 + read_len])
write_fastq(fastq / "sample_R1.fastq", r1_reads)
write_fastq(fastq / "sample_R2.fastq", r2_reads)

(calibration / "release_gate.json").write_text(json.dumps({
    "backend_gate": {"status": "pass", "details": "CI generated audit fixture backend gate"},
    "reference_gate": {"status": "pass", "details": "CI generated audit fixture reference gate"},
    "specificity_gate": {"status": "pass", "details": "CI generated audit fixture specificity gate", "negative_samples": 1, "false_positives": 0},
    "sensitivity_gate": {"status": "pass", "details": "CI generated audit fixture sensitivity gate", "spike_in_detected": 1, "spike_in_total": 1},
}, indent=2) + "\n")
PY

cargo run --locked -- ref build \
  --host-fasta "${fixture_root}/inputs/host.fa" \
  --virus-fasta "${fixture_root}/inputs/virus.fa" \
  --decoy-fasta "${fixture_root}/inputs/decoy.fa" \
  --manifest "${fixture_root}/inputs/manifest.json" \
  --taxonomy "${fixture_root}/inputs/taxonomy.tsv" \
  --out "${fixture_root}/reference-bundle"

reference_version="$(python3 - "${fixture_root}/reference-bundle/bundle.toml" <<'PY'
import pathlib
import re
import sys
body = pathlib.Path(sys.argv[1]).read_text()
match = re.search(r'^version\s*=\s*"([^"]+)"', body, re.MULTILINE)
if not match:
    raise SystemExit("reference bundle version not found")
print(match.group(1))
PY
)"

cat > "${fixture_root}/calibration-profile/profile.toml" <<EOF
profile_id = "rvscreen_calib_ci_audit_fixture"
status = "release_candidate"
reference_bundle = "${reference_version}"
backend = "minimap2"
preset = "sr-conservative"
seed = 20260507
supported_input = ["fastq", "fastq.gz", "bam", "ubam", "cram"]
supported_read_type = ["illumina_pe_shortread"]
negative_control_required = false

[sampling]
mode = "representative"
round_mode = "absolute"
rounds = [50, 100]
max_rounds = 2

[fragment_rules]
min_mapq = 0
min_as_diff = 0
max_nm = 100
require_pair_consistency = true

[candidate_rules]
min_nonoverlap_fragments = 1
min_breadth = 0.0
max_background_ratio = 0.0
theta_pos_absolute = 3
max_ambiguous_fraction_for_positive = 0.20
min_positive_evidence_strength = "high"
weak_positive_enabled = true

[decision_rules]
theta_pos = 0.01
theta_neg = 0.0001
allow_indeterminate = true
theta_neg_fixed = 0.00001
positive_alpha_global = 0.05
positive_statistic_method = "exact_binomial_survival_bonferroni"
negative_cp_lod_max = 0.00001
noise_floor_fraction = 0.000001
EOF

cargo run --locked -- screen \
  --input "${fixture_root}/fastq/sample_R1.fastq" "${fixture_root}/fastq/sample_R2.fastq" \
  --reference-bundle "${fixture_root}/reference-bundle" \
  --calibration-profile "${fixture_root}/calibration-profile" \
  --out "${fixture_root}/report-bundle" \
  --mode representative \
  --threads 2
