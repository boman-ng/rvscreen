#!/usr/bin/env bash

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT_DIR="${1:-${ROOT_DIR}/.local/reference-data}"

GRCH38_BASE="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38"
GRCH38_ALIGN_BASE="${GRCH38_BASE}/seqs_for_alignment_pipelines.ucsc_ids"
GRCH38_FASTA_URL="${GRCH38_ALIGN_BASE}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
GRCH38_REPORT_URL="${GRCH38_BASE}/GCA_000001405.15_GRCh38_assembly_report.txt"
GRCH38_SOURCE_RELEASE="GCA_000001405.15_GRCh38_no_alt_analysis_set"

VIRAL_BASE_URL="https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/"
VIRAL_SOURCE_RELEASE="RefSeq_viral_current"
UNIVEC_URL="https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec_Core"
UNIVEC_SOURCE_RELEASE="UniVec_Core_current"
PHIX_URL="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_001422.1&rettype=fasta&retmode=text"
PHIX_SOURCE_RELEASE="NC_001422.1"
ADAPTER_SOURCE_RELEASE="Illumina_adapter_set_embedded"

HOST_DIR="${OUT_DIR}/host"
VIRAL_DIR="${OUT_DIR}/viral"
DECOY_DIR="${OUT_DIR}/decoy"
META_DIR="${OUT_DIR}/metadata"
INPUT_DIR="${OUT_DIR}/bundle-inputs"
INPUT_META_DIR="${INPUT_DIR}/metadata"

HOST_FASTA_OUT="${INPUT_DIR}/host.fa"
VIRUS_FASTA_OUT="${INPUT_DIR}/virus.fa"
DECOY_FASTA_OUT="${INPUT_DIR}/decoy.fa"
MANIFEST_OUT="${INPUT_DIR}/manifest.json"
TAXONOMY_OUT="${INPUT_DIR}/taxonomy.tsv"

HOST_ROWS="${INPUT_META_DIR}/host_rows.tsv"
VIRUS_ROWS="${INPUT_META_DIR}/virus_rows.tsv"
DECOY_ROWS="${INPUT_META_DIR}/decoy_rows.tsv"

mkdir -p "${HOST_DIR}" "${VIRAL_DIR}" "${DECOY_DIR}" "${META_DIR}" "${INPUT_DIR}" "${INPUT_META_DIR}"

if [[ -t 1 && -z "${NO_COLOR:-}" ]]; then
  RESET='\033[0m'
  BOLD='\033[1m'
  DIM='\033[2m'
  BLUE='\033[34m'
  CYAN='\033[36m'
  GREEN='\033[32m'
  YELLOW='\033[33m'
  RED='\033[31m'
else
  RESET=''
  BOLD=''
  DIM=''
  BLUE=''
  CYAN=''
  GREEN=''
  YELLOW=''
  RED=''
fi

log_rule() {
  printf '%b\n' "${DIM}============================================================${RESET}"
}

log_banner() {
  log_rule
  printf '%b\n' "${BOLD}${BLUE}RVScreen Reference Test Input Prep${RESET}"
  printf '%b\n' "${DIM}One-stop workflow: discover -> download -> organize${RESET}"
  printf '%b\n' "${DIM}Output root: ${OUT_DIR}${RESET}"
  log_rule
}

log_section() {
  printf '\n%b[%s]%b %s\n' "${BOLD}${CYAN}" "$1" "${RESET}" "$2"
}

log_info() {
  printf '%b\n' "  ${DIM}- ${RESET}$1"
}

log_ok() {
  printf '%b\n' "  ${GREEN}ok${RESET}  $1"
}

log_warn() {
  printf '%b\n' "  ${YELLOW}warn${RESET} $1"
}

log_error() {
  printf '%b\n' "  ${RED}error${RESET} $1" >&2
}

fail() {
  log_error "$1"
  exit 1
}

remote_content_length() {
  local url="$1"
  local headers

  if command -v curl >/dev/null 2>&1; then
    headers="$(curl --silent --show-error --location --head "${url}" 2>/dev/null || true)"
  elif command -v wget >/dev/null 2>&1; then
    headers="$(wget --server-response --spider "${url}" 2>&1 || true)"
  else
    printf '\n'
    return 0
  fi

  printf '%s\n' "${headers}" \
    | awk '
        BEGIN {IGNORECASE=1}
        /^HTTP/ {status=$2}
        /^content-length:/ {gsub("\r", "", $2); len=$2}
        END {if (status ~ /^(200|206)$/ && len ~ /^[1-9][0-9]*$/) print len}
      '
}

is_download_complete() {
  local destination="$1"
  local expected_size="$2"

  if [[ ! -f "${destination}" ]]; then
    return 1
  fi

  if [[ -z "${expected_size}" || "${expected_size}" == "0" ]]; then
    [[ -s "${destination}" ]]
    return
  fi

  local actual_size
  actual_size="$(wc -c < "${destination}" | tr -d '[:space:]')"
  [[ "${actual_size}" == "${expected_size}" ]]
}

is_phix_complete() {
  local destination="$1"
  [[ -s "${destination}" ]] && grep -q '^>NC_001422\.1 ' "${destination}"
}

download() {
  local url="$1"
  local destination="$2"
  local expected_size

  expected_size="$(remote_content_length "${url}")"

  if is_download_complete "${destination}" "${expected_size}"; then
    log_ok "reuse complete download: ${destination}"
    return 0
  fi

  log_info "fetch ${url}"

  if command -v curl >/dev/null 2>&1; then
    curl --fail --location --retry 3 --retry-delay 2 --continue-at - --output "${destination}" "${url}"
  elif command -v wget >/dev/null 2>&1; then
    wget -c -O "${destination}" "${url}"
  else
    fail "curl or wget is required"
  fi

  if ! is_download_complete "${destination}" "${expected_size}"; then
    fail "download incomplete for ${destination}"
  fi

  log_ok "ready: ${destination}"
}

download_phix() {
  local destination="$1"

  if is_phix_complete "${destination}"; then
    log_ok "reuse complete download: ${destination}"
    return 0
  fi

  log_info "fetch ${PHIX_URL}"

  if command -v curl >/dev/null 2>&1; then
    curl --fail --location --retry 3 --retry-delay 2 --output "${destination}" "${PHIX_URL}"
  elif command -v wget >/dev/null 2>&1; then
    wget -O "${destination}" "${PHIX_URL}"
  else
    fail "curl or wget is required"
  fi

  if ! is_phix_complete "${destination}"; then
    fail "invalid PhiX FASTA content in ${destination}"
  fi

  log_ok "ready: ${destination}"
}

require_tool() {
  local tool="$1"
  command -v "${tool}" >/dev/null 2>&1 || fail "required tool not found: ${tool}"
}

needs_refresh() {
  local output="$1"
  shift

  if [[ ! -s "${output}" ]]; then
    return 0
  fi

  local source
  for source in "$@"; do
    if [[ ! -e "${source}" || "${source}" -nt "${output}" ]]; then
      return 0
    fi
  done

  return 1
}

write_adapter_fasta() {
  local destination="$1"

  cat > "${destination}" <<'EOF'
>Illumina_TruSeq_Adapter_Read1
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
>Illumina_TruSeq_Adapter_Read2
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
>Illumina_Nextera_Transposase_Read1
CTGTCTCTTATACACATCT
>Illumina_Nextera_Transposase_Read2
AGATGTGTATAAGAGACAG
EOF
}

discover_viral_urls() {
  local index_html
  index_html="$(mktemp)"
  trap 'rm -f "${index_html}"' RETURN

  download "${VIRAL_BASE_URL}" "${index_html}"

  mapfile -t viral_urls < <(
    grep -Eo 'viral\.[0-9]+\.[0-9]+\.genomic\.fna\.gz' "${index_html}" \
      | sort -u \
      | sed "s#^#${VIRAL_BASE_URL}#"
  )

  if [[ ${#viral_urls[@]} -eq 0 ]]; then
    fail "no viral genomic FASTA URLs discovered from ${VIRAL_BASE_URL}"
  fi

  printf '%s\n' "${viral_urls[@]}" > "${META_DIR}/viral_genomic_urls.txt"
  log_ok "discovered ${#viral_urls[@]} viral genomic FASTA volume(s)"
}

normalize_host_fasta() {
  local source="$1"
  local fasta_out="$2"
  local rows_out="$3"

  : > "${fasta_out}"
  : > "${rows_out}"

  gzip -dc "${source}" | awk -v fasta="${fasta_out}" -v rows="${rows_out}" -v release="${GRCH38_SOURCE_RELEASE}" '
    function flush_record() {
      if (contig == "") {
        return
      }
      printf "%s\t%s\t9606\thuman_background\thuman\t%d\t%s\thost_backbone\n", contig, accession, seq_len, release >> rows
    }
    /^>/ {
      flush_record()
      raw = substr($0, 2)
      split(raw, parts, /[[:space:]]+/)
      contig = parts[1]
      accession = contig
      if (match(raw, /AC:([^[:space:]]+)/, capture)) {
        accession = capture[1]
      }
      print ">" contig >> fasta
      seq_len = 0
      next
    }
    {
      gsub(/[[:space:]]/, "", $0)
      if ($0 == "") {
        next
      }
      line = toupper($0)
      print line >> fasta
      seq_len += length(line)
    }
    END {
      flush_record()
    }
  '
}

normalize_viral_fasta() {
  local fasta_out="$1"
  local rows_out="$2"
  shift 2

  : > "${fasta_out}"
  : > "${rows_out}"

  local source
  for source in "$@"; do
    gzip -dc "${source}" | awk -v fasta="${fasta_out}" -v rows="${rows_out}" -v release="${VIRAL_SOURCE_RELEASE}" '
      function flush_record() {
        if (contig == "") {
          return
        }
        printf "%s\t%s\t10239\t%s\t%s\t%d\t%s\trefseq-virus\n", contig, accession, virus_name, group_name, seq_len, release >> rows
      }
      /^>/ {
        flush_record()
        raw = substr($0, 2)
        split(raw, parts, /[[:space:]]+/)
        accession = parts[1]
        contig = accession
        virus_name = raw
        sub(/^[^[:space:]]+[[:space:]]*/, "", virus_name)
        if (virus_name == "") {
          virus_name = accession
        }
        group_name = accession
        print ">" accession >> fasta
        seq_len = 0
        next
      }
      {
        gsub(/[[:space:]]/, "", $0)
        if ($0 == "") {
          next
        }
        line = toupper($0)
        print line >> fasta
        seq_len += length(line)
      }
      END {
        flush_record()
      }
    '
  done
}

normalize_plain_fasta() {
  local source="$1"
  local fasta_out="$2"
  local rows_out="$3"
  local taxid="$4"
  local group_name="$5"
  local source_release="$6"
  local source_type="$7"
  local append_mode="$8"

  if [[ "${append_mode}" == "truncate" ]]; then
    : > "${fasta_out}"
    : > "${rows_out}"
  fi

  awk -v fasta="${fasta_out}" -v rows="${rows_out}" -v taxid="${taxid}" -v group_name="${group_name}" -v release="${source_release}" -v source_type="${source_type}" '
    function flush_record() {
      if (contig == "") {
        return
      }
      printf "%s\t%s\t%s\t%s\t%s\t%d\t%s\t%s\n", contig, accession, taxid, virus_name, group_name, seq_len, release, source_type >> rows
    }
    /^>/ {
      flush_record()
      raw = substr($0, 2)
      split(raw, parts, /[[:space:]]+/)
      accession = parts[1]
      contig = accession
      virus_name = raw
      sub(/^[^[:space:]]+[[:space:]]*/, "", virus_name)
      if (virus_name == "") {
        virus_name = accession
      }
      print ">" contig >> fasta
      seq_len = 0
      next
    }
    {
      gsub(/[[:space:]]/, "", $0)
      if ($0 == "") {
        next
      }
      line = toupper($0)
      print line >> fasta
      seq_len += length(line)
    }
    END {
      flush_record()
    }
  ' "${source}"
}

dedupe_organized_layer() {
  local fasta_path="$1"
  local rows_path="$2"
  local label="$3"
  local tmp_rows
  local tmp_fasta
  local tmp_dups
  local removed_count
  local removed_preview

  tmp_rows="$(mktemp)"
  tmp_fasta="$(mktemp)"
  tmp_dups="$(mktemp)"

  awk -F '\t' -v dup_file="${tmp_dups}" '
    !seen[$1]++ {
      print
      next
    }
    {
      removed += 1
      if (!dup_seen[$1]++) {
        print $1 >> dup_file
      }
    }
    END {
      printf "%d\n", removed + 0 > "/dev/stderr"
    }
  ' "${rows_path}" > "${tmp_rows}" 2> "${tmp_dups}.count"

  removed_count="$(tr -d '[:space:]' < "${tmp_dups}.count")"

  awk -v keep_rows="${tmp_rows}" '
    BEGIN {
      while ((getline < keep_rows) > 0) {
        split($0, fields, "\t")
        keep[fields[1]] = 1
      }
    }
    /^>/ {
      contig = substr($0, 2)
      split(contig, parts, /[[:space:]]+/)
      contig = parts[1]
      emit = keep[contig] && !written[contig]++
    }
    emit { print }
  ' "${fasta_path}" > "${tmp_fasta}"

  mv "${tmp_rows}" "${rows_path}"
  mv "${tmp_fasta}" "${fasta_path}"

  if [[ "${removed_count}" != "0" ]]; then
    removed_preview="$(paste -sd ', ' "${tmp_dups}" | cut -c1-180)"
    log_warn "removed ${removed_count} duplicate contig(s) from ${label}: ${removed_preview}"
  fi

  rm -f "${tmp_dups}" "${tmp_dups}.count"
}

exclude_contigs_from_layer() {
  local fasta_path="$1"
  local rows_path="$2"
  local excluded_rows_path="$3"
  local label="$4"
  local tmp_rows
  local tmp_fasta
  local tmp_removed
  local removed_count
  local removed_preview

  tmp_rows="$(mktemp)"
  tmp_fasta="$(mktemp)"
  tmp_removed="$(mktemp)"

  awk -F '\t' -v exclude_rows="${excluded_rows_path}" -v removed_file="${tmp_removed}" '
    BEGIN {
      while ((getline < exclude_rows) > 0) {
        split($0, fields, "\t")
        excluded[fields[1]] = 1
      }
    }
    excluded[$1] {
      removed += 1
      if (!seen_removed[$1]++) {
        print $1 >> removed_file
      }
      next
    }
    {
      print
    }
    END {
      printf "%d\n", removed + 0 > "/dev/stderr"
    }
  ' "${rows_path}" > "${tmp_rows}" 2> "${tmp_removed}.count"

  removed_count="$(tr -d '[:space:]' < "${tmp_removed}.count")"

  awk -v keep_rows="${tmp_rows}" '
    BEGIN {
      while ((getline < keep_rows) > 0) {
        split($0, fields, "\t")
        keep[fields[1]] = 1
      }
    }
    /^>/ {
      contig = substr($0, 2)
      split(contig, parts, /[[:space:]]+/)
      contig = parts[1]
      emit = keep[contig] && !written[contig]++
    }
    emit { print }
  ' "${fasta_path}" > "${tmp_fasta}"

  mv "${tmp_rows}" "${rows_path}"
  mv "${tmp_fasta}" "${fasta_path}"

  if [[ "${removed_count}" != "0" ]]; then
    removed_preview="$(paste -sd ', ' "${tmp_removed}" | cut -c1-180)"
    log_warn "removed ${removed_count} contig(s) from ${label} because they are handled by higher-priority layers: ${removed_preview}"
  fi

  rm -f "${tmp_removed}" "${tmp_removed}.count"
}

validate_unique_contigs() {
  local duplicates
  duplicates="$(cut -f1 "${HOST_ROWS}" "${VIRUS_ROWS}" "${DECOY_ROWS}" | sort | uniq -d | head -n 1 || true)"
  if [[ -n "${duplicates}" ]]; then
    fail "duplicate contig detected across organized inputs: ${duplicates}"
  fi
}

write_manifest_json() {
  awk -F '\t' '
    function escape_json(value, copy) {
      copy = value
      gsub(/\\/, "\\\\", copy)
      gsub(/\"/, "\\\"", copy)
      gsub(/\t/, " ", copy)
      return copy
    }
    BEGIN {
      print "["
    }
    {
      if (count > 0) {
        print ","
      }
      printf "  {\n"
      printf "    \"contig\": \"%s\",\n", escape_json($1)
      printf "    \"accession\": \"%s\",\n", escape_json($2)
      printf "    \"taxid\": %s,\n", $3
      printf "    \"virus_name\": \"%s\",\n", escape_json($4)
      printf "    \"segment\": null,\n"
      printf "    \"group\": \"%s\",\n", escape_json($5)
      printf "    \"genome_length\": %s,\n", $6
      printf "    \"source_release\": \"%s\",\n", escape_json($7)
      printf "    \"source_type\": \"%s\",\n", escape_json($8)
      printf "    \"masked_regions\": []\n"
      printf "  }"
      count += 1
    }
    END {
      if (count > 0) {
        printf "\n"
      }
      print "]"
    }
  ' "${HOST_ROWS}" "${VIRUS_ROWS}" "${DECOY_ROWS}" > "${MANIFEST_OUT}"
}

write_taxonomy_tsv() {
  {
    printf 'taxid\tname\n'
    awk -F '\t' '!seen[$3 FS $4]++ { printf "%s\t%s\n", $3, $4 }' "${HOST_ROWS}" "${VIRUS_ROWS}" "${DECOY_ROWS}"
  } > "${TAXONOMY_OUT}"
}

write_download_manifest() {
  cat > "${META_DIR}/download_manifest.tsv" <<EOF
dataset\tpath\turl
host_grch38_no_alt\t${HOST_DIR}/GRCh38_no_alt_analysis_set.fna.gz\t${GRCH38_FASTA_URL}
host_grch38_assembly_report\t${META_DIR}/GRCh38_assembly_report.txt\t${GRCH38_REPORT_URL}
decoy_univec_core\t${DECOY_DIR}/UniVec_Core.fa\t${UNIVEC_URL}
decoy_phix174\t${DECOY_DIR}/phix174.fa\t${PHIX_URL}
decoy_illumina_adapters\t${DECOY_DIR}/illumina_adapters.fa\tembedded_in_script
decoy_combined\t${DECOY_DIR}/test_decoy_combined.fa\tlocal_concat
bundle_input_host\t${HOST_FASTA_OUT}\tlocal_normalize
bundle_input_virus\t${VIRUS_FASTA_OUT}\tlocal_normalize
bundle_input_decoy\t${DECOY_FASTA_OUT}\tlocal_normalize
bundle_input_manifest\t${MANIFEST_OUT}\tlocal_generate
bundle_input_taxonomy\t${TAXONOMY_OUT}\tlocal_generate
EOF

  local url
  local file_name
  for url in "${viral_urls[@]}"; do
    file_name="$(basename "${url}")"
    printf 'viral_refseq\t%s\t%s\n' "${VIRAL_DIR}/${file_name}" "${url}" >> "${META_DIR}/download_manifest.tsv"
  done
}

organize_inputs() {
  local adapter_fasta="${DECOY_DIR}/illumina_adapters.fa"
  local combined_decoy="${DECOY_DIR}/test_decoy_combined.fa"
  local viral_sources=("${VIRAL_DIR}"/viral.*.genomic.fna.gz)

  write_adapter_fasta "${adapter_fasta}"
  cat "${DECOY_DIR}/phix174.fa" "${DECOY_DIR}/UniVec_Core.fa" "${adapter_fasta}" > "${combined_decoy}"
  log_ok "refreshed decoy source bundle: ${combined_decoy}"

  if needs_refresh "${HOST_FASTA_OUT}" "${HOST_DIR}/GRCh38_no_alt_analysis_set.fna.gz"; then
    log_info "normalize host FASTA -> ${HOST_FASTA_OUT}"
    normalize_host_fasta "${HOST_DIR}/GRCh38_no_alt_analysis_set.fna.gz" "${HOST_FASTA_OUT}" "${HOST_ROWS}"
    log_ok "wrote host bundle input"
  else
    log_ok "reuse organized host input: ${HOST_FASTA_OUT}"
  fi

  if needs_refresh "${VIRUS_FASTA_OUT}" "${viral_sources[@]}"; then
    log_info "normalize viral FASTA volumes -> ${VIRUS_FASTA_OUT}"
    normalize_viral_fasta "${VIRUS_FASTA_OUT}" "${VIRUS_ROWS}" "${viral_sources[@]}"
    dedupe_organized_layer "${VIRUS_FASTA_OUT}" "${VIRUS_ROWS}" "viral input"
    log_ok "wrote viral bundle input"
  else
    log_ok "reuse organized viral input: ${VIRUS_FASTA_OUT}"
  fi

  if needs_refresh "${DECOY_FASTA_OUT}" "${DECOY_DIR}/phix174.fa" "${DECOY_DIR}/UniVec_Core.fa" "${adapter_fasta}"; then
    log_info "normalize decoy FASTA sources -> ${DECOY_FASTA_OUT}"
    normalize_plain_fasta "${DECOY_DIR}/phix174.fa" "${DECOY_FASTA_OUT}" "${DECOY_ROWS}" 32630 decoy "${PHIX_SOURCE_RELEASE}" phix-decoy truncate
    normalize_plain_fasta "${DECOY_DIR}/UniVec_Core.fa" "${DECOY_FASTA_OUT}" "${DECOY_ROWS}" 32630 decoy "${UNIVEC_SOURCE_RELEASE}" vector-decoy append
    normalize_plain_fasta "${adapter_fasta}" "${DECOY_FASTA_OUT}" "${DECOY_ROWS}" 32630 decoy "${ADAPTER_SOURCE_RELEASE}" adapter-decoy append
    dedupe_organized_layer "${DECOY_FASTA_OUT}" "${DECOY_ROWS}" "decoy input"
    log_ok "wrote decoy bundle input"
  else
    log_ok "reuse organized decoy input: ${DECOY_FASTA_OUT}"
  fi

  exclude_contigs_from_layer "${VIRUS_FASTA_OUT}" "${VIRUS_ROWS}" "${DECOY_ROWS}" "viral input"

  validate_unique_contigs
  write_manifest_json
  write_taxonomy_tsv
  log_ok "wrote manifest: ${MANIFEST_OUT}"
  log_ok "wrote taxonomy: ${TAXONOMY_OUT}"
}

print_summary() {
  log_section "6/6" "summary"
  log_ok "host raw:    ${HOST_DIR}/GRCh38_no_alt_analysis_set.fna.gz"
  log_ok "viral raw:   ${VIRAL_DIR}"
  log_ok "decoy raw:   ${DECOY_DIR}/test_decoy_combined.fa"
  log_ok "host input:  ${HOST_FASTA_OUT}"
  log_ok "virus input: ${VIRUS_FASTA_OUT}"
  log_ok "decoy input: ${DECOY_FASTA_OUT}"
  log_ok "manifest:    ${MANIFEST_OUT}"
  log_ok "taxonomy:    ${TAXONOMY_OUT}"
  log_ok "download log: ${META_DIR}/download_manifest.tsv"
  printf '\n%b\n' "${BOLD}rvscreen ref build example${RESET}"
  printf '%s\n' "cargo run -- ref build --host-fasta \"${HOST_FASTA_OUT}\" --virus-fasta \"${VIRUS_FASTA_OUT}\" --decoy-fasta \"${DECOY_FASTA_OUT}\" --manifest \"${MANIFEST_OUT}\" --taxonomy \"${TAXONOMY_OUT}\" --out ./reference_bundle"
}

require_tool awk
require_tool cut
require_tool grep
require_tool gzip
require_tool sed
require_tool sort
require_tool uniq
require_tool wc

log_banner

log_section "1/6" "download host reference"
download "${GRCH38_FASTA_URL}" "${HOST_DIR}/GRCh38_no_alt_analysis_set.fna.gz"
download "${GRCH38_REPORT_URL}" "${META_DIR}/GRCh38_assembly_report.txt"

log_section "2/6" "discover viral reference volumes"
discover_viral_urls

log_section "3/6" "download viral reference volumes"
for url in "${viral_urls[@]}"; do
  file_name="$(basename "${url}")"
  download "${url}" "${VIRAL_DIR}/${file_name}"
done

log_section "4/6" "download and assemble decoy sources"
download "${UNIVEC_URL}" "${DECOY_DIR}/UniVec_Core.fa"
download_phix "${DECOY_DIR}/phix174.fa"

log_section "5/6" "organize bundle inputs"
organize_inputs
write_download_manifest

print_summary