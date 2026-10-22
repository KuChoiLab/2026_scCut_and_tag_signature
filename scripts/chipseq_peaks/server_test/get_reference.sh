#!/bin/bash
# =============================================================================
# get_reference.sh — download the hg38 files the pipeline uses (ENCODE portal) and check md5
#
# Usage: bash server_test/get_reference.sh <REF_DIR>
#   then run the pipeline with REF=<REF_DIR>
#
# Files (same as reference/hg38 on the Yale cluster, same as ENCODE chip-seq-pipeline2 v2.2.2 hg38):
#   <REF_DIR>/bowtie2_index/ENCFF110MCL.tar.gz   bowtie2 index, GRCh38 no-alt analysis set (3.5 GB)
#   <REF_DIR>/GRCh38_EBV.chrom.sizes.tsv         chromosome sizes
#   <REF_DIR>/ENCFF356LFX.bed.gz                 ENCODE blacklist v2
# The bowtie2 index is extracted by chip_run.sh on the first run (<REF_DIR>/bowtie2_index/ENCFF110MCL_untar, ~4 GB).
# Re-running skips files that are already complete (md5 OK).
# =============================================================================
set -euo pipefail
REF=${1:?usage: bash get_reference.sh <REF_DIR>}
mkdir -p "$REF/bowtie2_index"; REF=$(cd "$REF" && pwd)
E=https://www.encodeproject.org/files

get() {  # url dest md5
  local url=$1 dst=$2 md5=$3
  if [[ -f $dst && $(md5sum < "$dst" | cut -c1-32) == "$md5" ]]; then echo "ok (already there): $dst"; return; fi
  echo "downloading $url"
  curl -fL --retry 5 -C - -o "$dst.part" "$url"
  [[ $(md5sum < "$dst.part" | cut -c1-32) == "$md5" ]] || { echo "md5 mismatch: $dst.part (delete it and retry)" >&2; exit 1; }
  mv "$dst.part" "$dst"; echo "ok: $dst"
}
get "$E/ENCFF110MCL/@@download/ENCFF110MCL.tar.gz"                  "$REF/bowtie2_index/ENCFF110MCL.tar.gz" f23c2cabf3a91cc3ca953969ce979bec
get "$E/GRCh38_EBV.chrom.sizes/@@download/GRCh38_EBV.chrom.sizes.tsv" "$REF/GRCh38_EBV.chrom.sizes.tsv"       c95303fb77cc3e11d50e3c3a4b93b3fb
get "$E/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz"                  "$REF/ENCFF356LFX.bed.gz"                 393688b4f06c9ce26165d47433dd8c37
echo "reference ready: REF=$REF"
