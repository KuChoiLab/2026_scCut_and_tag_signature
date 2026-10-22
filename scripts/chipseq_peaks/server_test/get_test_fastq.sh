#!/bin/bash
# =============================================================================
# get_test_fastq.sh — download the test samples' FASTQs (ENCODE, md5 checked) and write their sample sheet
#
# Usage: bash server_test/get_test_fastq.sh <DATA_DIR> [SAMPLE[,SAMPLE]]
#   default: every sample in test_samples.tsv
#   -> <DATA_DIR>/<accession>.fastq.gz
#   -> <DATA_DIR>/samples_test.txt (input for chip_run_local.sh)
#
# Test samples (test_samples.tsv; all were run on the Yale cluster, verify.sh compares against them):
#   Stam_B_028NGF         H3K4me3 ChIP-seq, CD19+ B, donor 028NGF, single-end 36 bp, 4 files, 2.7 GB
#                         (also run with the ENCODE pipeline itself: caper/WDL v2.2.2)
#   Mint_naiveCD4_836ICQ  H3K4me3 MINT-ChIP, naive CD4+ T, donor 836ICQ, paired-end 76 bp, 40 files (20 lane pairs), 8.0 GB
# Lanes are listed in the order the Yale run concatenated them; the sheet keeps that order.
# Files already present with the right md5 are not downloaded again.
# =============================================================================
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
D=${1:?usage: bash get_test_fastq.sh <DATA_DIR> [SAMPLE[,SAMPLE]]}
T=$HERE/test_samples.tsv
ALL=$(awk -F'\t' '!/^#/ && !s[$1]++ {printf "%s%s", (n++ ? "," : ""), $1}' "$T")
WANT=${2:-$ALL}
mkdir -p "$D"; D=$(cd "$D" && pwd)
for s in ${WANT//,/ }; do
  awk -F'\t' -v s="$s" '$1==s{f=1} END{exit !f}' "$T" || { echo "unknown test sample '$s' (choose from: $ALL)" >&2; exit 1; }
done

get() {  # accession md5
  local f=$D/$1.fastq.gz
  if [[ -f $f && $(md5sum < "$f" | cut -c1-32) == "$2" ]]; then echo "ok (already there): $f"; return; fi
  echo "downloading $1"
  curl -fsSL --retry 5 -C - -o "$f.part" "https://www.encodeproject.org/files/$1/@@download/$1.fastq.gz"
  [[ $(md5sum < "$f.part" | cut -c1-32) == "$2" ]] || { echo "md5 mismatch: $f.part (delete it and retry)" >&2; exit 1; }
  mv "$f.part" "$f"; echo "ok: $f"
}

S=$D/samples_test.txt
printf "# sample_id\trole\tbio_rep(donor)\tfastq_R1\tfastq_R2\n" > "$S"
for s in ${WANT//,/ }; do
  while IFS=$'\t' read -r id role rep r1 m1 r2 m2; do
    get "$r1" "$m1"
    [[ $r2 == - ]] || get "$r2" "$m2"
    printf "%s\t%s\t%s\t%s\t%s\n" "$id" "$role" "$rep" "$D/$r1.fastq.gz" "$([[ $r2 == - ]] && echo - || echo "$D/$r2.fastq.gz")" >> "$S"
  done < <(awk -F'\t' -v s="$s" '$1==s' "$T")
done
echo "sample sheet: $S ($(grep -vc '^#' "$S") rows: ${WANT//,/ })"
