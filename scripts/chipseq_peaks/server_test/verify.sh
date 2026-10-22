#!/bin/bash
# =============================================================================
# verify.sh — compare a new run of a test sample with the expected results
#
# Usage: CONDA_ENV=<ENV_PREFIX> bash server_test/verify.sh <OUTDIR>/results/<sample>
#        (<sample> = Stam_B_028NGF or Mint_naiveCD4_836ICQ: expected_<sample>.tsv)
#        (CONDA_ENV gives samtools for the BAM rows; without it those rows are skipped)
#
# expected_<sample>.tsv: item, how the file is hashed, path, expected md5 ('-' = not compared), source
#   gz     md5 of the decompressed content (gzip headers differ between machines)
#   raw    md5 of the file itself (bigwig)
#   nohash md5 without '#' lines (Picard writes the date and command line there)
#   bam    md5 of the reads (samtools view, no header: the header holds paths and command lines)
#   source yale   = Yale cluster run (2026-10-01), itself identical to the ENCODE pipeline for every file it shares
#          encode = ENCODE pipeline (caper/WDL v2.2.2) on the same FASTQs, used where the Yale file predates a fix
#          none   = no reference: the Yale file predates a fix and ENCODE was not run on this sample
#          info   = not compared on purpose, values printed only: xcor/rep1.cc.qc. Its column 4 (cross-correlation at
#                   the fragment length) changes in the last digits from run to run (multithreaded bowtie2 output order
#                   -> order of equal-coordinate reads -> which reads shuf puts in the 15M subsample), so the file md5
#                   is not reproducible. fraglen, NSC and RSC are still checked through qc/align_qc.tsv, and a different
#                   fraglen would change every peak file.
# Exit code 0 = every file identical.
# =============================================================================
set -uo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
R=${1:?usage: bash verify.sh <results>/Stam_B_028NGF [--write expected.tsv]}
WRITE=${3:-}; [[ ${2:-} == --write ]] || WRITE=""
SAMTOOLS=""
if [[ -n ${CONDA_ENV:-} && -x $CONDA_ENV/bin/samtools ]]; then SAMTOOLS=$CONDA_ENV/bin/samtools
elif command -v samtools >/dev/null; then SAMTOOLS=$(command -v samtools); fi

h() {  # mode file -> md5 | MISSING | SKIP
  local m=$1 f=$2
  [[ -f $f ]] || { echo MISSING; return; }
  case $m in
    gz)     zcat -f "$f" | md5sum | cut -c1-32 ;;
    raw)    md5sum < "$f" | cut -c1-32 ;;
    nohash) grep -v '^#' "$f" | md5sum | cut -c1-32 ;;
    bam)    [[ -n $SAMTOOLS ]] && "$SAMTOOLS" view "$f" | md5sum | cut -c1-32 || echo SKIP ;;
  esac
}

# ---- --write: make the expected table from a finished run (used once, on the Yale results)
if [[ -n $WRITE ]]; then
  { printf "item\tmode\tpath\tmd5\tsource\n"
    while IFS=$'\t' read -r item mode path; do
      if [[ $path == xcor/rep1.cc.qc ]]; then printf "%s\t%s\t%s\t-\tinfo\n" "$item" "$mode" "$path"; continue; fi
      printf "%s\t%s\t%s\t%s\tyale\n" "$item" "$mode" "$path" "$(h "$mode" "$R/$path")"
    done <<'EOF'
ChIP reads (tagAlign)	gz	bam/chip_rep1.tagAlign.gz
control reads (tagAlign)	gz	bam/ctl_rep1.tagAlign.gz
ChIP BAM reads	bam	bam/chip_rep1.nodup.bam
control BAM reads	bam	bam/ctl_rep1.nodup.bam
ChIP bowtie2 log	raw	qc/chip_rep1.bowtie2.log
control bowtie2 log	raw	qc/ctl_rep1.bowtie2.log
ChIP flagstat	raw	qc/chip_rep1.flagstat.qc
ChIP nodup flagstat	raw	qc/chip_rep1.nodup.flagstat.qc
ChIP duplicates (Picard)	nohash	qc/chip_rep1.dup.qc
ChIP library complexity	raw	qc/chip_rep1.lib_complexity.qc
control flagstat	raw	qc/ctl_rep1.flagstat.qc
control nodup flagstat	raw	qc/ctl_rep1.nodup.flagstat.qc
control duplicates (Picard)	nohash	qc/ctl_rep1.dup.qc
control library complexity	raw	qc/ctl_rep1.lib_complexity.qc
cross-correlation (fraglen, NSC, RSC)	raw	xcor/rep1.cc.qc
control choice	raw	qc/chosen_ctl.tsv
peaks rep1	gz	peaks/rep1.bfilt.narrowPeak.gz
peaks pseudo-rep 1	gz	peaks/rep1-pr1.bfilt.narrowPeak.gz
peaks pseudo-rep 2	gz	peaks/rep1-pr2.bfilt.narrowPeak.gz
peaks pseudo-rep overlap	gz	peaks/rep1-pr1_vs_rep1-pr2.overlap.bfilt.narrowPeak.gz
peaks optimal	gz	peaks/optimal_peak.narrowPeak.gz
peaks conservative	gz	peaks/conservative_peak.narrowPeak.gz
reproducibility	raw	qc/reproducibility.qc
alignment QC table	raw	qc/align_qc.tsv
peak QC table (FRiP, sizes)	raw	qc/peak_qc.tsv
signal fold enrichment	raw	signal/rep1.fc.signal.bigwig
signal -log10 p	raw	signal/rep1.pval.signal.bigwig
EOF
  } > "$WRITE"
  echo "wrote $WRITE"; exit 0
fi

# ---- compare
R=${R%/}; E=$HERE/expected_$(basename "$R").tsv
[[ -d $R ]] || { echo "no such directory: $R" >&2; exit 1; }
[[ -f $E ]] || { echo "no expected table for $(basename "$R"): $E" >&2; exit 1; }
echo "comparing $R"
echo "with      $E${SAMTOOLS:+   (samtools: $SAMTOOLS)}"
n=0; nd=0; ns=0; nn=0
while IFS=$'\t' read -r item mode path md5 src; do
  [[ $item == item ]] && continue
  n=$((n+1))
  if [[ $md5 == - && $src == info ]]; then
    v=$(awk -F'\t' '{printf "fraglen %s, NSC %s, RSC %s", $3, $9, $10}' "$R/$path" 2>/dev/null)
    printf "  %-4s  %-38s %s\n" "info" "$item" "$path (not compared: varies in the last digits between runs;"
    printf "  %-4s  %-38s %s\n" "" "" " fraglen/NSC/RSC are checked in qc/align_qc.tsv. This run: ${v:-file missing})"
    nn=$((nn+1)); continue
  fi
  if [[ $md5 == - ]]; then printf "  %-4s  %-38s %s\n" "n/a" "$item" "$path (no reference: $src)"; nn=$((nn+1)); continue; fi
  got=$(h "$mode" "$R/$path")
  if [[ $got == SKIP ]]; then st=SKIP; ns=$((ns+1))
  elif [[ $got == "$md5" ]]; then st=SAME
  else st=DIFF; nd=$((nd+1)); fi
  printf "  %-4s  %-38s %s\n" "$st" "$item" "$path"
done < "$E"
echo
echo "$(basename "$R"): $((n-nd-ns-nn)) same, $nd different, $ns skipped, $nn not compared (of $n)"
if (( nd == 0 && ns == 0 )); then echo "PASS: identical to the Yale/ENCODE results"
elif (( nd == 0 )); then echo "PASS (BAM rows skipped: set CONDA_ENV to check them)"
else echo "FAIL: send the DIFF lines and qc/align_qc.tsv, qc/peak_qc.tsv, xcor/rep1.cc.qc of this run"; fi
(( nd == 0 ))
