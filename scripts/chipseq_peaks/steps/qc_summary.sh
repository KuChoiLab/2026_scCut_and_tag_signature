#!/bin/bash
# =============================================================================
# steps/qc_summary.sh — collect per-sample QC into two TSV tables
#   qc/align_qc.tsv : one row per replicate/control (reads, mapping, duplication, complexity, xcor)
#   qc/peak_qc.tsv  : one row per peak set (number of peaks, FRiP, peak width distribution)
#   qc/<set>.peak_region_size.{qc,pdf} : ENCODE peak_region_size for replicate peaks and final sets
# args: <pipeline_dir> <sample_dir>
# =============================================================================
PIPE=$1; SD=$2
source "$PIPE/config.sh"; source "$PIPE/lib.sh"; source "$SD/sample.env"
Q=$SD/qc; mkdir -p "$Q"

# peak width distribution for every blacklist-filtered peak set and the final optimal/conservative sets
for f in "$SD"/peaks/*.bfilt.narrowPeak.gz "$SD"/overlap/*.overlap.bfilt.narrowPeak.gz "$SD"/final/*_peak.narrowPeak.gz; do
  [[ -e $f ]] || continue
  n=$(basename "$f"); n=${n%.narrowPeak.gz}; n=${n%.bfilt}; n=${n%.overlap}
  [[ $f == */overlap/* ]] && n=overlap:$n
  [[ $f == */final/* ]] && n=final:$n
  Rscript "$PIPE/steps/peak_region_size.R" "$f" "$Q/${n//:/_}" || log "WARNING: peak_region_size failed for $n"   # QC only: never fail the peaks step
done
psize() { [[ -s $Q/${1//:/_}.peak_region_size.qc ]] && cut -f2 "$Q/${1//:/_}.peak_region_size.qc" | paste -sd'\t' || printf "NA\tNA\tNA\tNA\tNA\tNA"; }

{
  printf "sample\tunit\ttotal_reads\tmapped_reads\tmapped_pct\tfilt_nodup_reads\tdup_pct\tNRF\tPBC1\tPBC2\tfraglen\tNSC\tRSC\n"
  for d in "$SD"/align/*/; do
    u=$(basename "$d")
    tot=$(awk 'NR==1{print $1}' "$d/$u.flagstat.qc")
    map=$(grep -m1 " mapped (" "$d/$u.flagstat.qc" | awk '{print $1}')
    pct=$(awk -v a="$map" -v b="$tot" 'BEGIN{printf "%.2f", 100*a/b}')
    nod=$(awk 'NR==1{print $1}' "$d/$u.nodup.flagstat.qc")
    dup=$(awk -F'\t' '/^LIBRARY/{getline; print $9*100}' "$d/$u.dup.qc")   # PERCENT_DUPLICATION
    read -r _ _ _ _ nrf pbc1 pbc2 < "$d/$u.lib_complexity.qc"
    fl=NA; nsc=NA; rsc=NA
    if [[ $u == chip_* && -s $SD/xcor/${u#chip_}/xcor.cc.qc ]]; then
      read -r _ _ fl _ _ _ _ _ nsc rsc _ < "$SD/xcor/${u#chip_}/xcor.cc.qc"
    fi
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$SAMPLE" "$u" "$tot" "$map" "$pct" "$nod" "$dup" "$nrf" "$pbc1" "$pbc2" "$fl" "$nsc" "$rsc"
  done
} > "$Q/align_qc.tsv"

{
  printf "sample\tpeak_set\tnum_peaks\tFRiP\tsize_min\tsize_25pct\tsize_median\tsize_75pct\tsize_max\tsize_mean\n"
  for f in "$SD"/peaks/*.num_peaks.qc; do
    n=$(basename "$f" .num_peaks.qc)
    printf "%s\t%s\t%s\t%s\t%s\n" "$SAMPLE" "$n" "$(cat "$f")" "$(cat "$SD/peaks/$n.frip.qc")" "$(psize "$n")"
  done
  for f in "$SD"/overlap/*.overlap.num_peaks.qc; do
    n=$(basename "$f" .overlap.num_peaks.qc)
    printf "%s\toverlap:%s\t%s\t%s\t%s\n" "$SAMPLE" "$n" "$(cat "$f")" "$(cat "$SD/overlap/$n.overlap.frip.qc")" "$(psize "overlap:$n")"
  done
  for k in optimal conservative; do                       # final sets = copies of overlap sets (FRiP in that row)
    f=$SD/final/${k}_peak.narrowPeak.gz
    [[ -e $f ]] && printf "%s\tfinal:%s\t%s\tNA\t%s\n" "$SAMPLE" "$k" "$(nlines "$f")" "$(psize "final:${k}_peak")"
  done
} > "$Q/peak_qc.tsv"
cp "$SD/final/reproducibility.qc" "$Q/" 2>/dev/null || true
cp "$SD/ctl_choice/chosen_ctl.tsv" "$Q/" 2>/dev/null || true
log "QC tables written to $Q"
