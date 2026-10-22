#!/bin/bash
# =============================================================================
# steps/finalize.sh — copy the files needed downstream (DiffBind, browser, QC)
# from scratch ($OUT_ROOT/<sample>) to backed-up project storage ($PROJECT_RESULTS/<sample>)
# args: <pipeline_dir> <sample_dir>
# =============================================================================
set -euo pipefail
PIPE=$1; SD=$2
source "$PIPE/config.sh"; source "$SD/sample.env"
DST=$PROJECT_RESULTS/$SAMPLE
mkdir -p "$DST"/{bam,peaks,signal,qc,xcor}

for d in "$SD"/align/*/; do
  u=$(basename "$d")
  cp -p "$d/$u.nodup.bam" "$d/$u.nodup.bam.bai" "$DST/bam/"         # filtered, deduplicated BAMs (DiffBind input)
  cp -p "$d/$u.tagAlign.gz" "$DST/bam/"
  cp -p "$d"/*.qc "$d/$u.bowtie2.log" "$DST/qc/" 2>/dev/null || true
done
for d in "$SD"/xcor/*/; do
  r=$(basename "$d"); cp -p "$d/xcor.cc.qc" "$DST/xcor/$r.cc.qc"; cp -p "$d/xcor.cc.plot.pdf" "$DST/xcor/$r.cc.plot.pdf"
done
cp -p "$SD"/peaks/*.bfilt.narrowPeak.gz "$DST/peaks/"                  # per-rep, pseudo-rep, pooled peaks
cp -p "$SD"/overlap/*.overlap.bfilt.narrowPeak.gz "$DST/peaks/"
cp -p "$SD"/final/*.narrowPeak.gz "$DST/peaks/"                        # optimal / conservative (replicated) peaks
cp -p "$SD"/signal/*.bigwig "$DST/signal/"
cp -p "$SD"/qc/* "$DST/qc/"
cp -p "$SD"/sample.env "$SD"/units.tsv "$DST/"
echo "copied final outputs to $DST"
