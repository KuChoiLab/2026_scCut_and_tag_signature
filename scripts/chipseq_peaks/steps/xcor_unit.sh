#!/bin/bash
# =============================================================================
# steps/xcor_unit.sh — cross-correlation (phantompeakqualtools) for one ChIP replicate
# Runs in the tool environment (conda env or container).  args: <pipeline_dir> <sample_dir> <unit_line_number>
#
# ENCODE logic (chip.wdl: align_R1 -> filter_R1 -> bam2ta_no_dedup_R1 -> xcor):
#   R1 only, trimmed to 50 bp, aligned single-end, MAPQ-filtered but NOT deduplicated,
#   15M non-mito reads subsampled, run_spp.R -> NSC, RSC and the fragment length
#   used later as MACS2 --extsize.
# =============================================================================
PIPE=$1; SD=$2; LINE=$3
source "$PIPE/config.sh"; source "$PIPE/lib.sh"; source "$SD/sample.env"
NTH=${NTH:-$XCOR_CPU}

read -r ROLE REP PE R1L R2L < <(sed -n "${LINE}p" "$SD/units.tsv")
[[ $ROLE == chip ]] || die "xcor runs on ChIP units only"
D=$SD/xcor/$REP; mkdir -p "$D"; cd "$D"
IDXF=("$BOWTIE2_IDX_DIR"/*.rev.1.bt2*); IDX=$(sed -E 's/\.rev\.1\.bt2l?$//' <<< "${IDXF[0]}")   # bowtie2 index prefix

# 1) merge R1 files and trim reads to XCOR_TRIM_BP (trimfastq.py, same script ENCODE uses)
trim_r1() {
  zcat -f $(cat "$R1L") | gzip -nc > R1.merged.fastq.gz
  python "$(command -v trimfastq.py)" R1.merged.fastq.gz "$XCOR_TRIM_BP" | gzip -nc > R1.trim.fastq.gz
  # reads shorter than XCOR_TRIM_BP: trimfastq.py only prints a message -> use untrimmed R1 (ENCODE behaviour)
  if (( $(zcat -f R1.trim.fastq.gz | { grep -c 'sequences shorter than desired length' || true; }) > 0 )); then
    log "reads shorter than ${XCOR_TRIM_BP} bp: using untrimmed R1"
    mv R1.merged.fastq.gz R1.trim.fastq.gz
  fi
  rm -f R1.merged.fastq.gz
}

# 2) single-end bowtie2 alignment of trimmed R1
align_r1() {
  bowtie2 --mm --threads "$NTH" -x "$IDX" -U R1.trim.fastq.gz 2> R1.bowtie2.log \
    | samtools view -1 -S /dev/stdin > R1.tmp.bam
  st_sort R1.tmp.bam -o R1.bam -T R1.srt
  rm -f R1.tmp.bam R1.trim.fastq.gz
}

# 3) filter (no duplicate removal) and convert to tagAlign
filter_ta() {
  samtools view -F 1804 -q "$MAPQ_THRESH" -u R1.bam | st_sort /dev/stdin -o R1.filt.bam -T R1.fsrt
  bedtools bamtobed -i R1.filt.bam | awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}' | gzip -nc > R1.tagAlign.gz
  rm -f R1.bam R1.filt.bam
}

# 4) subsample 15M non-mito reads and run SPP cross-correlation
xcor() {
  subsample_ta_se R1.tagAlign.gz "$XCOR_SUBSAMPLE_READS" 1 R1.sub.tagAlign.gz
  local rl xmax
  rl=$(zcat -f R1.sub.tagAlign.gz | awk 'function abs(v){return v<0?-v:v} NR<=100{s+=abs($3-$2); n=NR} END{print int(s/n)}')  # mean length of first 100 reads
  xmax=$(( rl + 10 > 100 ? rl + 10 : 100 ))          # histone: max(read_len+10, 100)
  Rscript --max-ppsize=500000 "$(command -v run_spp.R)" -rf -c=R1.sub.tagAlign.gz -p="$NTH" \
    -filtchr="$MITO_CHR" -savp=xcor.cc.plot.pdf -out=xcor.cc.qc -x="${XCOR_EXCLUSION_RANGE_MIN}:${xmax}"
  sed -r 's/,[^\t]+//g' -i xcor.cc.qc                  # keep the top peak only
  # cc.qc columns: file numReads estFragLen corr_estFragLen phantomPeak corr_phantomPeak argmin_corr min_corr NSC RSC QualityTag
  awk '{print $3}' xcor.cc.qc > fraglen.txt
  rm -f R1.sub.tagAlign.gz R1.tagAlign.gz
  log "fraglen $(cat fraglen.txt)  NSC $(awk '{print $9}' xcor.cc.qc)  RSC $(awk '{print $10}' xcor.cc.qc)"
}

step "$D" trim_r1   trim_r1
step "$D" align_r1  align_r1
step "$D" filter_ta filter_ta
step "$D" xcor      xcor
touch "$D/.done"
