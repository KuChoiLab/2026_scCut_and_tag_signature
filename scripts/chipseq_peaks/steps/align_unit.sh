#!/bin/bash
# =============================================================================
# steps/align_unit.sh — one replicate (ChIP or control): FASTQ -> filtered BAM -> tagAlign
# Runs in the tool environment (conda env or container; launched by steps/job.sh).
#   args: <pipeline_dir> <sample_dir> <unit_line_number>
# Unit line in <sample_dir>/units.tsv:  role  rep  paired(1/0)  R1_list  R2_list
#
# ENCODE tasks re-implemented here:
#   merge_fastq -> align (bowtie2) -> filter (samtools/picard) -> bam2ta -> spr (ChIP only)
# =============================================================================
PIPE=$1; SD=$2; LINE=$3
source "$PIPE/config.sh"; source "$PIPE/lib.sh"; source "$SD/sample.env"
NTH=${NTH:-$ALIGN_CPU}

read -r ROLE REP PE R1L R2L < <(sed -n "${LINE}p" "$SD/units.tsv")
U=${ROLE}_${REP}                      # e.g. chip_rep1 / ctl_rep2
D=$SD/align/$U; mkdir -p "$D"; cd "$D"
log "unit $U  paired=$PE  R1 files=$(wc -l < "$R1L")"

IDXF=("$BOWTIE2_IDX_DIR"/*.rev.1.bt2*); IDX=$(sed -E 's/\.rev\.1\.bt2l?$//' <<< "${IDXF[0]}")   # bowtie2 index prefix
JAVA_HEAP=${JAVA_HEAP:-8G}

# 1) merge technical replicates / lanes  (encode_task_merge_fastq.py: zcat -f | gzip -nc)
merge_fastq() {
  zcat -f $(cat "$R1L") | gzip -nc > R1.merged.fastq.gz
  if [[ $PE == 1 ]]; then zcat -f $(cat "$R2L") | gzip -nc > R2.merged.fastq.gz; fi
}

# 2) bowtie2 alignment  (encode_task_bowtie2.py)
align() {
  if [[ $PE == 1 ]]; then
    bowtie2 -X2000 --mm --threads "$NTH" -x "$IDX" -1 R1.merged.fastq.gz -2 R2.merged.fastq.gz \
      2> "$U.bowtie2.log" | samtools view -1 -S /dev/stdin > "$U.tmp.bam"
  else
    bowtie2 --mm --threads "$NTH" -x "$IDX" -U R1.merged.fastq.gz \
      2> "$U.bowtie2.log" | samtools view -1 -S /dev/stdin > "$U.tmp.bam"
  fi
  st_sort "$U.tmp.bam" -o "$U.bam" -T "$U.srt"
  rm -f "$U.tmp.bam"
  samtools index "$U.bam"
  samtools flagstat "$U.bam" > "$U.flagstat.qc"
  [[ $(samtools view -c "$U.bam") -gt 0 ]] || die "BAM is empty"
  rm -f R1.merged.fastq.gz R2.merged.fastq.gz
}

# 3) filter: unmapped / low MAPQ / not properly paired / secondary / QC-fail
#    -F 1804 = unmapped(4) + mate unmapped(8) + secondary(256) + QC fail(512) + duplicate(1024)
filter() {
  if [[ $PE == 1 ]]; then
    samtools view -F 1804 -f 2 -q "$MAPQ_THRESH" -u "$U.bam" | st_sort -n /dev/stdin -o "$U.tmp_filt.bam" -T "$U.nsrt"
    samtools fixmate -r "$U.tmp_filt.bam" "$U.fixmate.bam"           # -r: drop secondary/unmapped, fix mate info
    rm -f "$U.tmp_filt.bam"
    samtools view -F 1804 -f 2 -u "$U.fixmate.bam" | st_sort /dev/stdin -o "$U.filt.bam" -T "$U.fsrt"  # drop orphans
    rm -f "$U.fixmate.bam"
  else
    samtools view -F 1804 -q "$MAPQ_THRESH" -u "$U.bam" | st_sort /dev/stdin -o "$U.filt.bam" -T "$U.fsrt"
  fi
  [[ $(samtools view -c "$U.filt.bam") -gt 0 ]] || die "no reads left after filtering"
}

# 4) mark duplicates (Picard) and 5) remove them
dedup() {
  java -Xmx"$JAVA_HEAP" -XX:ParallelGCThreads=1 -jar "$PICARD_JAR" MarkDuplicates \
    INPUT="$U.filt.bam" OUTPUT="$U.dupmark.bam" METRICS_FILE="$U.dup.qc" \
    VALIDATION_STRINGENCY=LENIENT USE_JDK_DEFLATER=TRUE USE_JDK_INFLATER=TRUE \
    ASSUME_SORTED=TRUE REMOVE_DUPLICATES=FALSE
  if [[ $PE == 1 ]]; then
    samtools view -F 1804 -f 2 -b "$U.dupmark.bam" -@ "$((NTH-1))" > "$U.nodup.bam"
  else
    samtools view -F 1804 -b "$U.dupmark.bam" -@ "$((NTH-1))" > "$U.nodup.bam"
  fi
  if [[ -n $FILTER_CHRS ]]; then   # optional: remove chromosomes (e.g. chrM) from final BAM
    local keep; keep=$(awk -v re="^($(echo $FILTER_CHRS | tr ' ' '|'))$" '$1!~re{printf "%s:1-%s ",$1,$2}' "$CHRSZ")
    samtools index "$U.nodup.bam"
    samtools view -b "$U.nodup.bam" $keep | st_sort /dev/stdin -o "$U.nodup.tmp.bam" -T "$U.csrt"
    mv "$U.nodup.tmp.bam" "$U.nodup.bam"
  fi
  [[ $(samtools view -c "$U.nodup.bam") -gt 0 ]] || die "no reads left after dedup"
  samtools index "$U.nodup.bam"
  samtools flagstat "$U.nodup.bam" > "$U.nodup.flagstat.qc"
  rm -f "$U.filt.bam"
}

# 6) library complexity (PBC1, PBC2, NRF) from the duplicate-marked BAM
pbc_qc() {
  local calc='BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1}
    END{m1_m2=-1.0; if(m2>0) m1_m2=m1/m2; m0_mt=0; if (mt>0) m0_mt=m0/mt; m1_m0=0; if (m0>0) m1_m0=m1/m0;
    printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0_mt,m1_m0,m1_m2}'
  if [[ $PE == 1 ]]; then
    st_sort -n "$U.dupmark.bam" -o "$U.dupmark.nsrt.bam" -T "$U.psrt"
    bedtools bamtobed -bedpe -i "$U.dupmark.nsrt.bam" | awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$6,$9,$10}' \
      | { grep -v "^${MITO_CHR}\s" || true; } | sort -S "$(sort_mem)" | uniq -c | awk "$calc" > "$U.lib_complexity.qc"
    rm -f "$U.dupmark.nsrt.bam"
  else
    bedtools bamtobed -i "$U.dupmark.bam" | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6}' \
      | { grep -v "^${MITO_CHR}\s" || true; } | sort -S "$(sort_mem)" | uniq -c | awk "$calc" > "$U.lib_complexity.qc"
  fi
  # columns: TotalReadPairs DistinctReadPairs OneReadPair TwoReadPairs NRF=m0/mt PBC1=m1/m0 PBC2=m1/m2
  rm -f "$U.dupmark.bam"
}

# 7) BAM -> tagAlign (BED6: chr start end N 1000 strand); PE writes both mates
bam2ta() {
  if [[ $PE == 1 ]]; then
    st_sort -n "$U.nodup.bam" -o "$U.nodup.nsrt.bam" -T "$U.tsrt"
    bedtools bamtobed -bedpe -mate1 -i "$U.nodup.nsrt.bam" | gzip -nc > "$U.bedpe.gz"
    rm -f "$U.nodup.nsrt.bam"
    zcat -f "$U.bedpe.gz" | awk 'BEGIN{OFS="\t"}{printf "%s\t%s\t%s\tN\t1000\t%s\n%s\t%s\t%s\tN\t1000\t%s\n",$1,$2,$3,$9,$4,$5,$6,$10}' \
      | gzip -nc > "$U.tagAlign.gz"
    rm -f "$U.bedpe.gz"
  else
    bedtools bamtobed -i "$U.nodup.bam" | awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}' | gzip -nc > "$U.tagAlign.gz"
  fi
  [[ -s $U.tagAlign.gz ]] || die "empty tagAlign"
}

# 8) self pseudo-replicates (ChIP only): shuffle reads (pairs kept together) and split in half
spr() {
  local seed nl
  seed=$PSEUDOREP_SEED
  [[ $seed == 0 ]] && seed=$(zcat -f "$U.tagAlign.gz" | wc -c)
  if [[ $PE == 1 ]]; then
    nl=$(awk -v n="$(nlines "$U.tagAlign.gz")" 'BEGIN{print int((n/2+1)/2)}')
    zcat -f "$U.tagAlign.gz" | sed 'N;s/\n/\t/' | shuf --random-source=<(rand_src "$seed") | split -d -l "$nl" - "$U.spr."
    for k in 1 2; do
      zcat -f "$U.spr.0$((k-1))" | awk 'BEGIN{OFS="\t"}{printf "%s\t%s\t%s\t%s\t%s\t%s\n%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' \
        | gzip -nc > "$U.pr$k.tagAlign.gz"
    done
  else
    nl=$(awk -v n="$(nlines "$U.tagAlign.gz")" 'BEGIN{print int((n+1)/2)}')
    zcat -f "$U.tagAlign.gz" | shuf --random-source=<(rand_src "$seed") | split -d -l "$nl" - "$U.spr."
    gzip -nc "$U.spr.00" > "$U.pr1.tagAlign.gz"; gzip -nc "$U.spr.01" > "$U.pr2.tagAlign.gz"
  fi
  rm -f "$U.spr.0"*
}

step "$D" merge_fastq merge_fastq
step "$D" align       align
step "$D" filter      filter
step "$D" dedup       dedup
step "$D" pbc_qc      pbc_qc
step "$D" bam2ta      bam2ta
if [[ $ROLE == chip ]]; then step "$D" spr spr; fi
touch "$D/.done"
log "unit $U finished"
