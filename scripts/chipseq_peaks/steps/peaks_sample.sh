#!/bin/bash
# =============================================================================
# steps/peaks_sample.sh — everything after alignment for one sample
# Runs in the tool environment (conda env or container).  args: <pipeline_dir> <sample_dir>
#
# ENCODE tasks re-implemented here (histone ChIP-seq, MACS2):
#   pool_ta -> fraglen_mean -> choose_ctl (+ subsample_ctl) -> call_peak (rep, pr1/pr2,
#   pooled, ppr1/ppr2) -> blacklist filter + FRiP -> macs2_signal_track -> naive overlap
#   (rep pairs, self pseudo-reps, pooled pseudo-reps) -> reproducibility (optimal/conservative)
# =============================================================================
PIPE=$1; SD=$2
source "$PIPE/config.sh"; source "$PIPE/lib.sh"; source "$SD/sample.env"
NTH=${NTH:-$PEAK_CPU}
N=$NUM_REP; M=$NUM_CTL
A=$SD/align; P=$SD/peaks; S=$SD/signal; O=$SD/overlap; F=$SD/final; C=$SD/ctl_choice
mkdir -p "$P" "$S" "$O" "$F" "$C"

ta()   { echo "$A/chip_rep$1/chip_rep$1.tagAlign.gz"; }
pr()   { echo "$A/chip_rep$1/chip_rep$1.pr$2.tagAlign.gz"; }
cta()  { echo "$A/ctl_rep$1/ctl_rep$1.tagAlign.gz"; }
fl()   { cat "$SD/xcor/rep$1/fraglen.txt"; }

# ---------------------------------------------------------------- pooling
# Pooled ChIP / pooled pseudo-reps / pooled control are only made when >1 replicate.
pool_all() {
  if (( N > 1 )); then
    pool_ta "$P/pooled.tagAlign.gz" $(for i in $(seq 1 $N); do ta $i; done)
    pool_ta "$P/ppr1.tagAlign.gz"   $(for i in $(seq 1 $N); do pr $i 1; done)
    pool_ta "$P/ppr2.tagAlign.gz"   $(for i in $(seq 1 $N); do pr $i 2; done)
  fi
  if (( M > 1 )); then
    pool_ta "$C/ctl_pooled.tagAlign.gz" $(for j in $(seq 1 $M); do cta $j; done)
  elif (( M == 1 )); then
    ln -sf "$(cta 1)" "$C/ctl_pooled.tagAlign.gz"      # single control = "pooled" control
  fi
  # mean fragment length (python2 round: half away from zero) -> used for pooled calls/overlaps
  for i in $(seq 1 $N); do fl $i; done | awk '{s+=$1} END{printf "%d\n", int(s/NR+0.5)}' > "$P/fraglen_mean.txt"
}

# ---------------------------------------------------------------- choose control
# Reproduces encode_task_choose_ctl.py; writes C/rep<i>.ctl.tagAlign.gz (+ C/pooled.ctl.tagAlign.gz)
choose_ctl() {
  local i j d dc lim idx
  declare -a DEP DCTL IDX
  for i in $(seq 1 $N); do DEP[$i]=$(nlines "$(ta $i)"); done
  for j in $(seq 1 $M); do DCTL[$j]=$(nlines "$(cta $j)"); done
  local dep_pool=0 dctl_pool=0
  for i in $(seq 1 $N); do dep_pool=$((dep_pool + DEP[i])); done
  for j in $(seq 1 $M); do dctl_pool=$((dctl_pool + DCTL[j])); done

  if (( M == 1 )); then
    for i in $(seq 1 $N); do IDX[$i]=1; done                     # one control for every replicate
  elif [[ $CTL_MODE == pooled ]]; then
    for i in $(seq 1 $N); do IDX[$i]=0; done                     # 0 = pooled control
  else
    local use_pooled=0
    for j in $(seq 1 $M); do for k in $(seq $((j+1)) $M); do     # controls differ too much in depth?
      awk -v a="${DCTL[j]}" -v b="${DCTL[k]}" -v r="$CTL_DEPTH_RATIO" 'BEGIN{exit !(a/b>r || b/a>r)}' && use_pooled=1
    done; done
    for i in $(seq 1 $N); do
      if (( use_pooled )) || (( i > M )) || (( DCTL[i] < DEP[i] )); then IDX[$i]=0; else IDX[$i]=$i; fi
    done
  fi

  : > "$C/chosen_ctl.tsv"
  for i in $(seq 1 $N); do
    if (( IDX[i] == 0 )); then src=$C/ctl_pooled.tagAlign.gz; dc=$dctl_pool; else src=$(cta ${IDX[i]}); dc=${DCTL[IDX[i]]}; fi
    lim=$(awk -v d="${DEP[i]}" -v r="$EXP_CTL_DEPTH_RATIO_LIMIT" -v l="$CTL_DEPTH_LIMIT" 'BEGIN{x=d*r; if (l>x) x=l; printf "%d", x}')
    if (( dc > lim )); then
      if [[ $CTL_PE == 1 ]]; then subsample_ta_pe "$src" "$lim" "$C/rep$i.ctl.tagAlign.gz"
      else subsample_ta_se "$src" "$lim" 0 "$C/rep$i.ctl.tagAlign.gz"; fi
    else
      ln -sf "$src" "$C/rep$i.ctl.tagAlign.gz"; lim=0
    fi
    printf "rep%d\tchip_reads=%d\tcontrol=%s\tcontrol_reads=%d\tsubsampled_to=%d\n" \
      "$i" "${DEP[i]}" "$([[ ${IDX[i]} == 0 ]] && echo pooled || echo ctl_rep${IDX[i]})" "$dc" "$lim" >> "$C/chosen_ctl.tsv"
  done
  if (( N > 1 )); then
    lim=$(awk -v d="$dep_pool" -v r="$EXP_CTL_DEPTH_RATIO_LIMIT" -v l="$CTL_DEPTH_LIMIT" 'BEGIN{x=d*r; if (l>x) x=l; printf "%d", x}')
    if (( dctl_pool > lim )); then
      if [[ $CTL_PE == 1 ]]; then subsample_ta_pe "$C/ctl_pooled.tagAlign.gz" "$lim" "$C/pooled.ctl.tagAlign.gz"
      else subsample_ta_se "$C/ctl_pooled.tagAlign.gz" "$lim" 0 "$C/pooled.ctl.tagAlign.gz"; fi
    else
      ln -sf "$C/ctl_pooled.tagAlign.gz" "$C/pooled.ctl.tagAlign.gz"; lim=0
    fi
    printf "pooled\tchip_reads=%d\tcontrol=pooled\tcontrol_reads=%d\tsubsampled_to=%d\n" "$dep_pool" "$dctl_pool" "$lim" >> "$C/chosen_ctl.tsv"
  fi
  cat "$C/chosen_ctl.tsv" >&2
}

# ---------------------------------------------------------------- MACS2 peak call (+ post-processing)
call_peak() {   # name chip_ta ctl_ta fraglen
  local name=$1 t=$2 c=$3 f=$4 pre=$P/$1
  macs2 callpeak -t "$t" -c "$c" -f BED -n "$pre" -g "$GENSZ" -p "$PVAL_THRESH" \
    --nomodel --shift 0 --extsize "$f" --keep-dup all -B --SPMR 2> "$pre.macs2.log"
  # rank by -log10 p (col 8), rename Peak_N, fix negative coords / missing summit, cap, clip
  # Missing summit ($10 == -1): written as start + half length, i.e. an absolute position, not the
  # narrowPeak offset. Kept exactly as ENCODE does it (chip-seq-pipeline2 v2.2.2 encode_task_macs2_chip.py)
  # so outputs stay identical. MACS2 callpeak always writes a valid offset, so the branch does not run:
  # 0 of ~1.9M peaks in the 12 Yale samples had a summit outside the peak (2026-10-02).
  sort -k 8gr,8gr -S "$(sort_mem)" "${pre}_peaks.narrowPeak" \
    | awk 'BEGIN{OFS="\t"}{$4="Peak_"NR; if ($2<0) $2=0; if ($3<0) $3=0; if ($10==-1) $10=$2+int(($3-$2+1)/2.0); print $0}' \
    | awk -v n="$CAP_NUM_PEAK" 'NR<=n' > "$pre.tmp"       # cap (awk reads to the end: no SIGPIPE under pipefail)
  bed_clip_gz "$pre.tmp" "$pre.narrowPeak.gz"
  rm -f "$pre.tmp" "${pre}_peaks.narrowPeak" "${pre}_peaks.xls" "${pre}_summits.bed" "${pre}_treat_pileup.bdg" "${pre}_control_lambda.bdg"
  blacklist_filter "$pre.narrowPeak.gz" "$pre.bfilt.narrowPeak.gz"
  [[ $(nlines "$pre.bfilt.narrowPeak.gz") -gt 0 ]] || die "$name: no peaks after blacklist filtering"
  frip_shifted "$t" "$pre.bfilt.narrowPeak.gz" "$f" "$pre.frip.qc"
  nlines "$pre.bfilt.narrowPeak.gz" > "$pre.num_peaks.qc"
}

# ---------------------------------------------------------------- MACS2 signal tracks (FE + -log10 p)
signal_track() {  # name chip_ta ctl_ta fraglen
  local name=$1 t=$2 c=$3 f=$4 pre=$S/$1 sval
  macs2 callpeak -t "$t" -c "$c" -f BED -n "$pre" -g "$GENSZ" -p "$PVAL_THRESH" \
    --nomodel --shift 0 --extsize "$f" --keep-dup all -B --SPMR 2> "$pre.macs2.log"
  sval=$(awk -v n="$(nlines "$t")" 'BEGIN{printf "%.6f", n/1000000.0}')   # full precision as ENCODE (awk print rounds to 6 sig. digits)
  for m in FE ppois; do
    if [[ $m == FE ]]; then macs2 bdgcmp -t "${pre}_treat_pileup.bdg" -c "${pre}_control_lambda.bdg" --o-prefix "$pre" -m FE
    else macs2 bdgcmp -t "${pre}_treat_pileup.bdg" -c "${pre}_control_lambda.bdg" --o-prefix "$pre" -m ppois -S "$sval"; fi
    local out; [[ $m == FE ]] && out=$pre.fc.signal || out=$pre.pval.signal
    bedtools slop -i "${pre}_${m}.bdg" -g "$CHRSZ" -b 0 | awk '{if ($3 != -1) print $0}' | bedClip stdin "$CHRSZ" "$out.bedgraph"
    sort -k1,1 -k2,2n -S "$(sort_mem)" "$out.bedgraph" \
      | awk 'BEGIN{OFS="\t"}{if (NR==1 || NR>1 && (prev_chr!=$1 || prev_chr==$1 && prev_chr_e<=$2)) {print $0}; prev_chr=$1; prev_chr_e=$3;}' \
      > "$out.srt.bedgraph"
    bedGraphToBigWig "$out.srt.bedgraph" "$CHRSZ" "$out.bigwig"
    rm -f "$out.bedgraph" "$out.srt.bedgraph" "${pre}_${m}.bdg"
  done
  rm -f "${pre}"_peaks.* "${pre}"_summits.bed "${pre}"_treat_pileup.bdg "${pre}"_control_lambda.bdg
}

# ---------------------------------------------------------------- naive overlap
# Keep pooled peaks that overlap a peak in BOTH inputs by >=50% of either peak's length.
naive_overlap() {  # name peak1 peak2 peak_pooled ta fraglen
  local name=$1 p1=$2 p2=$3 pp=$4 t=$5 f=$6 pre=$O/$1 tmpd
  tmpd=$(mktemp -d -p "$O")
  zcat -f "$p1" > "$tmpd/1"; zcat -f "$p2" > "$tmpd/2"; zcat -f "$pp" > "$tmpd/p"
  local aw='{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}'
  intersectBed -nonamecheck -wo -a "$tmpd/p" -b "$tmpd/1" | awk 'BEGIN{FS="\t";OFS="\t"} '"$aw" \
    | cut -f 1-10 | sort -S "$(sort_mem)" | uniq \
    | intersectBed -nonamecheck -wo -a stdin -b "$tmpd/2" | awk 'BEGIN{FS="\t";OFS="\t"} '"$aw" \
    | cut -f 1-10 | sort -S "$(sort_mem)" | uniq | gzip -nc > "$pre.overlap.narrowPeak.gz"
  rm -rf "$tmpd"
  blacklist_filter "$pre.overlap.narrowPeak.gz" "$pre.overlap.bfilt.narrowPeak.gz"
  frip_shifted "$t" "$pre.overlap.bfilt.narrowPeak.gz" "$f" "$pre.overlap.frip.qc"
  nlines "$pre.overlap.bfilt.narrowPeak.gz" > "$pre.overlap.num_peaks.qc"
}

# ---------------------------------------------------------------- run a list of jobs in parallel
run_parallel() {   # function_name  then lines "arg1 arg2 ..." on stdin
  # Every started job is waited for exactly once, by PID, oldest first, so no exit status is lost
  # (the earlier wait -n / jobs -rp version could miss a job that had already failed and exited).
  # After a failure no new job starts; the running ones are waited for, then the step fails.
  local fn=$1 max=$NTH line k=0 failed=0
  local -a pids=()
  while read -r line; do
    [[ -z $line ]] && continue
    (( failed )) && continue
    if (( ${#pids[@]} - k >= max )); then wait "${pids[k]}" || failed=1; k=$((k+1)); (( failed )) && continue; fi
    ( $fn $line ) & pids+=($!)
  done
  for (( ; k < ${#pids[@]}; k++ )); do wait "${pids[k]}" || failed=1; done
  (( failed == 0 )) || die "$fn: a job failed (see the messages above)"
}

peaks() {
  local FM; FM=$(cat "$P/fraglen_mean.txt")
  {
    for i in $(seq 1 $N); do
      echo "rep$i     $(ta $i)   $C/rep$i.ctl.tagAlign.gz $(fl $i)"
      echo "rep$i-pr1 $(pr $i 1) $C/rep$i.ctl.tagAlign.gz $(fl $i)"
      echo "rep$i-pr2 $(pr $i 2) $C/rep$i.ctl.tagAlign.gz $(fl $i)"
    done
    if (( N > 1 )); then
      echo "pooled $P/pooled.tagAlign.gz $C/pooled.ctl.tagAlign.gz $FM"
      echo "ppr1   $P/ppr1.tagAlign.gz   $C/pooled.ctl.tagAlign.gz $FM"
      echo "ppr2   $P/ppr2.tagAlign.gz   $C/pooled.ctl.tagAlign.gz $FM"
    fi
  } | run_parallel call_peak
}

signals() {
  local FM; FM=$(cat "$P/fraglen_mean.txt")
  {
    for i in $(seq 1 $N); do echo "rep$i $(ta $i) $C/rep$i.ctl.tagAlign.gz $(fl $i)"; done
    (( N > 1 )) && echo "pooled $P/pooled.tagAlign.gz $C/pooled.ctl.tagAlign.gz $FM"
    true
  } | run_parallel signal_track
}

overlaps() {
  local FM; FM=$(cat "$P/fraglen_mean.txt")
  {
    for i in $(seq 1 $N); do for j in $(seq $((i+1)) $N); do
      echo "rep${i}_vs_rep${j} $P/rep$i.narrowPeak.gz $P/rep$j.narrowPeak.gz $P/pooled.narrowPeak.gz $P/pooled.tagAlign.gz $FM"
    done; done
    for i in $(seq 1 $N); do
      echo "rep$i-pr1_vs_rep$i-pr2 $P/rep$i-pr1.narrowPeak.gz $P/rep$i-pr2.narrowPeak.gz $P/rep$i.narrowPeak.gz $(ta $i) $(fl $i)"
    done
    (( N > 1 )) && echo "pooled-pr1_vs_pooled-pr2 $P/ppr1.narrowPeak.gz $P/ppr2.narrowPeak.gz $P/pooled.narrowPeak.gz $P/pooled.tagAlign.gz $FM"
    true
  } | run_parallel naive_overlap
}

# ---------------------------------------------------------------- reproducibility (encode_task_reproducibility.py)
reproducibility() {
  local i j Nt=0 Nt_lab="" Np=0 n cons opt Nopt lab maxN=0 minN=0 rescue self repro
  declare -a NPR
  for i in $(seq 1 $N); do
    NPR[$i]=$(cat "$O/rep$i-pr1_vs_rep$i-pr2.overlap.num_peaks.qc")
    (( i == 1 || NPR[i] > maxN )) && maxN=${NPR[i]}
    (( i == 1 || NPR[i] < minN )) && minN=${NPR[i]}
  done
  if (( N > 1 )); then
    for i in $(seq 1 $N); do for j in $(seq $((i+1)) $N); do
      n=$(cat "$O/rep${i}_vs_rep${j}.overlap.num_peaks.qc")
      if (( n > Nt )); then Nt=$n; Nt_lab=rep${i}_vs_rep${j}; fi
    done; done
    Np=$(cat "$O/pooled-pr1_vs_pooled-pr2.overlap.num_peaks.qc")
    cons=$Nt_lab
    if (( Nt > Np )); then opt=$Nt_lab; Nopt=$Nt; else opt=pooled-pr1_vs_pooled-pr2; Nopt=$Np; fi
    rescue=$(awk -v a="$Nt" -v b="$Np" 'BEGIN{mx=(a>b)?a:b; mn=(a<b)?a:b; printf "%.6f", (mn>0)?mx/mn:0}')
  else
    cons=rep1-pr1_vs_rep1-pr2; opt=$cons; Nt=0; Np=0; Nopt=${NPR[1]}; rescue=0
  fi
  self=$(awk -v a="$maxN" -v b="$minN" 'BEGIN{printf "%.6f", (b>0)?a/b:0}')
  repro=pass
  awk -v r="$rescue" -v s="$self" 'BEGIN{exit !(r>2.0 || s>2.0)}' && repro=borderline
  awk -v r="$rescue" -v s="$self" 'BEGIN{exit !(r>2.0 && s>2.0)}' && repro=fail
  cp "$O/$opt.overlap.bfilt.narrowPeak.gz"  "$F/optimal_peak.narrowPeak.gz"
  cp "$O/$cons.overlap.bfilt.narrowPeak.gz" "$F/conservative_peak.narrowPeak.gz"
  {
    printf "Nt\tNp\tN1..N%d\tN_opt\tN_consv\topt_set\tconsv_set\trescue_ratio\tself_consistency_ratio\treproducibility\n" "$N"
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$Nt" "$Np" "$(IFS=,; echo "${NPR[*]:1}")" "$Nopt" \
      "$(nlines "$F/conservative_peak.narrowPeak.gz")" "$opt" "$cons" "$rescue" "$self" "$repro"
  } > "$F/reproducibility.qc"
  cat "$F/reproducibility.qc" >&2
}

step "$SD" pool_all        pool_all
step "$SD" choose_ctl      choose_ctl
step "$SD" peaks           peaks
step "$SD" signals         signals
step "$SD" overlaps        overlaps
step "$SD" reproducibility reproducibility
bash "$PIPE/steps/qc_summary.sh" "$PIPE" "$SD"
touch "$SD/.done.peaks_sample"
