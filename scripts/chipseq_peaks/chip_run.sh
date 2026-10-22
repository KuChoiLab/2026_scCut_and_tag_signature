#!/bin/bash
# =============================================================================
# chip_run.sh — ENCODE histone ChIP-seq pipeline logic in bash + Slurm (no WDL/caper)
#
# Usage:
#   bash chip_run.sh <samples.txt> [--only ID[,ID]] [--ctl-mode pooled|matched] [--runtime conda|sif] [--max-parallel N] [--out DIR] [--results DIR]
#                     [--executor slurm|local] [--threads N] [--dry-run]
#   (chip_run_local.sh is the short form for --executor local)
#
# samples.txt: tab-separated, '#' lines ignored
#   sample_id  role(chip|ctl)  bio_rep  fastq_R1  fastq_R2('-' = single-end)
#
#   * rows with the same (sample_id, role, bio_rep) = technical replicates / lanes
#     -> merged into one replicate (FASTQs concatenated in row order)
#   * distinct bio_rep labels = biological replicates -> rep1..repN (sorted labels)
#   * control bio_rep labels must equal the ChIP labels, or give one control only
#
# --executor slurm (default) submits 4 dependent Slurm jobs per sample:
#   align (array: every ChIP + control replicate) ─┐
#   xcor  (array: every ChIP replicate)          ──┴─> peaks ──> finalize
# --executor local runs the same 4 steps one after another in this shell, each step with --threads CPUs;
# --max-parallel N samples run at the same time (local default 1 = one sample after another,
# total CPUs = N x --threads). A failed sample is reported and the other samples still run.
# At the end the QC table of all samples is written: <results>/QC_summary_<sheet name>.tsv (qc_aggregate.py).
# Slurm: at most --max-parallel samples (default MAX_PARALLEL in config.sh) run at once: sample k+N
# starts only after sample k has finished (or failed), so one failure never blocks the queue.
# Re-running the same command resumes: finished steps (.done.* markers) are skipped.
# =============================================================================
set -euo pipefail
export LC_ALL=C
PIPE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# absolute paths: steps run from other directories and record these in sample.env
[[ -n ${REF:-} && -d $REF ]] && REF=$(cd "$REF" && pwd)
[[ -n ${CONDA_ENV:-} && -d $CONDA_ENV ]] && CONDA_ENV=$(cd "$CONDA_ENV" && pwd)
source "$PIPE/config.sh"

SAMPLES=${1:?usage: bash chip_run.sh <samples.txt> [--only ID] [--ctl-mode pooled|matched] [--runtime conda|sif] [--executor slurm|local] [--threads N] [--max-parallel N] [--out DIR] [--results DIR] [--dry-run]}
shift
ONLY=""; DRY=0; MODE=$CTL_MODE; MAXP=""; EXECUTOR=slurm; THREADS=${THREADS:-8}
while [[ $# -gt 0 ]]; do
  case $1 in
    --only) ONLY=$2; shift ;;
    --ctl-mode) MODE=$2; shift ;;
    --runtime) RUNTIME=$2; shift ;;
    --max-parallel) MAXP=$2; shift ;;
    --out) OUT_ROOT=$2; shift ;;
    --results) PROJECT_RESULTS=$2; shift ;;
    --executor) EXECUTOR=$2; shift ;;
    --threads) THREADS=$2; shift ;;
    --dry-run) DRY=1 ;;
    *) echo "unknown option $1" >&2; exit 1 ;;
  esac; shift
done
[[ $MODE == pooled || $MODE == matched ]] || { echo "--ctl-mode must be pooled|matched" >&2; exit 1; }
[[ -n $MAXP ]] || { [[ $EXECUTOR == local ]] && MAXP=1 || MAXP=$MAX_PARALLEL; }   # local: one sample at a time unless asked
[[ $MAXP =~ ^[1-9][0-9]*$ ]] || { echo "--max-parallel must be a positive integer" >&2; exit 1; }
[[ $THREADS =~ ^[1-9][0-9]*$ ]] || { echo "--threads must be a positive integer" >&2; exit 1; }
[[ $EXECUTOR == slurm || $EXECUTOR == local ]] || { echo "--executor must be slurm|local" >&2; exit 1; }
[[ $EXECUTOR == local ]] || command -v sbatch >/dev/null || { echo "sbatch not found: use --executor local (chip_run_local.sh)" >&2; exit 1; }
mkdir -p "$OUT_ROOT" "$PROJECT_RESULTS"
OUT_ROOT=$(cd "$OUT_ROOT" && pwd); PROJECT_RESULTS=$(cd "$PROJECT_RESULTS" && pwd)   # absolute: steps run from other dirs
case $RUNTIME in
  conda) [[ -x $CONDA_ENV/bin/macs2 && -x $CONDA_ENV/bin/run_spp.R && -f $CONDA_ENV/share/picard-2.20.7-0/picard.jar ]] \
           || { echo "conda env incomplete: $CONDA_ENV (see env/chipseq_v2.2.2.yml)" >&2; exit 1; } ;;
  sif)   [[ -f $SIF ]] || { echo "missing $SIF" >&2; exit 1; } ;;
  *)     echo "--runtime must be conda|sif" >&2; exit 1 ;;
esac

# ---- one-time: extract the bowtie2 index tarball (into a temporary folder, renamed when complete,
#      so an interrupted extraction is redone instead of being taken for a finished one)
if ! bt2_index_ok; then
  [[ -f $BOWTIE2_IDX_TAR ]] || { echo "missing $BOWTIE2_IDX_TAR (REF=$REF; see server_test/get_reference.sh)" >&2; exit 1; }
  echo "extracting bowtie2 index -> $BOWTIE2_IDX_DIR"
  rm -rf "$BOWTIE2_IDX_DIR.tmp"; mkdir -p "$BOWTIE2_IDX_DIR.tmp"
  tar -xzf "$BOWTIE2_IDX_TAR" -C "$BOWTIE2_IDX_DIR.tmp"
  rm -rf "$BOWTIE2_IDX_DIR"; mv "$BOWTIE2_IDX_DIR.tmp" "$BOWTIE2_IDX_DIR"
  bt2_index_ok || { echo "bowtie2 index incomplete after extraction: $BOWTIE2_IDX_DIR" >&2; exit 1; }
fi

# ---- sample ids in file order
mapfile -t IDS < <(awk -F'\t' '!/^#/ && NF>=5 && $1!="sample_id" && !seen[$1]++ {print $1}' "$SAMPLES")
if [[ -n $ONLY ]]; then
  IFS=, read -ra keep <<< "$ONLY"
  for k in "${keep[@]}"; do printf "%s\n" "${IDS[@]}" | grep -Fxq -- "$k" || { echo "--only: sample '$k' is not in $SAMPLES" >&2; exit 1; }; done
  mapfile -t IDS < <(printf "%s\n" "${IDS[@]}" | grep -Fxf <(printf "%s\n" "${keep[@]}"))
fi
[[ ${#IDS[@]} -gt 0 ]] || { echo "no samples to run" >&2; exit 1; }

# ---- check the whole sheet before anything runs: one bad row stops the run, nothing is half done
NAME_RE='^[A-Za-z0-9][A-Za-z0-9._-]*$'      # sample_id, bio_rep: letters, digits, . _ - (no spaces, no /)
PATH_RE='^/[A-Za-z0-9._/+=,@:-]+$'          # FASTQ: absolute path, no spaces or shell characters
sheet_bad=0
while IFS= read -r ln; do   # every data line must have exactly 5 tab-separated columns
  echo "sheet line ${ln%%:*}: needs 5 tab-separated columns (found ${ln#*:}); spaces instead of tabs?" >&2; sheet_bad=1
done < <(awk -F'\t' '!/^#/ && NF>0 && $1!="sample_id" && NF!=5 {print NR":"NF}' "$SAMPLES")
check_sample() {   # ID -> prints every problem, returns 1 if any
  local ID=$1 rows bad=0 role lab r1 r2 cl tl n m pe_l prev
  [[ $ID =~ $NAME_RE ]] || { echo "[$ID] sample_id may only use letters, digits, . _ - (and must start with a letter or digit)" >&2; return 1; }
  rows=$(awk -F'\t' -v id="$ID" '!/^#/ && $1==id' "$SAMPLES")
  while IFS=$'\t' read -r _ role lab r1 r2; do
    [[ $role == chip || $role == ctl ]] || { echo "[$ID] role must be chip or ctl, not '$role'" >&2; bad=1; }
    [[ $lab =~ $NAME_RE ]] || { echo "[$ID] bio_rep '$lab' may only use letters, digits, . _ -" >&2; bad=1; }
    for f in "$r1" "$r2"; do
      [[ $f == - && $f == "$r2" ]] && continue
      if [[ ! $f =~ $PATH_RE ]]; then echo "[$ID] FASTQ path must be absolute, without spaces or special characters: '$f'" >&2; bad=1
      elif [[ ! -s $f ]]; then echo "[$ID] missing or empty FASTQ: $f" >&2; bad=1; fi
    done
  done <<< "$rows"
  cl=$(awk -F'\t' '$2=="chip"{print $3}' <<< "$rows" | sort -u | tr '\n' ' ')
  tl=$(awk -F'\t' '$2=="ctl"{print $3}'  <<< "$rows" | sort -u | tr '\n' ' ')
  n=$(wc -w <<< "$cl"); m=$(wc -w <<< "$tl")
  (( n >= 1 )) || { echo "[$ID] no chip rows" >&2; bad=1; }
  (( m >= 1 )) || { echo "[$ID] no control (ctl) rows: a control is required" >&2; bad=1; }
  (( n <= 10 && m <= 10 )) || { echo "[$ID] max 10 replicates" >&2; bad=1; }
  if (( m > 1 )) && [[ $cl != "$tl" ]]; then
    echo "[$ID] control bio_rep labels ($tl) must match the chip labels ($cl), or use one control" >&2; bad=1
  fi
  # a replicate is all paired-end or all single-end; all controls too
  while read -r role lab; do echo "[$ID] $role $lab mixes paired-end and single-end rows" >&2; bad=1
  done < <(awk -F'\t' '{k=$2" "$3; pe=($5=="-")?"SE":"PE"; if(k in t && t[k]!=pe) bad[k]=1; t[k]=pe} END{for(k in bad) print k}' <<< "$rows")
  pe_l=$(awk -F'\t' '$2=="ctl"{print ($5=="-")?"SE":"PE"}' <<< "$rows" | sort -u | wc -l)
  (( pe_l <= 1 )) || { echo "[$ID] controls mix paired-end and single-end" >&2; bad=1; }
  # the reference recorded by an earlier run of this sample must not change (steps would mix references)
  prev=$( (source "$OUT_ROOT/$ID/sample.env" 2>/dev/null && echo "${REF:-}") || true)
  if [[ -n $prev && $prev != "$REF" ]]; then
    echo "[$ID] this sample was started with REF=$prev, now REF=$REF: use the same REF, or delete $OUT_ROOT/$ID to start over" >&2; bad=1
  fi
  return $bad
}
BAD_IDS=()
for ID in "${IDS[@]}"; do check_sample "$ID" || BAD_IDS+=("$ID"); done
if (( sheet_bad || ${#BAD_IDS[@]} )); then
  echo "sample sheet has problems${BAD_IDS[*]:+ in ${#BAD_IDS[@]} sample(s): ${BAD_IDS[*]}}; nothing was run. Fix $SAMPLES and run again." >&2
  exit 1
fi
LOCAL_RUN=(); LOCAL_PIDS=()   # local executor: samples started, background pids (--max-parallel > 1)
STATUS_DIR=$OUT_ROOT/.local_status   # local executor: <sample> -> ok | failed
# local executor: run align -> xcor -> peaks -> finalize for one sample, one step at a time
local_sample() {   # ID SD CODE
  local ID=$1 SD=$2 CODE=$3 na nx ok=1 t
  na=$(wc -l < "$SD/align_units.txt"); nx=$(wc -l < "$SD/xcor_units.txt")
  run() {   # step [task] -> logs/<step>[_task].log
    local lg=$SD/logs/$1${2:+_$2}.log
    echo "[$(date '+%F %T')] [$ID] $1${2:+ $2} ... (log: $lg)"
    NTH=$THREADS bash "$CODE/steps/job.sh" "$CODE" "$1" "$SD" ${2:+"$2"} >> "$lg" 2>&1 \
      || { echo "[$(date '+%F %T')] [$ID] FAILED at $1${2:+ $2}: see $lg" >&2; return 1; }
  }
  for t in $(seq 1 "$na"); do run align "$t" || { ok=0; break; }; done
  if (( ok )); then for t in $(seq 1 "$nx"); do run xcor "$t" || { ok=0; break; }; done; fi
  (( ok )) && run peaks || ok=0
  (( ok )) && run finalize || ok=0
  if (( ok )); then echo ok > "$STATUS_DIR/$ID"; echo "[$(date '+%F %T')] [$ID] done -> $PROJECT_RESULTS/$ID"
  else echo failed > "$STATUS_DIR/$ID"; fi
  return 0
}
n_running() { local p n=0; for p in "${LOCAL_PIDS[@]}"; do kill -0 "$p" 2>/dev/null && n=$((n+1)); done; echo "$n"; }
[[ $EXECUTOR == local ]] && mkdir -p "$STATUS_DIR"
DONE_JOBS=()   # finalize job id of every submitted sample, in submission order (parallel-slot chain)

for ID in "${IDS[@]}"; do
  SD=$OUT_ROOT/$ID
  mkdir -p "$SD/inputs" "$SD/logs"
  rows=$(awk -F'\t' -v id="$ID" '!/^#/ && $1==id' "$SAMPLES")

  # (the sheet was checked above, before anything ran)
  mapfile -t CHIP_LABS < <(awk -F'\t' '$2=="chip"{print $3}' <<< "$rows" | sort -u)
  mapfile -t CTL_LABS  < <(awk -F'\t' '$2=="ctl"{print $3}'  <<< "$rows" | sort -u)
  N=${#CHIP_LABS[@]}; M=${#CTL_LABS[@]}

  # ---- per-replicate FASTQ lists + unit table
  : > "$SD/units.tsv"; : > "$SD/align_units.txt"; : > "$SD/xcor_units.txt"
  line=0; declare -A PE_OF=()
  for role in chip ctl; do
    if [[ $role == chip ]]; then labs=("${CHIP_LABS[@]}"); else labs=("${CTL_LABS[@]}"); fi
    k=0
    for lab in "${labs[@]}"; do
      k=$((k+1)); rep=rep$k
      awk -F'\t' -v r="$role" -v l="$lab" '$2==r && $3==l {print $4}' <<< "$rows" > "$SD/inputs/${role}_$rep.R1.txt"
      awk -F'\t' -v r="$role" -v l="$lab" '$2==r && $3==l {print $5}' <<< "$rows" > "$SD/inputs/${role}_$rep.R2.txt"
      nse=$(grep -cx -- '-' "$SD/inputs/${role}_$rep.R2.txt" || true); ntot=$(wc -l < "$SD/inputs/${role}_$rep.R2.txt")
      if (( nse == 0 )); then pe=1; elif (( nse == ntot )); then pe=0; else echo "[$ID] $role $lab mixes PE and SE rows" >&2; exit 1; fi
      PE_OF[$role]+="$pe"
      line=$((line+1))
      printf "%s\t%s\t%s\t%s\t%s\n" "$role" "$rep" "$pe" "$SD/inputs/${role}_$rep.R1.txt" "$SD/inputs/${role}_$rep.R2.txt" >> "$SD/units.tsv"
      echo "$line" >> "$SD/align_units.txt"
      [[ $role == chip ]] && echo "$line" >> "$SD/xcor_units.txt"
      printf "  %-4s %-5s <- bio_rep %-14s %3d FASTQ unit(s)  %s\n" "$role" "$rep" "$lab" "$ntot" "$([[ $pe == 1 ]] && echo PE || echo SE)"
    done
  done
  [[ ${PE_OF[ctl]} =~ ^(1+|0+)$ ]] || { echo "[$ID] controls mix PE and SE" >&2; exit 1; }
  {   # sourced by every step: each value shell-quoted (printf %q), never pasted in as code
    printf "SAMPLE=%q\nNUM_REP=%q\nNUM_CTL=%q\nCHIP_LABELS=%q\nCTL_LABELS=%q\nCTL_PE=%q\nCTL_MODE=%q\nRUNTIME=%q\n" \
      "$ID" "$N" "$M" "${CHIP_LABS[*]}" "${CTL_LABS[*]}" "${PE_OF[ctl]:0:1}" "$MODE" "$RUNTIME"
    printf "PROJECT_RESULTS=%q\nREF=%q\nCONDA_ENV=%q\n" "$PROJECT_RESULTS" "$REF" "$CONDA_ENV"
  } > "$SD/sample.env"
  echo "[$ID] $N ChIP rep(s), $M control rep(s), ctl_mode=$MODE, runtime=$RUNTIME -> $SD"
  (( DRY == 1 )) && continue

  # ---- snapshot the code: jobs run from $SD/code, so editing the pipeline later
  #      cannot break running jobs, and the exact code used stays with the results
  mkdir -p "$SD/code"; cp -r "$PIPE/config.sh" "$PIPE/lib.sh" "$PIPE/steps" "$SD/code/"
  CODE=$SD/code

  # ---- local: in this shell (--max-parallel 1) or in the background, at most --max-parallel samples at once
  if [[ $EXECUTOR == local ]]; then
    rm -f "$STATUS_DIR/$ID"; LOCAL_RUN+=("$ID")
    if (( MAXP == 1 )); then local_sample "$ID" "$SD" "$CODE"
    else
      while (( $(n_running) >= MAXP )); do wait -n 2>/dev/null || true; done
      local_sample "$ID" "$SD" "$CODE" &
      LOCAL_PIDS+=("$!")
      echo "[$(date '+%F %T')] [$ID] started ($(n_running)/$MAXP running, $THREADS threads each)"
    fi
    continue
  fi

  # ---- submit (align ∥ xcor) -> peaks -> finalize
  common=(--account="$ACCOUNT" --parsable)
  # slot: wait for the sample submitted MAXP places earlier to end (afterany: ok, failed or cancelled)
  slot=(); k=${#DONE_JOBS[@]}
  (( k >= MAXP )) && slot=(--dependency=afterany:${DONE_JOBS[k-MAXP]})
  na=$(wc -l < "$SD/align_units.txt"); nx=$(wc -l < "$SD/xcor_units.txt")
  j1=$(sbatch "${common[@]}" -J "al_$ID" -p "$ALIGN_PARTITION" -c "$ALIGN_CPU" --mem="$ALIGN_MEM" -t "$ALIGN_TIME" \
        --array=1-"$na" "${slot[@]}" -o "$SD/logs/align_%a_%A.log" "$CODE/steps/job.sh" "$CODE" align "$SD")
  j2=$(sbatch "${common[@]}" -J "xc_$ID" -p "$XCOR_PARTITION" -c "$XCOR_CPU" --mem="$XCOR_MEM" -t "$XCOR_TIME" \
        --array=1-"$nx" "${slot[@]}" -o "$SD/logs/xcor_%a_%A.log" "$CODE/steps/job.sh" "$CODE" xcor "$SD")
  j3=$(sbatch "${common[@]}" -J "pk_$ID" -p "$PEAK_PARTITION" -c "$PEAK_CPU" --mem="$PEAK_MEM" -t "$PEAK_TIME" \
        --dependency=afterok:$j1:$j2 --kill-on-invalid-dep=yes -o "$SD/logs/peaks_%j.log" "$CODE/steps/job.sh" "$CODE" peaks "$SD")
  j4=$(sbatch "${common[@]}" -J "fin_$ID" -p "$FINAL_PARTITION" -c "$FINAL_CPU" --mem="$FINAL_MEM" -t "$FINAL_TIME" \
        --dependency=afterok:$j3 --kill-on-invalid-dep=yes -o "$SD/logs/finalize_%j.log" "$CODE/steps/job.sh" "$CODE" finalize "$SD")
  DONE_JOBS+=("$j4")
  echo "[$ID] submitted align=$j1 xcor=$j2 peaks=$j3 finalize=$j4${slot:+  (starts after job ${DONE_JOBS[k-MAXP]})}  (logs: $SD/logs)"
done

if [[ $EXECUTOR == local && $DRY == 0 ]]; then
  wait || true   # background samples
  LOCAL_OK=(); LOCAL_FAIL=()
  for ID in "${LOCAL_RUN[@]}"; do
    if [[ $(cat "$STATUS_DIR/$ID" 2>/dev/null) == ok ]]; then LOCAL_OK+=("$ID"); else LOCAL_FAIL+=("$ID"); fi
  done
  echo "finished: ${#LOCAL_OK[@]} sample(s) ok${LOCAL_OK[*]:+ (${LOCAL_OK[*]})}, ${#LOCAL_FAIL[@]} failed${LOCAL_FAIL[*]:+ (${LOCAL_FAIL[*]})}"
  echo "results: $PROJECT_RESULTS   intermediates: $OUT_ROOT"
  # QC table of every sample in the sheet (python of the conda env, so the host needs no python)
  PY=python3; [[ $RUNTIME == conda ]] && PY=$CONDA_ENV/bin/python
  if command -v "$PY" >/dev/null; then
    env -i PATH=/usr/bin:/bin LC_ALL=C "$PY" "$PIPE/qc_aggregate.py" "$SAMPLES" --results "$PROJECT_RESULTS" \
      || echo "WARNING: qc_aggregate.py failed (results are not affected)" >&2
  else
    echo "WARNING: no python for qc_aggregate.py; run it by hand: python3 qc_aggregate.py $SAMPLES --results $PROJECT_RESULTS" >&2
  fi
  (( ${#LOCAL_FAIL[@]} == 0 ))
fi
