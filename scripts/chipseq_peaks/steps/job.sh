#!/bin/bash
# =============================================================================
# steps/job.sh — run one step script in the tool environment (conda env or container)
#   Slurm : sbatch ... steps/job.sh <pipeline_dir> <step: align|xcor|peaks|finalize> <sample_dir>
#   local : NTH=16 bash steps/job.sh <pipeline_dir> <step> <sample_dir> [task_number]
# align/xcor: the task number (4th arg, else SLURM_ARRAY_TASK_ID) picks the line of the step's unit list.
# =============================================================================
set -euo pipefail
PIPE=$1; STEP=$2; SD=$3; TASK=${4:-${SLURM_ARRAY_TASK_ID:-}}
source "$SD/sample.env"                                # RUNTIME chosen at submission
source "$PIPE/config.sh"

case $STEP in
  align|xcor)
    [[ -n $TASK ]] || { echo "$STEP needs a task number (4th argument or Slurm array)" >&2; exit 1; }
    LINE=$(sed -n "${TASK}p" "$SD/${STEP}_units.txt")   # line number in units.tsv
    ARGS=("$PIPE" "$SD" "$LINE") ;;
  peaks|finalize)
    ARGS=("$PIPE" "$SD") ;;
  *) echo "unknown step $STEP" >&2; exit 1 ;;
esac

NTH=${NTH:-${SLURM_CPUS_PER_TASK:-1}}                 # local runs pass NTH, Slurm jobs use their CPU count

# ---- check the reference and tool files this step needs (REF / CONDA_ENV come from sample.env, written by chip_run.sh)
if [[ $STEP != finalize ]]; then
  miss=()
  for f in "$CHRSZ" "$BLACKLIST"; do [[ -s $f ]] || miss+=("$f"); done
  [[ $STEP == align || $STEP == xcor ]] && { bt2_index_ok || miss+=("bowtie2 index (6 files) in $BOWTIE2_IDX_DIR"); }
  if [[ $RUNTIME == conda ]]; then
    for f in "$CONDA_ENV/bin/macs2" "$PICARD_JAR"; do [[ -e $f ]] || miss+=("$f"); done
  fi
  if (( ${#miss[@]} )); then
    echo "ERROR: $STEP cannot start, missing:" >&2; printf "  %s\n" "${miss[@]}" >&2
    echo "  (REF=$REF  CONDA_ENV=$CONDA_ENV  from $SD/sample.env)" >&2
    exit 1
  fi
fi
echo "[$(date '+%F %T')] $STEP ${TASK:+task $TASK }on $(hostname), $NTH CPUs, runtime=$RUNTIME"
SCRIPT=$PIPE/steps/${STEP}_$( [[ $STEP == peaks ]] && echo sample || echo unit ).sh
if [[ $STEP == finalize ]]; then
  bash "$PIPE/steps/finalize.sh" "${ARGS[@]}"          # plain copy, no tool environment needed
elif [[ $RUNTIME == conda ]]; then
  # env -i keeps host modules/conda/.bashrc out (like --cleanenv); only the conda env + system base tools on PATH.
  # REF and CONDA_ENV must be passed on: the step re-reads config.sh, whose defaults are this cluster's paths.
  env -i HOME="$HOME" USER="$USER" TMPDIR="${TMPDIR:-/tmp}" LC_ALL=C NTH="$NTH" RUNTIME=conda \
    REF="$REF" CONDA_ENV="$CONDA_ENV" ${SORT_MEM:+SORT_MEM="$SORT_MEM"} ${JAVA_HEAP:+JAVA_HEAP="$JAVA_HEAP"} \
    CONDA_PREFIX="$CONDA_ENV" PATH="$CONDA_ENV/bin:/usr/bin:/bin" \
    bash --noprofile --norc "$SCRIPT" "${ARGS[@]}"
elif [[ $RUNTIME == sif ]]; then
  # --cleanenv keeps host modules/conda out of the container; NTH passes the CPU count in, REF the reference location
  singularity exec --cleanenv --env NTH="$NTH",RUNTIME=sif,REF="$REF"${SORT_MEM:+,SORT_MEM=$SORT_MEM}${JAVA_HEAP:+,JAVA_HEAP=$JAVA_HEAP} \
    -B "$BINDS" "$SIF" bash --noprofile --norc "$SCRIPT" "${ARGS[@]}"
else
  echo "unknown RUNTIME $RUNTIME (conda|sif)" >&2; exit 1
fi
echo "[$(date '+%F %T')] $STEP finished"
