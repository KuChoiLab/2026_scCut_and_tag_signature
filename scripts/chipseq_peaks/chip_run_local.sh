#!/bin/bash
# =============================================================================
# chip_run_local.sh — run the pipeline without Slurm, in this shell (one machine / one allocation)
#
# Usage:
#   bash chip_run_local.sh <samples.txt> -o OUTDIR [-t THREADS] [-j SAMPLES] [--runtime conda|sif]
#                          [--only ID[,ID]] [--ctl-mode pooled|matched] [--dry-run]
#
#   -o OUTDIR   everything goes here:
#                 OUTDIR/work/<sample>/     intermediates, logs, code snapshot (resumable)
#                 OUTDIR/results/<sample>/  final BAMs, peaks, bigwigs, QC (same layout as the Slurm run)
#   -t THREADS  CPUs per step (default 8). Measured on MINT samples (25-90M read pairs):
#               align up to ~55G with 16 threads (samtools sort uses 1G per thread + Picard 8G),
#               xcor up to ~18G and peaks up to ~36G with 8 threads -> use a machine with >= 64G RAM.
#               One sample takes ~4-5 h with 16 threads; intermediates need up to ~40G per sample.
#   -j SAMPLES  samples run at the same time (default 1). Uses SAMPLES x THREADS CPUs and up to
#               SAMPLES x ~55G RAM (all aligning at once with 16 threads): e.g. -j 3 -t 16 -> 48 CPUs, >= 180G RAM.
#   --runtime   conda (default: env built from env/chipseq_v2.2.2.yml at CONDA_ENV) or sif (Singularity image SIF)
#
# Within a sample the steps run in order
#   align (ChIP, control) -> xcor -> peaks -> finalize
# With -j 1 samples run one after another; with -j N a new sample starts as soon as one finishes.
# Each sample logs to OUTDIR/work/<sample>/logs/. A failed sample is reported and the remaining
# samples still run; the exit code is non-zero if any sample failed. Running the same command again resumes from the first unfinished step.
#
# On another system, point config.sh at local installs through the environment, e.g.
#   REF=/data/hg38 CONDA_ENV=/opt/envs/chipseq_v2.2.2 bash chip_run_local.sh samples.txt -o out -t 16
# =============================================================================
set -euo pipefail
PIPE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
usage() { sed -n '4,18p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//' >&2; exit 1; }

[[ $# -ge 1 && $1 != -* ]] || usage
SAMPLES=$1; shift
OUTDIR=""; ARGS=()
while [[ $# -gt 0 ]]; do
  case $1 in
    -o|--outdir) OUTDIR=$2; shift ;;
    -t|--threads) ARGS+=(--threads "$2"); shift ;;
    -j|--jobs) ARGS+=(--max-parallel "$2"); shift ;;
    --runtime|--only|--ctl-mode) ARGS+=("$1" "$2"); shift ;;
    --dry-run) ARGS+=(--dry-run) ;;
    -h|--help) usage ;;
    *) echo "unknown option $1" >&2; usage ;;
  esac; shift
done
[[ -n $OUTDIR ]] || { echo "-o OUTDIR is required" >&2; usage; }

exec bash "$PIPE/chip_run.sh" "$SAMPLES" --executor local \
  --out "$OUTDIR/work" --results "$OUTDIR/results" "${ARGS[@]}"
