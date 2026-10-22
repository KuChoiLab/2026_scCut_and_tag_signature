#!/bin/bash
# Usage: count_run.sh samplesheet.csv peak_counts.tsv [threads]
set -euo pipefail
SHEET=$1
OUT=$2
THREADS=${3:-${SLURM_CPUS_PER_TASK:-1}}
DIR=$(cd "$(dirname "$0")" && pwd)
Rscript "$DIR/steps/count_reads.R" "$SHEET" "$OUT" "$THREADS"
