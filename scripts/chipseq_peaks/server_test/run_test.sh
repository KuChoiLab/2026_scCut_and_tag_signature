#!/bin/bash
# =============================================================================
# run_test.sh — full test on a new server: env check -> reference -> test FASTQs -> pipeline -> verify
#
# Usage: bash server_test/run_test.sh --env ENV_PREFIX --ref REF_DIR --data DATA_DIR --out OUTDIR
#                                     [--samples S[,S]] [-t THREADS] [-j SAMPLES]
#
#   --env      conda env built with env/setup_env.sh (see env/README_env.md)
#   --ref      where the hg38 files go (downloaded once, ~3.5 GB + ~4 GB extracted index)
#   --data     where the test FASTQs go (Stam_B_028NGF 2.7 GB, Mint_naiveCD4_836ICQ 8.0 GB)
#   --out      pipeline output (~40 GB intermediates per sample in OUTDIR/work, finals in OUTDIR/results)
#   --samples  test samples (default both): Stam_B_028NGF (single-end, also checked against the ENCODE
#              pipeline itself), Mint_naiveCD4_836ICQ (paired-end, MINT-ChIP)
#   -t         CPUs per step (default 16); -j samples at the same time (default 1)
#              RAM needed ~ j x 55 GB, e.g. -t 16 -j 2 -> 32 CPUs, >= 120 GB RAM. ~4-5 h per sample with 16 CPUs.
# Every step resumes: running the same command again skips what is already done.
# Long run: start it with nohup (or inside tmux/screen) so it survives logging out, e.g.
#   nohup bash server_test/run_test.sh --env ... > run_test.log 2>&1 &
# =============================================================================
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd); PIPE=$(dirname "$HERE")
usage() { sed -n '4,19p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//' >&2; exit 1; }
ENV=""; REFD=""; DATA=""; OUT=""; T=16; J=1
SAMPLES=$(awk -F'\t' '!/^#/ && !s[$1]++ {printf "%s%s", (n++ ? "," : ""), $1}' "$HERE/test_samples.tsv")
while [[ $# -gt 0 ]]; do
  case $1 in
    --env) ENV=$2; shift ;; --ref) REFD=$2; shift ;; --data) DATA=$2; shift ;; --out) OUT=$2; shift ;;
    --samples) SAMPLES=$2; shift ;;
    -t|--threads) T=$2; shift ;;
    -j|--jobs) J=$2; shift ;;
    *) usage ;;
  esac; shift
done
[[ -n $ENV && -n $REFD && -n $DATA && -n $OUT ]] || usage
ENV=$(cd "$ENV" && pwd)

echo "######## 1/5 conda env check";   bash "$PIPE/env/check_env.sh" "$ENV"
echo "######## 2/5 reference";          bash "$HERE/get_reference.sh" "$REFD"; REFD=$(cd "$REFD" && pwd)
echo "######## 3/5 test FASTQs";        bash "$HERE/get_test_fastq.sh" "$DATA" "$SAMPLES"; DATA=$(cd "$DATA" && pwd)
echo "######## 4/5 pipeline (${SAMPLES//,/ }; $T CPUs per step, $J sample(s) at a time)"
CONDA_ENV=$ENV REF=$REFD bash "$PIPE/chip_run_local.sh" "$DATA/samples_test.txt" -o "$OUT" -t "$T" -j "$J" --runtime conda || true
echo "######## 5/5 verify"
fail=0
for s in ${SAMPLES//,/ }; do
  CONDA_ENV=$ENV bash "$HERE/verify.sh" "$(cd "$OUT" && pwd)/results/$s" || fail=1
  echo
done
(( fail == 0 )) && echo "ALL TEST SAMPLES PASS" || echo "SOME TEST SAMPLES FAILED (see above)"
exit $fail
