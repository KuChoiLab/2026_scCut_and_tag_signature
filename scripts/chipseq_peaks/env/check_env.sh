#!/bin/bash
# =============================================================================
# check_env.sh — check that a conda env can run the pipeline (run on the new server after installing)
#
# Usage: bash env/check_env.sh <ENV_PREFIX>
#
# Checks: tool versions (same as the env used for this project's results), picard.jar path,
#         run_spp.R / trimfastq.py present and identical (md5), R packages for run_spp.R,
#         host tools the steps take from /usr/bin and /bin (the pipeline runs with env -i,
#         PATH=<ENV_PREFIX>/bin:/usr/bin:/bin).
# Exit code 0 = everything OK.
# =============================================================================
set -uo pipefail
ENV=${1:?usage: bash check_env.sh <ENV_PREFIX>}
ENV=$(cd "$ENV" 2>/dev/null && pwd) || { echo "no such directory: $1" >&2; exit 1; }
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
fail=0
ok()  { printf "  OK    %s\n" "$*"; }
bad() { printf "  FAIL  %s\n" "$*"; fail=1; }
run() { env -i HOME="$HOME" LC_ALL=C PATH="$ENV/bin:/usr/bin:/bin" bash --noprofile --norc -c "$1" 2>&1; }

echo "== tool versions (expected = project env)"
chk() {  # name command expected-substring
  local out; out=$(run "$2" | head -1)
  [[ $out == *"$3"* ]] && ok "$1: $out" || bad "$1: got '$out', expected '$3'"
}
chk bowtie2   "bowtie2 --version"          "version 2.3.4.3"
chk samtools  "samtools --version"         "samtools 1.9"
chk macs2     "macs2 --version"            "macs2 2.2.4"
chk bedtools  "bedtools --version"         "bedtools v2.29.0"
chk python    "python --version"           "Python 3.6.7"
chk Rscript   "Rscript --version"          "version 3.4.1"
chk java      "java -version"              "1.8.0"
chk openssl   "openssl version"            "OpenSSL 1.0.2u"
for t in bedClip bedGraphToBigWig bedToBigBed; do
  [[ -x $ENV/bin/$t ]] && ok "$t" || bad "$t missing in $ENV/bin"
done

echo "== picard"
[[ -f $ENV/share/picard-2.20.7-0/picard.jar ]] && ok "share/picard-2.20.7-0/picard.jar" \
  || bad "missing $ENV/share/picard-2.20.7-0/picard.jar (picard=2.20.7 not installed?)"

echo "== scripts not on conda (see README_env.md)"
while read -r md5 f; do
  if [[ ! -f $ENV/bin/$f ]]; then bad "$f missing in $ENV/bin"
  elif [[ $(md5sum < "$ENV/bin/$f" | cut -c1-32) != "$md5" ]]; then bad "$f differs from the project copy (md5)"
  else ok "$f (md5 $md5)"; fi
done < "$HERE/extra_scripts.md5"

echo "== R packages used by run_spp.R"
out=$(run "Rscript -e 'suppressMessages({library(spp); library(caTools); library(snow); library(snowfall); library(bitops); library(Rsamtools)}); cat(\"loaded\")'")
[[ $out == *loaded* ]] && ok "spp caTools snow snowfall bitops Rsamtools" || bad "R packages: $out"

echo "== host tools from /usr/bin:/bin"
for t in awk sort shuf gzip zcat grep sed cut head tail wc mktemp tr paste find; do
  p=$(run "command -v $t") && ok "$t -> $p" || bad "$t not found"
done
run "echo a1 | grep -P 'a\d'" >/dev/null && ok "grep -P (PCRE)" || bad "grep -P not supported (GNU grep needed)"
run "seq 10 | shuf -n 2 --random-source=/dev/zero" >/dev/null && ok "shuf --random-source" || bad "shuf --random-source not supported (GNU coreutils needed)"
echo "  info  $(run 'sort --version' | head -1) | $(run 'awk --version' | head -1) | $(run 'gzip --version' | head -1)"

echo
(( fail == 0 )) && echo "ALL OK: $ENV" || echo "SOME CHECKS FAILED: $ENV"
exit $fail
