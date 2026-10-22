#!/bin/bash
# =============================================================================
# setup_env.sh — build the pipeline's conda env on a new server (one command)
#
# Usage: bash env/setup_env.sh [ENV_PREFIX] [--miniconda DIR]
#   ENV_PREFIX   where the env goes (default ~/conda_envs/chipseq_v2.2.2)
#   --miniconda  where to install Miniconda if no conda is found (default ~/miniconda3)
#
# Steps (each is skipped when already done, so re-running is safe):
#   1. server check: linux x86_64, CPUs, RAM, free disk (warnings only)
#   2. conda: use the one on PATH (or in --miniconda DIR); install Miniconda there if none
#   3. env from chipseq_v2.2.2.lock.txt (exact packages, no solver)
#   4. run_spp.R and trimfastq.py (not on conda) downloaded into ENV/bin, md5 checked
#   5. check_env.sh
# =============================================================================
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ENV=""; MC=$HOME/miniconda3
while [[ $# -gt 0 ]]; do
  case $1 in
    --miniconda) MC=$2; shift ;;
    -h|--help) sed -n '4,15p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) ENV=$1 ;;
  esac; shift
done
ENV=${ENV:-$HOME/conda_envs/chipseq_v2.2.2}

echo "######## 1/5 server"
arch=$(uname -m); cpu=$(nproc); ram=$(awk '/MemTotal/{printf "%d", $2/1024/1024}' /proc/meminfo)
mkdir -p "$(dirname "$ENV")"
disk=$(df -Pk "$(dirname "$ENV")" | awk 'NR==2{printf "%d", $4/1024/1024}')
echo "  $(uname -s) $arch, $cpu CPUs, ${ram} GB RAM, ${disk} GB free at $(dirname "$ENV")"
[[ $(uname -s) == Linux && $arch == x86_64 ]] || { echo "  ERROR: needs Linux x86_64 (lock file is linux-64)" >&2; exit 1; }
(( ram >= 60 )) || echo "  WARNING: < 64 GB RAM: use fewer threads (-t 8) when running the pipeline"
(( cpu >= 8 ))  || echo "  WARNING: < 8 CPUs: the test sample will take much longer"
(( disk >= 5 )) || echo "  WARNING: < 5 GB free here: the env needs ~3 GB"

echo "######## 2/5 conda"
if command -v conda >/dev/null; then CONDA=$(command -v conda)
elif [[ -x $MC/bin/conda ]]; then CONDA=$MC/bin/conda
else
  echo "  no conda found -> installing Miniconda to $MC"
  inst=$(mktemp --suffix=.sh)
  curl -fL --retry 5 -o "$inst" https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
  bash "$inst" -b -p "$MC"; rm -f "$inst"
  CONDA=$MC/bin/conda
fi
echo "  using $CONDA ($("$CONDA" --version))"

echo "######## 3/5 env $ENV"
if [[ -x $ENV/bin/macs2 && -x $ENV/bin/bowtie2 && -x $ENV/bin/Rscript ]]; then
  echo "  already there (delete $ENV to rebuild)"
elif [[ -e $ENV && -n $(ls -A "$ENV" 2>/dev/null) ]]; then
  # never touch a folder this script did not create: it may be another env or someone's data
  cat >&2 <<EOF
  ERROR: $ENV already exists, is not empty, and is not a finished chipseq environment.
  Nothing was changed. Give a new path (bash env/setup_env.sh /path/to/new/chipseq_v2.2.2), or, if this is
  a half-built environment left by an interrupted setup_env.sh run, delete it yourself and run again.
EOF
  exit 1
else
  created=0; [[ -e $ENV ]] || created=1      # remove on failure only what this run created
  if ! "$CONDA" create -y -p "$ENV" --file "$HERE/chipseq_v2.2.2.lock.txt"; then
    cat >&2 <<EOF

  conda create failed. Most common cause on recent Miniconda: the Anaconda channel terms of service
  (3 packages in the lock file come from repo.anaconda.com: libopenblas, wget, libcurl). If the error mentions Terms of Service, run
  the commands it lists, usually both:
      $CONDA tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
      $CONDA tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r
  and run this script again. If that channel is blocked at your institution, send the error message.
EOF
    (( created )) && rm -rf "$ENV"
    exit 1
  fi
fi

echo "######## 4/5 ENCODE scripts not on conda"
want() { awk -v f="$1" '$2==f{print $1}' "$HERE/extra_scripts.md5"; }
have() { [[ -f $ENV/bin/$1 ]] && md5sum < "$ENV/bin/$1" | cut -c1-32 || echo none; }
if [[ $(have run_spp.R) != "$(want run_spp.R)" ]]; then
  curl -fL --retry 5 https://github.com/kundajelab/phantompeakqualtools/archive/refs/tags/1.2.1.tar.gz \
    | tar -xz -O phantompeakqualtools-1.2.1/run_spp.R > "$ENV/bin/run_spp.R"
fi
if [[ $(have trimfastq.py) != "$(want trimfastq.py)" ]]; then
  curl -fL --retry 5 -o "$ENV/bin/trimfastq.py" \
    https://raw.githubusercontent.com/ENCODE-DCC/chip-seq-pipeline2/v2.2.2/src/trimfastq.py
fi
chmod +x "$ENV/bin/run_spp.R" "$ENV/bin/trimfastq.py"
echo "  run_spp.R, trimfastq.py in $ENV/bin"

echo "######## 5/5 check"
bash "$HERE/check_env.sh" "$ENV"
echo
echo "next: bash server_test/run_test.sh --env $ENV --ref <REF_DIR> --data <DATA_DIR> --out <OUTDIR> -t <CPUs>"
