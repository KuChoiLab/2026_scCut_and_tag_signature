# =============================================================================
# lib.sh — helper functions shared by the step scripts (sourced inside the tool environment)
# Each function re-implements one helper from chip-seq-pipeline2/src/encode_lib_*.py
# =============================================================================
set -euo pipefail
export LC_ALL=C

log()  { echo "[$(date '+%F %T')] $*" >&2; }
die()  { log "ERROR: $*"; exit 1; }

# Run a named step once: skip it when <dir>/.done.<name> exists (resume support).
step() {
  local dir=$1 name=$2; shift 2
  if [[ -e $dir/.done.$name ]]; then log "skip $name (done)"; return 0; fi
  log "start $name"
  "$@"
  touch "$dir/.done.$name"
  log "done  $name"
}

nlines() { zcat -f "$1" | wc -l; }                 # get_num_lines()

sort_mem() { echo "${SORT_MEM:-4G}"; }             # GNU sort -S buffer

# openssl-seeded random stream used for every shuf (identical to ENCODE)
rand_src() { openssl enc -aes-256-ctr -pass pass:"$1" -nosalt </dev/zero 2>/dev/null; }

# samtools sort with ENCODE's -@/-m convention (nth threads, mem per thread)
st_sort() { samtools sort "$@" -@ "$((NTH>1 ? NTH-1 : 0))" -m 1G; }

# --- subsample_ta_se(): shuffle-sample N reads (optionally drop chrM)
subsample_ta_se() {   # in N non_mito(0/1) out
  local ta=$1 n=$2 nonmito=$3 out=$4 seed
  seed=$(zcat -f "$ta" | wc -c)
  if [[ $nonmito == 1 ]]; then
    zcat -f "$ta" | grep -v "^${MITO_CHR}\b" | shuf -n "$n" --random-source=<(rand_src "$seed") | gzip -nc > "$out"
  else
    zcat -f "$ta" | shuf -n "$n" --random-source=<(rand_src "$seed") | gzip -nc > "$out"
  fi
}

# --- subsample_ta_pe(): keep mate pairs together (two lines -> one), sample N reads = N/2 pairs
subsample_ta_pe() {   # in N(reads, even) out
  local ta=$1 n=$2 out=$3 seed
  seed=$(zcat -f "$ta" | wc -c)
  zcat -f "$ta" | sed 'N;s/\n/\t/' | shuf -n "$(( n / 2 ))" --random-source=<(rand_src "$seed") \
    | awk 'BEGIN{OFS="\t"}{printf "%s\t%s\t%s\t%s\t%s\t%s\n%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' \
    | gzip -nc > "$out"
}

# --- bed_clip(): clip intervals to [0, chrom size]
bed_clip_gz() {  # in.bed out.gz
  bedClip "$1" "$CHRSZ" "$2.clip_tmp" -truncate -verbose=2
  gzip -nc "$2.clip_tmp" > "$2"; rm -f "$2.clip_tmp"
}

# --- blacklist_filter(): drop peaks overlapping blacklist, cap score at 1000, keep chr1-22XY
blacklist_filter() {  # peak.gz out.gz
  local tmpd; tmpd=$(mktemp -d -p "$(dirname "$2")")
  zcat -f "$1" > "$tmpd/p"; zcat -f "$BLACKLIST" > "$tmpd/b"
  if [[ ! -s $tmpd/p ]]; then
    : | gzip -nc > "$2"
  else
    bedtools intersect -nonamecheck -v -a "$tmpd/p" -b "$tmpd/b" \
      | awk 'BEGIN{OFS="\t"}{if ($5>1000) $5=1000; print $0}' \
      | { grep -P "${REGEX_BFILT_PEAK_CHR}\b" || true; } | gzip -nc > "$2"
  fi
  rm -rf "$tmpd"
}

# --- frip_shifted(): fraction of reads (shifted by fraglen/2) falling in peaks
frip_shifted() {  # ta peak fraglen out.qc
  local ta=$1 peak=$2 fraglen=$3 out=$4 half n_in n_all tmp
  half=$(( (fraglen + 1) / 2 ))
  n_all=$(nlines "$ta")
  if [[ $(nlines "$peak") -eq 0 ]]; then n_in=0; else
    tmp=$(mktemp -p "$(dirname "$out")"); zcat -f "$peak" > "$tmp"
    n_in=$(bedtools slop -i "$ta" -g "$CHRSZ" -s -l "-$half" -r "$half" \
      | awk '{if ($2>=0 && $3>=0 && $2<=$3) print $0}' \
      | bedtools intersect -nonamecheck -a stdin -b "$tmp" -wa -u | wc -l)
    rm -f "$tmp"
  fi
  awk -v a="$n_in" -v b="$n_all" 'BEGIN{print a/b}' > "$out"
}

# --- pool_ta(): concatenate tagAligns
pool_ta() {  # out in1 in2 ...
  local out=$1; shift
  zcat -f "$@" | gzip -nc > "$out"
}
