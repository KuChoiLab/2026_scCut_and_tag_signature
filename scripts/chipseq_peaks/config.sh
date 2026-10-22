# =============================================================================
# config.sh — settings shared by every step (sourced, not executed)
# Defaults reproduce ENCODE chip-seq-pipeline2 v2.2.2 (pipeline_type=histone).
# =============================================================================

# ---- tool environment: all tools (bowtie2, samtools, picard, macs2, bedtools, spp ...)
#   conda: env built from env/chipseq_v2.2.2.yml (same tool versions as the sif), isolated with env -i
#   sif  : ENCODE chip-seq-pipeline v2.2.2 container
# Paths below are this cluster's; on another system set them in the environment (CONDA_ENV=... REF=... bash ...)
RUNTIME=${RUNTIME:-conda}      # conda | sif   (chip_run.sh --runtime, recorded in sample.env)
CONDA_ENV=${CONDA_ENV:-$HOME/conda_envs/chipseq_v2.2.2}
SIF=${SIF:-}                   # path to chip-seq-pipeline_v2.2.2.sif (needed only for --runtime sif)
BINDS=${BINDS:-$HOME}          # folders the container can see (--runtime sif)
case $RUNTIME in
  conda) PICARD_JAR=$CONDA_ENV/share/picard-2.20.7-0/picard.jar ;;
  sif)   PICARD_JAR=/software/picard.jar ;;
esac

# ---- genome (values from reference/hg38/hg38.tsv)
REF=${REF:-$HOME/scCUTTag_signature/data/raw/reference/hg38}
BOWTIE2_IDX_TAR=$REF/bowtie2_index/ENCFF110MCL.tar.gz
BOWTIE2_IDX_DIR=$REF/bowtie2_index/ENCFF110MCL_untar   # extracted once by chip_run.sh
CHRSZ=$REF/GRCh38_EBV.chrom.sizes.tsv
BLACKLIST=$REF/ENCFF356LFX.bed.gz
GENSZ=hs
MITO_CHR=chrM
bt2_index_ok() {   # 0 when all six bowtie2 index files are in BOWTIE2_IDX_DIR (a cut-off extraction leaves some out)
  local m=("$BOWTIE2_IDX_DIR"/*.rev.1.bt2*) idx s
  [[ -e ${m[0]} ]] || return 1
  idx=${m[0]%.rev.1.bt2*}
  for s in 1 2 3 4 rev.1 rev.2; do [[ -s $idx.$s.bt2 || -s $idx.$s.bt2l ]] || return 1; done
}
REGEX_BFILT_PEAK_CHR='chr[\dXY]+'          # peaks kept only on chr1-22,X,Y

# ---- output
OUT_ROOT=${OUT_ROOT:-$HOME/scCUTTag_signature/work}         # all intermediates (scratch)
PROJECT_RESULTS=${PROJECT_RESULTS:-$HOME/scCUTTag_signature/data/processed/chipseq} # finals copied here (backed up)

# ---- alignment / filtering (ENCODE defaults)
MAPQ_THRESH=30                 # samtools view -q
FILTER_CHRS=""                 # e.g. "chrM" to drop reads on these chrs from the final BAM (ENCODE chip default: none)
XCOR_TRIM_BP=50                # R1 trimmed to 50 bp and re-aligned SE for cross-correlation
XCOR_SUBSAMPLE_READS=15000000  # reads used by run_spp.R
XCOR_EXCLUSION_RANGE_MIN=-500

# ---- control selection (ENCODE choose_ctl)
CTL_MODE=pooled                # pooled (always_use_pooled_ctl=true) | matched
CTL_DEPTH_RATIO=1.2
CTL_DEPTH_LIMIT=200000000
EXP_CTL_DEPTH_RATIO_LIMIT=5.0

# ---- peak calling (ENCODE macs2 for histone)
PVAL_THRESH=0.01
CAP_NUM_PEAK=500000
PSEUDOREP_SEED=0               # 0 = use tagAlign size as seed (ENCODE default)

# ---- Slurm resources (each can be overridden from the environment, e.g. ALIGN_CPU=4 bash chip_run.sh ...)
MAX_PARALLEL=${MAX_PARALLEL:-3}  # samples running at the same time (chip_run.sh --max-parallel)
ACCOUNT=${ACCOUNT:-choilab}
ALIGN_PARTITION=${ALIGN_PARTITION:-debug}; ALIGN_CPU=${ALIGN_CPU:-16}; ALIGN_MEM=${ALIGN_MEM:-128G}; ALIGN_TIME=${ALIGN_TIME:-23:00:00}
XCOR_PARTITION=${XCOR_PARTITION:-debug};   XCOR_CPU=${XCOR_CPU:-8};    XCOR_MEM=${XCOR_MEM:-32G};   XCOR_TIME=${XCOR_TIME:-12:00:00}
PEAK_PARTITION=${PEAK_PARTITION:-debug};   PEAK_CPU=${PEAK_CPU:-8};    PEAK_MEM=${PEAK_MEM:-64G};   PEAK_TIME=${PEAK_TIME:-23:00:00}
FINAL_PARTITION=${FINAL_PARTITION:-debug}; FINAL_CPU=${FINAL_CPU:-1};  FINAL_MEM=${FINAL_MEM:-4G};  FINAL_TIME=${FINAL_TIME:-04:00:00}
