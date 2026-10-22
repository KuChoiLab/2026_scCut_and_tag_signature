# peak_count_matrix

Stage 2 of scCUTTag_signature. Builds one region × sample read count table
from the per-sample BAM and peak files made by stage 1 (`chipseq_peaks`).

DiffBind is used to:
1. merge all samples' peaks into consensus regions (kept if a peak in ≥ 2 samples)
2. count reads from every sample's BAM in every consensus region

All samples for one tissue × histone mark (e.g. all PBMC H3K4me3 samples)
go into a single run, so they share the same regions.

## Requirements

- conda (to build the environment below)
- stage 1 output for every sample

## Inputs (per sample, from stage 1 output)

| File | Path inside each sample's stage 1 output |
|---|---|
| BAM, with its `.bai` index in the same folder | `bam/chip_rep1.nodup.bam` |
| Peaks (the sample's own peaks) | `peaks/rep1.bfilt.narrowPeak.gz` |

Do not use `optimal_peak`, `conservative_peak`, or the `pr1`/`pr2` peak files.

## 1. Create the conda environment (once)

```bash
conda env create -f scripts/peak_count_matrix/environment.yml
conda activate peak_count_matrix
Rscript -e 'library(DiffBind); packageVersion("DiffBind")'   # should print 3.20.0
```

## 2. Write a sample sheet

CSV, one row per sample, all samples of the run in one file.
Column names must be exactly as below.

```
SampleID,Tissue,Factor,Condition,Replicate,bamReads,Peaks,PeakCaller
MINT_B_donor1,B,H3K4me3,MINT,1,/full/path/MINT_B_donor1/bam/chip_rep1.nodup.bam,/full/path/MINT_B_donor1/peaks/rep1.bfilt.narrowPeak.gz,narrow
MINT_B_donor2,B,H3K4me3,MINT,2,/full/path/MINT_B_donor2/bam/chip_rep1.nodup.bam,/full/path/MINT_B_donor2/peaks/rep1.bfilt.narrowPeak.gz,narrow
MINT_NK_donor1,NK,H3K4me3,MINT,1,/full/path/MINT_NK_donor1/bam/chip_rep1.nodup.bam,/full/path/MINT_NK_donor1/peaks/rep1.bfilt.narrowPeak.gz,narrow
```

| Column | Content |
|---|---|
| SampleID | unique sample name (use the stage 1 sample name); becomes the output column name |
| Tissue | cell type (e.g. B, NK, CD4T) |
| Factor | histone mark (e.g. H3K4me3) |
| Condition | ChIP platform: MINT (Mint-ChIP) or Stam (conventional ChIP) |
| Replicate | replicate number within a cell type and platform |
| bamReads | full path to the BAM |
| Peaks | full path to the peak file |
| PeakCaller | `narrow` |

Name the sheet `<tissue>_<mark>_samplesheet.csv`, e.g. `pbmc_h3k4me3_samplesheet.csv`.

## 3. Run

```bash
conda activate peak_count_matrix
scripts/peak_count_matrix/count_run.sh pbmc_h3k4me3_samplesheet.csv pbmc_h3k4me3_peak_counts.tsv 8 > pbmc_h3k4me3_peak_counts.log 2>&1
```

Arguments: sample sheet, output table, number of threads (optional; default
is `$SLURM_CPUS_PER_TASK` under Slurm, otherwise 1).

Run time and memory grow with the number and size of BAMs. A 2-sample test
ran within 16 GB. For many samples, run it on a machine or job with enough
memory (start around 32 GB) and check the log.

- With a job scheduler (Slurm, PBS, ...): put the `conda activate` and
  `count_run.sh` lines in a job script.
- Without a scheduler: run it so it survives logging out, e.g. add `nohup`
  in front of the command and `&` at the end.

## Output

- `<tissue>_<mark>_peak_counts.tsv`: tab-separated, one row per consensus
  region (CHR, START, END), one column of raw read counts per SampleID
  (not normalized)
- the log: run messages, ending with the R package versions

Please send back both files.

## Settings

- consensus regions: peak in ≥ 2 samples (DiffBind default)
- `summits = FALSE`: original peak widths kept
- `score = DBA_SCORE_READS`: raw counts
- single-end and paired-end BAMs are detected automatically
