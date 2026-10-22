# 2026_scCut_and_tag_signature

Choi Lab, Korea University College of Medicine.
Reference-guided cell-type annotation of scCUT&Tag data using bulk histone
ChIP-seq peak atlases (PBMC and liver). Current scope: PBMC, H3K4me3.

## Pipeline

| Stage | Folder | Input → Output | Run unit |
|---|---|---|---|
| 1 | `scripts/chipseq_peaks/` | FASTQ → BAM + peaks | per sample |
| 2 | `scripts/peak_count_matrix/` | stage 1 BAMs + peaks → region × sample read count table | all samples of one tissue × mark |

### Stage 1: `chipseq_peaks`

Bash reimplementation of ENCODE chip-seq-pipeline2 v2.2.2 (histone).
Runs locally (`chip_run_local.sh`) or with Slurm (`chip_run.sh`).
Reproduces ENCODE results file-for-file (md5) for single-end and paired-end
test samples; see `server_test/README_server_test.md`.

- Environment setup: `env/README_env.md`
- Defaults are in `config.sh`; on another server, override them with
  environment variables, e.g.
  `CONDA_ENV=/path/to/env REF=/path/to/hg38 bash chip_run_local.sh samples.txt -o OUTDIR -t 16`

### Stage 2: `peak_count_matrix`

DiffBind builds consensus regions from all samples' peaks and counts reads
from each sample's BAM. Own conda environment (`environment.yml`).
Full instructions, including the sample sheet format: `scripts/peak_count_matrix/README.md`.

## What is in this repository

Code, environment files, and documentation only. Data, intermediates, logs,
and results stay on each server and are excluded by `.gitignore`.

## Contributing

1. `git pull` before you start.
2. Add files by name (`git add <file>`), not `git add .`. Check
   `git status` and `git diff --staged` before committing.
3. To change a file someone else maintains, or this README, push a branch
   and open a pull request.
4. This repository is public: no data, personal paths, account names, or credentials.
5. Never `git push --force`.
6. A new stage gets its own folder `scripts/<stage>/` with a README and an `environment.yml`.
