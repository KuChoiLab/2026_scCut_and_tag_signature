# 다른 서버에서 테스트 (스케줄러 없이, conda)

Yale에서 이미 돌린 실제 샘플 2개를 새 서버에서 다시 돌리고, 결과가 **파일 단위(md5)로 같은지** 확인합니다.

| 샘플 | 종류 | FASTQ | 비교 기준 |
|---|---|---|---|
| `Stam_B_028NGF` | 일반 ChIP-seq, single-end 36 bp | 4개, 2.7 GB | Yale 결과 + ENCODE 원본 파이프라인(caper/WDL) 결과 |
| `Mint_naiveCD4_836ICQ` | MINT-ChIP, paired-end 76 bp | 40개(lane 20쌍), 8.0 GB | Yale 결과 |

## 준비물
- linux-64, CPU 16개 이상, RAM 64 GB 이상 (두 샘플 동시 실행 `-j 2`는 120 GB 이상)
- 디스크 약 100 GB (참조 8 + FASTQ 11 + 중간 파일 샘플당 ~40 + 결과)
- 인터넷 (ENCODE 포털에서 참조·FASTQ 다운로드)
- conda 환경: `bash env/setup_env.sh` 결과가 `ALL OK`

## 실행 (한 줄)
연구실 HPC에서는 head node(node00)에서 직접 돌리지 말고 `sbatch ~/scCUTTag_signature/slurm/server_test.sbatch <샘플>`을 사용합니다. 아래는 일반 서버용 예시입니다.
```bash
cd ~/scCUTTag_signature/scripts/chipseq_peaks
nohup bash server_test/run_test.sh \
  --env  ~/conda_envs/chipseq_v2.2.2 \
  --ref  ~/scCUTTag_signature/data/raw/reference/hg38 \
  --data ~/scCUTTag_signature/data/raw/server_test_fastq \
  --out  ~/scCUTTag_signature/work/server_test \
  -t 16 -j 2 > run_test.log 2>&1 &
tail -f run_test.log
```
- `-j 2`: 두 샘플을 동시에 (RAM 120 GB 이상). RAM이 적으면 `-j 1`(하나씩, 시간 2배)
- 한 샘플만: `--samples Stam_B_028NGF` 또는 `--samples Mint_naiveCD4_836ICQ`
- 진행 순서: 환경 점검 → 참조 다운로드 → FASTQ 다운로드(md5 확인) → 파이프라인(샘플당 약 4–5시간) → 결과 비교
- 중간에 끊기면 **같은 명령을 다시 실행**하면 됩니다. 끝난 단계와 받은 파일은 건너뜁니다.

단계별로 따로 실행할 수도 있습니다.
```bash
bash server_test/get_reference.sh  ~/scCUTTag_signature/data/raw/reference/hg38
bash server_test/get_test_fastq.sh ~/scCUTTag_signature/data/raw/server_test_fastq            # -> samples_test.txt 생성 (두 샘플)
CONDA_ENV=~/conda_envs/chipseq_v2.2.2 REF=~/scCUTTag_signature/data/raw/reference/hg38 \
  bash chip_run_local.sh ~/scCUTTag_signature/data/raw/server_test_fastq/samples_test.txt -o ~/scCUTTag_signature/work/server_test -t 16 -j 2
CONDA_ENV=~/conda_envs/chipseq_v2.2.2 bash server_test/verify.sh ~/scCUTTag_signature/work/server_test/results/Stam_B_028NGF
CONDA_ENV=~/conda_envs/chipseq_v2.2.2 bash server_test/verify.sh ~/scCUTTag_signature/work/server_test/results/Mint_naiveCD4_836ICQ
```

## 결과 판정
`verify.sh`가 샘플마다 27개 항목을 확인합니다: 정렬 read(tagAlign, BAM), bowtie2·flagstat·중복·라이브러리 복잡도 QC,
cross-correlation, peak 6종(최종 `optimal_peak` 포함), reproducibility, QC 표, signal bigwig 2종.
26개는 md5로 비교하고, cross-correlation(`xcor/rep1.cc.qc`)은 값만 보여 줍니다(`info`, not compared로 집계).
- 정상 결과: 샘플마다 `26 same, 0 different, 0 skipped, 1 not compared (of 27)` + `PASS`, 마지막 줄 `ALL TEST SAMPLES PASS`
- cc.qc를 비교하지 않는 이유: 4번째 열(fragment length 위치의 상관값)이 실행마다 소수 7–9째 자리에서 달라집니다.
  bowtie2를 여러 스레드로 돌리면 출력 순서가 매번 다르고, 중복 제거 없이 정렬한 xcor용 read에서 좌표가 같은 read의 순서가 바뀌어
  `shuf`가 15M 중 몇 개를 다르게 뽑기 때문입니다(ENCODE 원본도 같음). fraglen·NSC·RSC는 `qc/align_qc.tsv`에서 그대로 비교되고,
  fraglen이 다르면 peak 파일이 모두 달라지므로 검증력은 줄지 않습니다.
- `DIFF` 항목이 있으면 → `run_test.log` 전체와 해당 샘플의 `qc/align_qc.tsv`, `qc/peak_qc.tsv`, `xcor/rep1.cc.qc`를 보내 주세요.
  어디서 갈라졌는지(정렬 / 중복 제거 / pseudo-replicate / peak) 앞에서부터 보면 알 수 있습니다.

기준값 출처 (`expected_<샘플>.tsv`의 source 열):
- `yale`: Yale 클러스터 실행 결과(2026-10-01). Stam_B는 ENCODE 원본 파이프라인 결과와도 peak·QC·fc bigwig가 모두 같습니다.
- `encode`: Stam_B와 Mint의 pval bigwig. Yale 파일은 2026-10-02 수정(`peaks_sample.sh` sval 자릿수) 전에 만들어져 ENCODE 값을 기준으로 씁니다.
- `info`: cross-correlation `xcor/rep1.cc.qc` (위 설명). fraglen, NSC, RSC 값만 출력
- `none` (n/a 표시): 이전 버전 expected 파일의 Mint pval bigwig (ENCODE 결과가 생기기 전). 지금 파일에는 ENCODE 기준이 들어 있습니다.

## 결과가 달라질 수 있는 곳 (서버 차이)
파이프라인은 conda 도구 외에 서버 기본 도구(GNU `sort`, `shuf`, `awk`, `gzip`)를 씁니다. pseudo-replicate는 `shuf`의 난수 순서에 의존하므로 coreutils 버전이 크게 다르면 pseudo-rep peak부터 달라질 수 있습니다.
`check_env.sh`가 출력하는 버전 줄(Yale: coreutils 8.30, gawk 4.2.1, gzip 1.9)을 함께 보내 주세요.
