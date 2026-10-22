# conda 환경 옮기기 (chipseq_v2.2.2)

이 프로젝트의 모든 결과는 이 환경으로 만들었습니다. 다른 서버에서 같은 결과를 내려면 같은 환경을 만들면 됩니다.

## 파일
| 파일 | 내용 | 용도 |
|---|---|---|
| `chipseq_v2.2.2.lock.txt` | 모든 패키지의 정확한 다운로드 주소 + md5 (`@EXPLICIT`, linux-64) | **권장.** 패키지 해석(solver) 없이 완전히 같은 환경 |
| `environment.yml` | 모든 패키지의 버전 + 빌드 고정 | lock을 쓸 수 없을 때 (conda가 패키지를 다시 해석) |
| `chipseq_v2.2.2.yml` | 주요 도구 버전만 | 참고용 (어떤 도구가 왜 들어갔는지 주석) |
| `extra_scripts.md5` | conda에 없는 스크립트 2개의 md5 | `check_env.sh`가 사용 |
| `check_env.sh` | 설치 후 점검 | 도구 버전, picard, 스크립트, R 패키지, 서버 기본 도구 |

## 한 번에 설치: `bash env/setup_env.sh [ENV_PREFIX]`
서버 확인 → conda 찾기(없으면 Miniconda 설치) → lock으로 환경 생성 → 스크립트 2개 받기 → `check_env.sh`.
다시 실행해도 끝난 단계는 건너뜁니다. 아래는 같은 과정을 손으로 하는 방법입니다.

## 설치 (새 서버, linux-64, 인터넷 연결)
```bash
ENV=/path/to/conda_envs/chipseq_v2.2.2        # 원하는 위치

# 1) 환경 생성 (둘 중 하나)
conda create -p $ENV --file env/chipseq_v2.2.2.lock.txt      # 권장
# conda env create -f env/environment.yml -p $ENV           # lock이 안 될 때

# 2) conda에 없는 ENCODE 스크립트 2개를 env/bin에 넣기 (ENCODE 컨테이너 v2.2.2와 같은 원본)
curl -L https://github.com/kundajelab/phantompeakqualtools/archive/refs/tags/1.2.1.tar.gz \
  | tar -xz -O phantompeakqualtools-1.2.1/run_spp.R > $ENV/bin/run_spp.R
curl -L -o $ENV/bin/trimfastq.py \
  https://raw.githubusercontent.com/ENCODE-DCC/chip-seq-pipeline2/v2.2.2/src/trimfastq.py
chmod +x $ENV/bin/run_spp.R $ENV/bin/trimfastq.py

# 3) 점검 -> 마지막 줄이 "ALL OK" 이어야 함
bash env/check_env.sh $ENV
```

## 파이프라인에 환경 위치 알려주기
`config.sh`의 기본값은 연구실 HPC 경로입니다 (2026-10 변경, 원래 Yale 경로). 다른 서버에서는 실행할 때 환경 변수로 바꿉니다.
```bash
CONDA_ENV=$ENV REF=/path/to/hg38 bash chip_run_local.sh samples.txt -o OUTDIR -t 16
```

## 주의
- **activate 하지 않아도 됩니다.** 파이프라인이 `env -i PATH=$ENV/bin:/usr/bin:/bin`으로 이 환경만 써서 실행합니다. 그래서 서버의 기본 도구(`awk`, `sort`, `shuf`, `gzip`, `grep -P` 등, GNU 버전)도 사용합니다. `check_env.sh`가 이 도구들도 확인합니다.
- lock 파일의 패키지 3개(`libcurl`, `libopenblas`, `wget`)는 `repo.anaconda.com`(defaults 채널)에서 받습니다. 기관 방화벽이나 Anaconda 이용 약관 때문에 막히면 알려 주세요. conda-forge 대체 패키지를 찾아보겠습니다.
- 환경 폴더를 통째로 복사하는 방법(`cp -r`, `conda-pack`)도 됩니다. 단, 경로가 바뀌면 일부 스크립트의 첫 줄(`#!`)이 옛 경로를 가리킬 수 있습니다. 위처럼 새로 만드는 쪽이 안전합니다.
- Yale에서 측정한 버전(기준): bowtie2 2.3.4.3, samtools 1.9, macs2 2.2.4, bedtools 2.29.0, picard 2.20.7, python 3.6.7, R 3.4.1, openjdk 1.8.0_472, openssl 1.0.2u.
  서버 기본 도구: GNU coreutils 8.30, GNU awk 4.2.1, gzip 1.9 (RHEL 8).
