> **논문 본문 범위는 `CORE_PAPER.md`, 후속 연구 후보는 `FOLLOW_UP_PAPERS.md`를 먼저 읽는다.**

# v2 — 270개 전수 계산 파이프라인 (대리물질·2단계 노출기준·3모델 비교)

v1 파이프라인은 이성질체 혼합물·석유계 혼합물(UVCB)이 하나라도 들어 있으면 제품 전체를 보류했다(270건 중 216건).
또한 모든 물성·노출기준이 "미승인"이어서 결과가 0건이었다.
v2는 **모든 입력값에 출처를 붙인 채** 270건 전부를 계산한다. 승인 대기 대신, 가정이 들어간 부분은 민감도 분석으로 다룬다.

## 실행 순서

> **용도 분류 감사(2026-10-24):** [USE_CLASSIFICATION_AUDIT.md](USE_CLASSIFICATION_AUDIT.md), [270건 감사표](data/use_classification_audit_270.csv). 분석은 용도가 확인되지 않았거나 충돌하는 제품을 사전 제외하고, 세부 사용군·계절형으로 비교합니다. 기존 `results/core/`의 성능 수치는 감사표를 적용해 다시 분석하기 전 **과거 결과**입니다.


```bash
cd research/korean_thinners_271/v2
python src/parse_korean_oel.py      # 고시 JSON -> data/korean_oel_table.csv (659 물질)
python src/build_properties.py      # 63 화학종 SMILES·MW·증기압(25/35°C) + 출처
python src/build_oel.py             # Tier A(국내 공식) / Tier B(보충) 노출기준 배정
python src/run_v2.py --samples 100 --seed 20261024 --hanna-home <HANNA 경로>
python src/analyze_v2.py            # RQ1–RQ5 표·그림·headline_numbers.json
python src/build_flammability.py    # 성분별 인화점·폭발하한 63종 (출처·추정 표시)
python src/extract_msds_flash.py <MSDS PDF 폴더들>   # MSDS 제9항 인화점·폭발범위 -> data/msds_flash.csv (pdftotext 필요)
python src/run_flammability.py --samples 100 --seed 20261024 --hanna-home <HANNA 경로>   # PYTHONPATH=src 필요
python src/analyze_flammability.py  # RQ6–RQ8, results/rq/headline_flammability.json (결과 해석: PAPER_DIRECTION_fire_health.md)
python src/analyze_pooled.py     # 본문 결과 P1(전체 방향성)·P2(같은 비교군 제품쌍)·P3(사업장용 목록), PYTHONPATH=src
python src/analyze_core.py       # 본문(core) 표·그림: results/core/ (PYTHONPATH=src). 범위는 CORE_PAPER.md, 확장은 FOLLOW_UP_PAPERS.md
python -m unittest discover -s src -p "test_v2.py"
```

필요 패키지: `thermo chemicals rdkit numpy pandas matplotlib`, HANNA 실행 시 `torch transformers`.
HANNA는 공개 다성분 버전(https://github.com/marco-hoffmann/HANNA, 가중치·ChemBERTa 포함)을 그대로 사용한다. 새로 학습하지 않는다.
Windows에서는 `--hanna-home D:\fire2\HANNA`처럼 지정한다.

## v1 대비 바뀐 점

| 항목 | v1 | v2 |
|---|---|---|
| 계산 대상 | 미승인으로 0건 (UVCB 포함 216건은 구조적 보류) | **270건 전수** |
| 크실렌 혼합 이성질체 (181건) | 차단 | o/m/p 대리 조성 + 구간 샘플링 |
| 솔벤트 나프타 등 UVCB (8종) | 차단 | 공개 근거 기반 대리 성분 구간 (`data/surrogates.csv`) |
| 고분자 (1건) | 차단 | 비휘발 성분으로 증기 계산에서 제외 (몰수 무시 가능) |
| 물성 | 수동 승인 대기 | `chemicals/thermo` 문헌 상관식 (방법명 기록) + QC |
| 노출기준 | 수동 승인 대기 | 고시 별표1 자동 파싱 (Tier A) + 표시된 보충값 (Tier B) |
| 기준 없는 성분 | 지수에서 조용히 제외 | **OEL 커버리지 공백으로 별도 정량** |
| 활동도계수 모델 | Raoult, HANNA | Raoult, **수정 UNIFAC(Dortmund)**, HANNA |
| 온도 | 25°C | 25°C, 35°C |

## 핵심 정의

- 성분 분압 `p_i = x_i γ_i P_i^sat(T)`, 총압 기준 `P = 101.325 kPa`
- 혼합물 증기위험지수 `I_mix = Σ_i 10^6 p_i / (P · OEL_i)` — 포화 상부공간 농도를 노출기준 배수로 합산한 값 (Popendorf 증기위험비의 혼합물 확장, 고시의 혼합물 가산식과 같은 형태)
- OEL 커버리지 `= Σ_{OEL 있는 i} p_i / Σ_i p_i` — 포화증기 중 노출기준이 있는 성분의 비율
- Tier A: 고용노동부 고시 별표1의 TWA만 사용. CAS 일치 또는 고시의 이성체 일괄 항목(크실렌, 트리메틸벤젠, 헥산 이성체)
- Tier B: Tier A + 출처가 표시된 보충값. PGMEA는 EU IOELV 50 ppm. 나프타 대리 성분에는 국내 유사물질 값(TMB 25, 나프탈렌 10, 헵탄 400, 스토다드 100 ppm)
- 스토다드 용제(고시 No.326)는 UVCB 자체에 공식 TWA가 있으므로 그 대리 성분 전부에 100 ppm을 적용

## 해석의 한계 (논문에 그대로 명시할 것)

1. 결과는 **평형 포화 상부공간의 노출 잠재력**이다. 실제 호흡영역 농도·개인노출·법 위반 판정이 아니다.
2. MSDS 제3항 함량은 범위이고 실제 배합과 다를 수 있다(윤충식 외 2000). 이 불확실성은 구간 샘플링(RQ3)으로만 다룬다.
3. UVCB 대리 조성은 가정이다(`surrogates.csv`의 basis 열 참조). 하한·상한 구간을 샘플링에 포함한다.
4. 증기압은 문헌 상관식이다. 상관식이 없는 5종(PGMEA, 젖산에틸, 숙신산디메틸, MHIB, 초산 t-부틸)은 단일 문헌값과 Clausius–Clapeyron 추정이다.
5. 법제처 2026-10-24 현행 행정규칙 목록상 노출기준은 고용노동부고시 제2020-48호(2020.1.16 시행)이다. 저장소의 고시 JSON도 2020년판이지만 **별표 CAS별 OEL 매핑과 UVCB·이성체 예외는 행 단위 확인을 완료한 것이 아니다**.
6. 같은 이름·같은 조성 범위의 문서는 `조성지문`으로 묶어 하나만 순위 비교에 쓴다.

## 용도검증을 반영한 논문용 분석

본래 `제품271건_용도재분류.csv`는 출처를 보존하고, `data/use_classification_audit_270.csv`를 `analyze_v2.load()`에서 문서키별 1:1 병합합니다. 감사표 누락, ID·원본 SHA256 불일치 시 분석을 중단합니다. `탐색비교허용=0`인 제품은 비교군 통계·탐색 후보 선정에서 제외합니다. `현장대체검증=0`은 해당 제품이 실제 대체 가능함을 승인하지 않았다는 뜻입니다.

```bash
cd research/korean_thinners_271/v2
python -m unittest discover -s src -p "test_use_audit.py" -v
python src/analyze_v2.py
PYTHONPATH=src python src/analyze_flammability.py
PYTHONPATH=src python src/analyze_core.py
```

기존 `results/core/` 출력은 재실행 전 결과이며 **최신 감사표가 자동으로 기존 CSV를 수정하지 않습니다**. 원본 MSDS 제1항의 권고용도 추출 필드가 비어 있는 경우, 실제 PDF/제조사 TDS 검수 없이 새로운 확정 용도를 생성하지 않습니다.

