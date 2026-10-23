# 국내 시판 신나 연구: 원본 MSDS 270건 기반 00–09 재현 파이프라인

## 무엇을 검증하는가

**연구 목표**: 국내 시판 신나 MSDS 구성성분 함량 범위와 25°C 기액평형 AI(HANNA)를 결합해 비이상적 증기상 조성·분압을 계산하고, **실제로 같은 도료·수지·공정에서 호환되는 신나 후보 안에서 유해성 잠재지표가 상대적으로 낮은 제품**을 탐색한다. 검증된 기능 제약조건이 확보된 경우에만 새로운 배합의 Pareto 가상 탐색을 별도로 수행한다.

**입력 기초 데이터**: 노루 183 + 삼화 87 = 270개의 *실제 신나* MSDS PDF; 원본 제3항 구성성분 기록 1,454행, 54종 CAS. 초기 검색 결과 271건 가운데 NR007은 ‘신나무’가 포함된 수성페인트 오탐으로 제외됐다. 파일명·동명제품·개정 이력상 270개의 *독립 시판 제품*으로 확정한 것은 아니다.

이 파이프라인은 기존 7 CAS 파일럿 물성표나 `제품271건_용도재분류.csv`를 사용해 6건만 선택하지 않는다. **검증된 54 CAS 중 입력 가능한 모든 PDF**에 HANNA를 적용하되, 미확인 혼합물·불완전 조성은 임의 대체하지 않는다. MSDS/기존 모델의 물성 정확성은 별도 실측 검증이 필요하다.

## 설치 및 비밀키 관리

저장소 `D:\fire2`에서 실행한다고 가정한다.

```powershell
conda activate chemprop
cd D:\fire2
git status --short
git pull --ff-only origin main
python -m pip install numpy pymupdf
```

기존 `D:\fire2\.env` 파일의 `LAW_OC`, `KOSHA_MSDS_SERVICE_KEY`는 그대로 보존된다. 이 파일을 GitHub에 업로드하지 마라. 저장소에는 빈 **`.env.example`**만 공개한다. `00`은 `.env`가 없을 때만 빈 예제파일을 복사한다.

```powershell
python research/korean_thinners_271/src/00_environment.py
# LAW_OC/KOSHA_MSDS_SERVICE_KEY=설정여부 True/False 만 출력, 값은 출력하지 않음
```

`HANNA_HOME=D:\fire2\HANNA` 경로를 실제 설치 위치로 지정한다. 현재 PowerShell에 같은 이름의 환경변수가 이미 설정되어 있다면 그 값이 우선된다.

## 00–09 순서별 실행

각 단계는 `research/korean_thinners_271/output/full270/` 아래 별도 폴더에 결과를 보관하며, 원본·검토 CSV를 덮어쓰지 않는다.

```powershell
# 00. 로컬 비밀키 환경 설정 검증 (GitHub에 실제 .env 업로드 안 함)
python research/korean_thinners_271/src/00_environment.py

# 01. 원본 ZIP 270 PDF 재추출 + SHA/CAS/함량 구간 대조 (권장)
python research/korean_thinners_271/src/01_verify_raw_msds.py `
  --noru-zip "D:\fire2\노루_신나(1).zip" `
  --samhwa-zip "D:\fire2\삼화_신나(1).zip"

# ZIP이 없어도 저장소에 커밋된 원본-대조 정제 DB만 검증할 수 있음
# (이 경우 SOURCE_PDF_SHA_VERIFIED가 아닌 SNAPSHOT_ONLY로 기록됨)
# python research/korean_thinners_271/src/01_verify_raw_msds.py

# 02. 270개 MSDS 제품검토, 54 CAS 물성/국내 TWA 검증서식 만들기
python research/korean_thinners_271/src/02_make_review_templates.py

# 03. LAW_OC로 법제처 행정규칙·시행법령 XML 내려받기, 원문 검수 대기
python research/korean_thinners_271/src/03_fetch_official_sources.py
# API 인증/IP 등록이 아직 완료되지 않았다면 네트워크만 건너뛰기
# python research/korean_thinners_271/src/03_fetch_official_sources.py --offline

# 04. 아래의 3개 검증서식을 실제 근거에 따라 채운 후 재평가
python research/korean_thinners_271/src/04_check_model_eligibility.py

# 05. 모든 성분 구조·물성·MSDS 함량이 승인된 PDF만 실제 HANNA 예측
python research/korean_thinners_271/src/05_predict_hanna_raoult.py `
  --samples 200 --batch-size 8 --sampler gibbs --seed 20261023

# 06. 공인 국내 TWA가 모든 성분에 승인된 문서만 잠재지표 계산
python research/korean_thinners_271/src/06_score_headspace_hazard.py

# 07. 같은 수지·공정에서 호환성이 승인된 독립 제품 2개 이상만 비교
python research/korean_thinners_271/src/07_compare_same_use.py

# 08. 최적화 제약 템플릿 생성 (미승인 자동 최적화 금지)
python research/korean_thinners_271/src/08_prepare_blend_optimization.py

# 09. 커버리지/보류/동일용도 비교 산출물을 연구요약으로 통합
python research/korean_thinners_271/src/09_integrated_paper_report.py
```

한 번에 순서대로 실행하려면 아래도 가능하다. 단, 승인되지 않은 단계는 **0건 산출이 정상**이고, 학술 검증이 완료됐다는 뜻은 아니다.

```powershell
python research/korean_thinners_271/src/run_full_pipeline.py --offline `
  --noru-zip "D:\fire2\노루_신나(1).zip" `
  --samhwa-zip "D:\fire2\삼화_신나(1).zip"
```

## 02 단계 후 사람이 검토해야 할 필드

`output/full270/02_review/`의 세 CSV는 처음 모두 `미검증/미승인`이다. 다시 `02`를 실행해도 **이미 작성한 승인사항을 덮어쓰지 않는다**.

| 검증서식 | 근거가 필요한 핵심 항목 | 자동 허용 여부 |
|---|---|---|
| `제품_원문TDS_검증.csv` | MSDS 제1항 용도, 제3항 함량, 구성성분의 모델상 완결성, MSDS 개정/중복판 검토, 독립제품, 실제 수지·공정 호환성, TDS 성능범위, 출처 | 자동 승인하지 않음 |
| `CAS_물성_검증.csv` | CAS 구조·SMILES, 순물질 분자량, **25°C 증기압**, 단위/온도·공식 출처, 혼합이성질체/UVCB 제외 | 자동 승인하지 않음 |
| `CAS_한국TWA_검증.csv` | 대한민국 현행 고시의 **CAS별 TWA(ppm)**, 기준일, 고시명·번호·공식 원문·별표 대조 | 자동 승인하지 않음 |

**제품별 ‘전체모델성분확정’**은 MSDS의 표시성분만으로 유효한 완전 혼합물 모델을 구성할 근거가 있을 때만 승인한다. MSDS에 함량 구간이 나온다고 실제 제조비율이 제공된 것은 아니며 비공개 성분이나 복합물은 자동 무시할 수 없다.

`04` 이후 검수 CSV가 수정되면 기존 예측을 재사용하지 않도록 해시 불일치를 검사한다. `04 → 05 → 06 → 07`을 재실행하라.

## API 키와 공식 정보의 역할

- 법제처 `LAW_OC`: 국가법령정보 공동활용 `eflaw`, `admrul` 목록의 공식 응답 XML·SHA-256을 로컬 저장한다. 이것만으로 CAS별 국내 TWA가 자동 승인되는 것은 아니다. 공식 별표·시행일·고시번호를 연구자가 직접 대조해야 한다.
- `KOSHA_MSDS_SERVICE_KEY`: [공공데이터포털 KOSHA MSDS 조회](https://www.data.go.kr/data/15157612/openapi.do) 사용. 공식 기본 URL `https://apis.data.go.kr/B552468/msdschem1/getChemList001`. 서비스별 검색 매개변수는 최신 API 명세에서 확인 후 03 단계에 `--kosha-query-param <확인된_파라미터명>`을 추가한다(기본 CAS 3종만 시험, `--kosha-limit 1..54`). 확인되지 않은 엔드포인트를 추측하거나 인증키를 로그에 출력하지 않는다. **KOSHA의 참고용 화학물질 MSDS는 제조사 원본 MSDS 및 법정 고시 TWA의 자동 대체물이 아니다.**
- 법제처 OpenAPI 연결이 실패한다면 인증/승인·IP 허용 여부를 확인한다. 오류 메세지에는 키를 포함하지 않는다.

## 06: 유해성 잠재지표와 07: 같은 용도 안에서 덜 위험한 신나 고르기

검증된 성분별 예측 평형분압 `p_i`에 대해 두 모델을 같은 조성 시나리오에 적용한다. 성분별 고시 TWA가 확인된 경우에만 아래 **스크리닝 지표**를 정의한다.

`H_s = max_i [ 10^6 p_{i,s} / (101.325 kPa × TWA_i(ppm)) ]`

- `07`은 **MSDS 제1항 실제 용도**, 같은 수지·공정의 검증된 비교군 ID, 독립제품 및 동명/개정판 처리, TDS 성능 제약에 승인 표시와 근거문서가 있는 **제품 2개 이상**에 대해서만 표를 만든다.
- 최솟값과 최댓값이 겹치는 제품쌍은 `가상조성범위분리=아니오`로 명시한다. 서로 다른 제품의 제조사 정확배합·분포를 알 수 없으므로 중앙값 순위의 통계적 안전성 확정은 아니다.
- **잠재지표가 낮은 제품 = 안전한 제품/작업장 노출이 낮은 제품/인화성·기능이 동등한 대체품**이라는 주장을 하지 않는다. 이는 논문 내 엄격히 제한된 후보 탐색 결과이다.
- 동일용도 비교가 불가능한 문서도 반드시 `비교보류_270건.csv`에 보류 사유를 남긴다.

## 08: 검증된 제약 아래에서만 Pareto 탐색

`08_optimization/검증전_Pareto_제약서식.json`은 비승인 양식이며 자동으로 실행되지 않는다. 실제 도료 용해력/호환성·도막품질·건조성·인화점, 물성·TWA 근거를 승인한 뒤에만 기존 `optimize_verified_blend.py --config <근거 JSON> --gamma-predictor hanna_gamma_adapter:predict_gammas`로 별도 실행한다.

## 산출물 및 논문 타당성 기준

| 출력 폴더 | 주요 출력 | 실증 수준 |
|---|---|---|
| `01_source` | ZIP 원본 SHA, 제3항 재추출 비교, **로컬 전용 원문 발췌** | 데이터 계보 |
| `02_review` | 원문·성분/물성·TWA 근거와 미검증 상태 | 수동 검증 절차 |
| `03_official` | 기관 공식 응답 XML, 증거 SHA | 법령·MSDS 증거 수집 |
| `04_eligibility` | 270건 중 계산 가능한 문서와 제외 사유 | 모델 적용 범위 |
| `05_predictions` | 실제 HANNA/Raoult 쌍대 분압·기상몰분율 | *예측*, 실측 정확도 아님 |
| `06_hazard` | 성분별 최대 TWA 정규화 평형분압 | 실제 노출지수 아님 |
| `07_same_use` | 같은 용도 낮은 잠재지표 후보, 제품쌍 중첩, 보류 270건 | 대체 후보 참고만 가능 |
| `08_optimization` | Pareto 제약 입력서식, 보류 조건 | 최적 실제품 아님 |
| `09_report` | 연구 진행 표·보류 사유·논문 해석범위 | 논문 준비상태 |

**중요**: 외부 실측 VLE 또는 headspace 참값이 없는 한 모델 예측의 실제 정확도가 검증된 것으로 주장하면 안 된다. 270개 MSDS는 통계적으로 독립된 임의 시장 표본이 아니다. 반복 시뮬레이션은 독립 제품 샘플 수를 늘리지 않는다. 모든 단계에서 실제 작업장 공기측정/직업노출기준 위반 여부/산업현장 대체 가능성은 판정하지 않는다.