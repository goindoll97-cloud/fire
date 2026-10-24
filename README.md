# 사업주의 화재 이익과 근로자의 흡입 보호는 같은 방향인가

> **2026-10-24 용도 분류 감사 완료:** MSDS 270건을 문서키·원본 SHA256 기준으로 대조하여 **5건 분류 보완, 108건 탐색비교 유보, 162건 조건부 탐색비교**로 분리했습니다. 용도검수 결과와 판단 근거는 [USE_CLASSIFICATION_AUDIT.md](research/korean_thinners_271/v2/USE_CLASSIFICATION_AUDIT.md), 전체 270행은 [use_classification_audit_270.csv](research/korean_thinners_271/v2/data/use_classification_audit_270.csv)에 있습니다. 아래 결과 요약은 감사 후 재실행한 값입니다.

### 국내 시판 신나 270종의 MSDS 조성 기반 화재·흡입 동시 평가 (AI 활동도계수 HANNA)

한국산업보건학회지 투고를 목표로 하는 연구 저장소입니다.

## 연구 질문
소규모 사업장은 근로자 건강 보호만으로는 움직이기 어렵습니다. 반면 **화재 손실과 위험물 허가 부담은 사업주 본인의 비용**입니다.
그래서 같은 용도의 신나 안에서 **사업주가 화재 때문에 고를 제품(인화점이 높은 제품)이 근로자의 흡입 노출도 낮추는지** 확인하고, 두 위험이 엇갈리는 제품을 찾아냅니다.

## 방법
1. 원본 신나 MSDS 270건(삼화 87·노루 183)의 제3항 조성에서 출발합니다. 석유계 혼합물은 대리 성분 구간으로 전개합니다.
2. **HANNA(기계학습 활동도계수)** 한 모델로 증기 분압을 구하고, 같은 계산에서 두 지표를 냅니다.
   - 흡입 노출 지수 `I_mix = Σ 10^6 p_i/(P·OEL_i)`
   - 혼합물 인화점(Liaw 모델)과 폭발하한 도달비
3. 같은 용도·같은 계절/건조 등급 안에서 비교합니다. 노출기준이 없는 성분이 증기의 대부분인 제품은 순위에서 빼고 따로 보고합니다.
4. AI의 역할은 기존 방법(Raoult, 수정 UNIFAC)이 계산하지 못하거나 결과가 달라지는 제품을 다루는 것입니다. 인화점 예측 정확도를 높인다고 주장하지 않습니다.

## 현재 결과 (비실험 스크리닝, HANNA, 25°C, 제품당 100개 조성 시나리오, 용도감사 반영)
- **일반 방향성:** 용도를 구분하지 않은 195종에서 인화점이 높을수록 노출 지수가 낮습니다(Spearman −0.85, 95% CI −0.88~−0.80). 제조사별로 따로 봐도 같습니다.
- **같은 용도 안에서도:** 감사를 통과한 같은 비교군의 117개 제품쌍 중 79.5%에서 두 지표가 같은 방향입니다(비교군 단위 클러스터 부트스트랩 95% CI 67.2~88.1%)(두 지표가 모두 확실히 분리되는 쌍에서는 90.5%). 일장일단인 군 내 제품은 3종입니다.
- **사업장용 목록:** 비교군별로 엄격 Pareto 비지배 29종과 동시 개선 상대가 있는 39종을 구분했습니다(`P3_field_list.csv`). 건조속도와 수지 호환성이 검증되지 않아 대체 승인은 0종입니다.
- 노출기준 공백 제품 10종("LV", "Zero VOC" 포함)은 25°C에서 포화증기가 폭발하한에 도달합니다(상한은 확인하지 않음). 저VOC 표시가 저화재나 저독성을 뜻하지 않습니다.
- 사업주가 보는 MSDS 인화점으로 비교하면 관계가 더 약합니다(전체 −0.64, 제품쌍 66%). MSDS 인화점과 조성의 불일치는 원문·배합·실측 재검토 신호이지 오기재의 증거는 아닙니다.

## 반드시 알아둘 한계
결과는 포화 상부공간의 **노출 잠재력**이지 실제 작업장 노출이 아닙니다. 혼합제품 인화점 실측은 본 연구에서 수행하지 않으며, MSDS 값 비교와 계산모델 간 일치는 실험 검증이 아닙니다. 법제처 2026-10-24 현행 행정규칙 목록은 노출기준을 2020-48호(2020-01-16 시행)로 표기합니다. 저장소의 별표 CAS별 TWA 매핑은 독립 재대조가 남아 있고, 21°C·70°C의 모형상 석유류 경계는 법적 위험물 분류 확정이 아닙니다.

## 어디부터 읽을까

- **[비실험 연구 보완 결과](research/korean_thinners_271/v2/RESULTS_NONEXPERIMENTAL_20261024.md)** — 의존적 제품쌍 CI(67.2~88.1%), Pareto 정정, 3모델·TWA·최신 MSDS 선택 민감도.
- **[비실험 논문 설계 및 한계](research/korean_thinners_271/v2/NONEXPERIMENTAL_MANUSCRIPT_READINESS.md)** — 실험 없이 가능한 검증과 해석 한계.

- **[용도 분류 감사](research/korean_thinners_271/v2/USE_CLASSIFICATION_AUDIT.md)** — 270건 제1항 대조, 분류 보완·충돌·108건 비교 유보, 감사 기준 재실행 방법.
- **[노벨티 및 검증 설계](research/korean_thinners_271/v2/MANUSCRIPT_NOVELTY_AND_VALIDATION.md)** — 선행연구 차별성, 가설, 심사 대응과 실측 검증 계획.
| 문서 | 내용 |
|---|---|
| [CORE_PAPER.md](research/korean_thinners_271/v2/CORE_PAPER.md) | 본문 범위, 결과 표 C1~C4, 정직하게 써야 할 한계 |
| [FOLLOW_UP_PAPERS.md](research/korean_thinners_271/v2/FOLLOW_UP_PAPERS.md) | 후속 논문 후보 8개 (AI 배합 탐색, MSDS 신뢰성, 노출기준 공백 등) |
| [v2/README.md](research/korean_thinners_271/v2/README.md) | 재현 명령, 입력 자료와 가정 |
| [PAPER_DIRECTION_fire_health.md](research/korean_thinners_271/v2/PAPER_DIRECTION_fire_health.md) | 화재·흡입 동시 평가의 근거와 남은 일 |
| [PAPER_DIRECTION_v2.md](research/korean_thinners_271/v2/PAPER_DIRECTION_v2.md) | 노출기준 공백 중심의 이전 방향 (확장 주제로 이동) |

결과 파일은 `research/korean_thinners_271/v2/results/core/`(본문용)와 `results/rq/`(확장 분석)에 있습니다.

## 재현
```bash
cd research/korean_thinners_271/v2
python src/parse_korean_oel.py && python src/build_properties.py && python src/build_oel.py
python src/build_flammability.py
python src/run_v2.py --samples 100 --seed 20261024 --hanna-home <HANNA 경로>
PYTHONPATH=src python src/run_flammability.py --samples 100 --seed 20261024 --hanna-home <HANNA 경로>
PYTHONPATH=src python src/analyze_core.py
python -m unittest discover -s src -p "test_v2.py"
```
HANNA는 공개 다성분 버전(https://github.com/marco-hoffmann/HANNA)을 그대로 사용하며 새로 학습하지 않습니다. MSDS PDF에서 인화점을 뽑는 방법은 v2/README.md를 보세요.

## 선행 연구와의 관계
혼합 용제 증기의 활동도계수 기반 예측(Hori 등 2012, Tischer·Roitzsch 2022), 혼합물 인화점 예측(Liaw 계열, 하동명·이성진 2014), 용제 선정 지침의 건강·안전 점수(CHEM21), 국내 MSDS 오기재(윤충식 외 2000, 안병준 외 2018)가 각각 있습니다. 제품 MSDS 조성에서 출발해 화재와 흡입을 같은 계산으로 함께 평가하고 사업주 이익과 연결한 연구는 검색에서 찾지 못했습니다. 투고 전 재확인이 필요합니다.

---

## 이전 자료 (v1·271건 파일럿)
아래 설명과 문서는 이전 단계의 것입니다. 현재 분석은 위 v2를 기준으로 합니다.

> **2026-10-24 v2:** [v2/README.md](research/korean_thinners_271/v2/README.md) · 결과 `research/korean_thinners_271/v2/results/rq/`
> **2026-10-24 RQ1~RQ3:** [PAPER_RQ1_RQ2_RQ3_270.md](research/korean_thinners_271/docs/PAPER_RQ1_RQ2_RQ3_270.md)
> **원본 MSDS 270건 00–09 파이프라인:** [PIPELINE_00_09_270.md](research/korean_thinners_271/docs/PIPELINE_00_09_270.md). `.env`의 LAW_OC/KOSHA_MSDS_SERVICE_KEY는 GitHub에 공개하지 않습니다.

### v1 파이프라인 설명 (참고용)

기존 검색 수집 271개 PDF 문서(삼화 87·노루 184)를 원본 재검증한 결과, 노루 자료의 `NR007`은 비신나 수성페인트 오탐이었습니다. **실제 신나 MSDS는 270건(삼화 87·노루 183)**이며 제3항에 1,454개 성분기록과 54종 CAS가 포함되어 있습니다. 초기 271건 기록은 추적용으로 유지합니다. 이는 고유 제품 270종이나 국내시장 대표표본을 의미하지 않습니다.

**[전체 연구 실행법과 검증 기준](research/korean_thinners_271/README.md)**

## 4개 연구 모듈

1. **기액평형 예측**: Raoult 기준모델과 사용자 PC의 다성분 HANNA 2026 예측을 동일 MSDS 조성에 적용 (`src/run_raoult_271.py`, `src/run_hanna_271.py`). HANNA는 신규 학습이 아닌 기존 모델 사용.
2. **증기상 유해성 잠재성**: 성분별 평형분압 및 공식 검증된 국내 TWA 기반 잠재지표. 이 지표는 실측 공기중 농도·직업노출·법규 위반 여부를 뜻하지 않음.
3. **동일 용도 후보 선별**: 제품군 검토 및 비교 가능성·중복/개정·공정 호환성 **승인 후에만** 참고순위 계산 (`src/screen_same_use_271.py`). 검증되지 않은 271개 문서에 순위를 부여하지 않음.
4. **AI 가상 배합비 탐색**: 다섯 검증조건을 명시적으로 충족한 경우만 HANNA 활동도계수를 이용한 제약조건부 Pareto 탐색 (`src/optimize_verified_blend.py`). 성능·안전성을 보장하는 최적제품 판정이 아님.

## 빠른 시작 (Windows / conda chemprop)

```powershell
cd D:\fire2
$env:HANNA_HOME = 'D:\fire2\HANNA'
python research/korean_thinners_271/src/validate_271.py
python research/korean_thinners_271/src/run_hanna_271.py --samples 20 --batch-size 8
python research/korean_thinners_271/src/screen_same_use_271.py --write-templates
python -m unittest discover -s research/korean_thinners_271/src -p 'test_*.py' -v
```

`HANNA_HOME`은 외부 HANNA 코드와 가중치가 설치된 로컬 폴더입니다. GitHub에는 가중치를 포함하지 않습니다. `screen_same_use_271.py`는 공식 출처·전문가 확인 없이 유해성 순위를 출력하지 않습니다.

## 자료·범위·주의

- `research/korean_thinners_271/data/`: 271 MSDS 문서 추적정보, 1,462 성분행, 62 CAS, 14 잠정 용도군
- `research/korean_thinners_271/src/`: 검증, Raoult, HANNA 연결, 비교/보류 감사, 검증조건부 가상 배합 탐색과 테스트
- `research/korean_thinners_271/output/`: 실행시 생성되는 분석결과 (Git 미추적)
- `thinner_vle_reference_25C.csv`: 7개 순물질 25°C **참고 물성(출처 검증 필요)**
- `domestic_product_inputs/`: 기존 KCC 등 **별도 참고자료 7건**, 271건에 자동 합산하지 않음

현재 검증 필요 없는 것은 **인터페이스 작동**뿐입니다. 실제 모델 예측 정확도, 현행 한국 TWA의 CAS별 확인, 271개 문서 중 구조 미상 성분과 중복 문서, 공정별 호환성은 후속 검증이 필요합니다. 실제 기화속도·작업장 노출농도·저유해성 제품 인증을 결과로 주장하지 않습니다.
## 논문 가치 검증을 위한 재현 실험

기존 MSDS 271문서 분석을 바탕으로 동일 조성 내 Raoult/HANNA 차이, 표본추출법·난수시드 민감도, 증기상 성분 전환, 외부 VLE 검증을 단계적으로 평가합니다. **모형 간 불일치만으로 모델 정확도 향상이나 안전성 우위를 주장하지 않습니다.**

- [실험 설계·Windows 실행법](research/korean_thinners_271/docs/논문_가치검증_전체실험.md)
- [기존 데이터·동일 용도 비교와 배합 탐색 가이드](research/korean_thinners_271/README.md)

## 원본 PDF 기반 전체 연구 DB

[270건 원본대조 및 1,454개 성분·54 CAS 연구 DB](research/korean_thinners_271/docs/원본MSDS_270건_데이터재구축.md) — 원본 문서의 해시·함량 범위·복합물 구조 차단·용도 및 중복 검토·물성 검증상태를 분리 관리합니다. 결측 물성, 공식 국내 TWA, 제조사 실제 배합비는 자동 생성하지 않습니다.

### 54 CAS·270개 제품 자동 근거수집 (02a–02d)

[PubChem/법제처/KOSHA 및 원본 MSDS 자동 수집 안내](research/korean_thinners_271/docs/PIPELINE_02a_02d_AUTO_ENRICH.md) — CAS 구조·분자량은 출처와 함께 미승인 상태로 자동 채우고, 제품 용도 원문 및 한국 TWA 공식 근거를 검토 대기 파일로 모읍니다. **현행 한국 TWA·물성·동일 용도 비교제품을 자동 승인하지 않습니다.**
