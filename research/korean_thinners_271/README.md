# 국내 시판 신나 271개 MSDS 문서: 기액평형·증기상 잠재지표·동일 용도 비교·배합 탐색

## 연구 목적과 현재 도달 단계

1. **기액평형 계산**: 25°C Raoult 기준모델과 설치된 **외부 HANNA 2026 다성분 예측 모델**의 활동도계수를 같은 MSDS 조성 시나리오에 적용한다. HANNA는 새로 학습한 모델이 아니며, 계산 정확도는 독립 실측 검증이 필요하다.
2. **증기상 유해성 잠재성 분석**: `p_i = x_i γ_i P_i^sat`로 계산된 분압을 기록한다. 국내 현행 TWA가 *CAS별로 모두* 원문 대조 승인된 경우에만 `max_i [10^6 p_i / (101.325 × TWA_i)]`를 **평형분압 기반 스크리닝 지표**로 사용한다. 가산 독성지수, 작업장 실제 노출지수, 법적 기준초과값은 아니다.
3. **동일 용도 내 후보 비교**: MSDS 제1항 용도, 같은 공정·수지 호환성, 중복/개정판 검토, 원문 근거를 **제품별 승인**하고 비교군 내 두 문서 이상이 있을 때만 지표상 참고순위를 산출한다. 이는 실제 대체 가능 판정이나 구매 권고가 아니다.
4. **AI 기반 가상 배합비 탐색**: 기존 `optimize_verified_blend.py`의 *검증된 공정 제약 입력*이 모두 충족되면 HANNA 어댑터로 예측한 활동도계수를 이용한 Pareto 가상 조성 탐색을 수행한다. 최적 배합의 물리적 타당성 검증은 별도다.

**분석 범위**: 삼화 87개, 노루 184개, 합계 **271개 PDF 문서** (독립 제품 271종이 아님); 1,462 성분기록, 62 CAS, 14 잠정 분류군. 현재 25°C Raoult/HANNA 완비 참고물성이 있는 문서는 **6개**이며, 그 외는 `271_제품별_분석가능성_감사.csv`에서 보류 처리한다. 제품명·용도 분류 오류, 중복·개정판, 자일렌 이성질체/나프타 UVCB가 남아 있다.

## Windows / chemprop 실행

로컬에서 `fire` 저장소 파일이 `D:\fire2`에 있고 HANNA checkout이 `D:\fire2\HANNA`에 있는 경우:

```powershell
conda activate chemprop
cd D:\fire2
$env:HANNA_HOME = 'D:\fire2\HANNA'

# 1. 데이터 진단 및 기존 Raoult 비교
python research/korean_thinners_271/src/validate_271.py
python research/korean_thinners_271/src/run_raoult_271.py --samples 2000

# 2. 실제 HANNA 예측 (PyTorch, ChemBERTa, 체크포인트 등은 로컬 HANNA 설치 필요)
python research/korean_thinners_271/src/run_hanna_271.py --samples 20 --batch-size 8

# 3. 271 문서 감사 + 설명적 증기압 표 + 아직 미승인인 검토서식 두 개 생성
python research/korean_thinners_271/src/screen_same_use_271.py --write-templates

# 4. 전문가 검토와 현행 국내 고시 확인 후에만 조건부 참고순위 산출
#    먼저 output/동일용도_검증입력_서식.csv 및 output/국내TWA_공식검증입력_서식.csv 를
#    복사하여 각 CAS/제품 검증 근거를 확인하고 '승인'으로 수정한다.
python research/korean_thinners_271/src/screen_same_use_271.py `
  --review-csv "D:\fire2\research\korean_thinners_271\output\검증완료_비교군.csv" `
  --verified-oel-csv "D:\fire2\research\korean_thinners_271\output\검증완료_한국TWA.csv"

# 5. 혼합비 탐색은 모든 검증조건을 근거문서로 채운 후에만 가능
python research/korean_thinners_271/src/optimize_verified_blend.py --write-template constraints_to_verify.json
# 아래 명령은 constraints_to_verify.json의 모든 조건을 검증한 후에만 실행
python research/korean_thinners_271/src/optimize_verified_blend.py `
  --config constraints_to_verify.json `
  --gamma-predictor hanna_gamma_adapter:predict_gammas --n 1000
```

Windows에서 저장소 폴더가 `D:\fire2\fire-main`처럼 다른 경로인 경우 **그 위치에서 실행**해야 한다. `HANNA_HOME`에는 실제 HANNA checkout 위치를 지정한다. Python 코드를 HANNA checkout으로 옮길 필요는 없다.

### 검증파일 규칙

- `동일용도_검증입력_서식.csv`: `제품ID`, `비교군ID`(같은 수지·공정·요구 성능일 때 동일), `용도검증`, `공정호환성검증`, `중복개정검토`, `근거문서`. 3개 검증 필드 모두 정확히 `승인`, 근거문서 비어 있지 않아야 한다. 제품명과 MSDS 제1항이 불일치한 문서는 별도 원문 정정 전 차단한다.
- `국내TWA_공식검증입력_서식.csv`: `CAS`, `한국_TWA_ppm`, `기준시점`, `공식출처`, `검증상태`. 마지막 필드가 `승인`이고 값·출처·시점이 있어야 한다. NIOSH/국외 OEL을 한국 TWA로 채우지 않는다. 출처의 유효성은 프로그램이 자동 판단하지 않으며 승인자의 확인 책임이다.
- 일부 성분이라도 미검증이면 **해당 문서 전체 유해성 잠재지표 계산을 보류**한다. 불명 성분을 0으로 처리하지 않는다.
- 용도분류가 불명확하거나 비신나 오탐인 문서는 자동 평가 제외. **동일 브랜드/제품명도 무조건 같은 제품으로 간주하지 않는다.**
- 모든 승인 데이터 파일은 별도 로컬 작업파일로 보관하고 원문 출처를 기록한다. GitHub에 검증되지 않은 실제 배합·TWA·등급을 게시하지 않는다.

### 파일별 역할

| 파일 | 역할 |
|---|---|
| `src/hanna_gamma_adapter.py` | 로컬 HANNA의 `ln γ`를 `γ`로 변환; 모델 실패 시 즉시 중단, Raoult 대체 금지 |
| `src/run_hanna_271.py` | 동일 가능한 MSDS 조성에서 두 VLE 모델 예측, 분압 및 실행 입력 SHA-256 기록 |
| `src/screen_same_use_271.py` | 271 문서 커버리지, 비교군 전문가 승인, 성분별 검증 TWA 기반 조건부 지표 및 참고순위 |
| `src/optimize_verified_blend.py` | 기능·호환성·물성·한국 TWA가 확인된 가상 조성의 제약조건부 Pareto 탐색 |
| `src/test_same_use_271.py` | 미검증 자료 순위 차단, TWA 누락 차단, 원본 해시 차단, 성분 합계 검사 테스트 |

### 신규 산출물 (`research/korean_thinners_271/output/`)

- `271_제품별_분석가능성_감사.csv`: 분류 신뢰도, 모델 출처, 동일 용도·중복·공정 검증, 보류 이유를 271개 문서별로 기록
- `제품군별_증기상_기술비교.csv`: 예측 가능한 문서의 Raoult/HANNA 압력값, **유해성/구매 순위 아님**
- `동일용도_저유해성_후보선별.csv`: 조건이 모두 승인된 문서의 **평형분압 기반 잠재지표 참고순위**; 기본적으로 행 0개가 정상
- `동일용도_제품쌍_잠재지표비교.csv`: 참고순위에 따른 중앙값 비율과 샘플링 최솟값/최댓값 범위의 중첩 여부
- `동일용도_시나리오_잠재지표.csv`: 성분별 최대 TWA 정규화 평형분압과 기여 CAS

## 해석·검증 한계

- 원본 MSDS는 **함량 범위**만 제공. 시나리오 5–95백분위는 임의 생성된 유효 조성들의 기술적 분포이지, 실제 제조 배합이나 신뢰구간이 아니다. 현 조성 생성기는 전역 조성 공간의 균일 표본추출을 보장하지 않는다.
- 평형분압 기반 TWA 정규화는 작업장 공기노출을 의미하지 않는다. 1보다 크더라도 노출기준 위반이 아니다. 여러 물질의 독성 가산성도 가정하지 않으며 최대 성분 지표만 사용한다.
- 유해성 잠재지표의 참고순위가 낮더라도 더 안전한 제품 또는 사용 가능한 대체제품으로 인증하지 않는다. 실측 VLE/헤드스페이스, 원료 동정, 온도별 순물성, 실제 공정 배출·노출, 대체재의 인화성·용해력·도막품질 검증이 필요하다.
- 기존 7종 물성표는 출처·온도·단위를 논문 분석 전 재검증해야 하고, 외부 HANNA 모델의 도료신나 적용 가능 영역도 별도 검증해야 한다.

참고: `docs/최적화_입력조건.md`, `data/모델검수_가이드.md`.
## 논문 타당성 검증 실험: 다중 시드·조성 추출 민감도·외부 VLE

**[전체 실험 설계 및 로컬 실행 가이드](docs/논문_가치검증_전체실험.md)**

기존 200시나리오 단일실험을 확장한 별도 파이프라인입니다. 두 조성추출법(legacy/Gibbs)과 여러 난수시드로 동일 MSDS 조성 범위 안의 HANNA/Raoult 차이를 평가합니다. 결과 파일의 질량합계, 몰분율, 분압, 활동도계수, 입력 SHA-256을 검사하고 검증 실패 시 중단합니다.

~~~powershell
cd D:\fire2
$env:HANNA_HOME = 'D:\fire2\HANNA'
python research/korean_thinners_271/src/paper_experiment_271.py run --samples 200 --batch-size 8
python research/korean_thinners_271/src/plot_paper_experiment_271.py --reports research/korean_thinners_271/output/paper_full/reports
python research/korean_thinners_271/src/validate_external_vle_271.py --write-templates research/korean_thinners_271/output/external_vle_evidence
~~~

**주의:** 실제 로컬 HANNA 체크포인트가 있어야 하며, 2개 추출법 × 3개 시드 × 현재 6개 문서 × 200조성으로 총 7,200개 가상 시나리오를 계산합니다. 독립 실측과 순물질 물성 승인 전에는 모델 정확도·산업안전·저유해성 제품 순위를 주장하지 않습니다. 기존 동일 용도 검증과 최적화 안전조건은 유지합니다.
