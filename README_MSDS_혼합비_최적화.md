# 국내 MSDS 정보공백 탐색용 혼합비 최적화 (1차 구현)

이 분석은 국내 MSDS 성분 기재 한계농도에 **근접하지만 미만인 미량 성분**의 평형 증기상 상대 지수가 원래 제품의 주성분 혼합비에 따라 어떻게 달라지는지 탐색합니다. `HANNA` 또는 `UNIFAC` 등의 활동도계수를 활용해 최적화할 수 있으며, **AI를 새로 학습시키지 않고 기존 사전학습 HANNA 예측기를 사용**합니다.

> **경고:** 모든 미량 성분은 원래 EPA ChemExpo 제품에서 실측된 성분이 아니라 **가상 첨가**입니다. 최적화된 혼합물도 실제 유통제품이 아닙니다. 예측 평형 증기 조성과 국내 TWA 가중 상대 지수는 작업장 실측 노출·인체 위해도·물질별 독성가산·법적 MSDS 기재 생략의 적법성을 뜻하지 않습니다. 해당 한계농도 미만이라도 실제 기재 여부는 제11조·별표 6·시행규칙 별표 18 및 다른 법적 의무를 별도로 검토해야 합니다.

## 파일과 기존 코드 연결

- `optimize_msds_gap.py`: 최적화 실행 파일
- `test_optimize_msds_gap.py`: 질량수지·중복·누락·제약조건·라울의 법칙 예비검사
- `compare_thinner_raoult_unifac_hanna.py`: 이미 사용한 라울의 법칙·UNIFAC·HANNA 계산기
- `analyze_thinner_multichem.py`: 기존 다물질 등록정보와 구성 관리
- `apply_korean_msds_oel.py`: 국내 노출기준 JSON의 별표 1을 읽는 검증된 CAS별 TWA 파서
- `multichem_screening_registry.csv`: 벤젠 0.1 wt%를 **예비 유해성 기준**으로만 등록한 CSV. 메탄올·아세톤은 법적 분류 미확정이므로 숫자 한계농도를 자동 부여하지 않음
- `law_data/행정규칙/화학물질 및 물리적 인자의 노출기준.json`: 고용노동부고시 2020-48 (다운로드본 메타데이터)
- `law_data/행정규칙/화학물질의 분류·표시 및 물질안전보건자료에 관한 기준.json`: 고용노동부고시 2026-26 (다운로드본 메타데이터)

GitHub 법령 JSON에 저장된 시행일과 고시번호는 실행마다 결과 메타데이터에 기록합니다. 국내 OEL 고시가 현재 최신본인지 여부는 사용 전 국가법령정보센터에서 재검토해야 합니다. 위 OEL 고시는 2020-48 다운로드본입니다.

## 실행 방법 (Windows PowerShell)

```powershell
cd D:\fire
python -m pip install scipy
python -m unittest -v test_optimize_msds_gap.py

# 1) 가장 단순한 알고리즘 검증: 라울의 법칙으로 배합비 최적화
python optimize_msds_gap.py --optimize-model raoult --check-models raoult

# 2) 핵심 연구: AI(HANNA) 목적함수 최적화 후 다른 모델과 비교
python optimize_msds_gap.py --optimize-model hanna --check-models raoult unifac hanna --strict-models

# 3) 단일 제품 + 좁은 조건으로 실행
python optimize_msds_gap.py --product-ids 1183746 --doses 0.099 --optimize-model hanna --check-models raoult unifac hanna
```

**주의:** 이미 설치한 HANNA의 공식 코드와 학습 파라미터는 기본적으로 `D:\fire\HANNA` 폴더에서 불러옵니다. 이 모델이 없거나 실행 실패하면 자동으로 라울의 법칙 결과를 대체하지 않고 **실패를 명시**합니다. 분석 온도는 순수 증기압 자료가 준비된 **25°C에 고정**되어 있습니다. `thermo`가 없는 경우 `python -m pip install thermo`로 설치할 수 있습니다.

### 최적화 조건 (기본 설정)

- EPA 제품: 세 제품 (`1208563`, `273605`, `1183746`)
- 대상 물질: 벤젠 CAS `71-43-2` (원래 제품에 없다면 별도 가상 첨가)
- 벤젠 액상 함량: `0.039`, `0.099 wt%` (두 값 모두 0.1% 미만)
- 목적함수: 대상 성분의 **국내 TWA 가중 상대지수(%) 최대화**
- 제품의 원래 주성분을 상대적으로 **±50% 범위** 안에서만 조정하고, 원래 주요 성분 간 **총 재배분량 ≤10%p**로 제한합니다. 예: 원래 20 wt%라면 허용 범위 10~30 wt%지만, 합계 보존과 재배분 제한도 함께 충족해야 합니다.
- 미량 성분 함량 고정, 총 액상 조성 100 wt%, 모든 성분 양수
- 결정적 초기화 시드 42, 다중 초기점 SLSQP. **전역 최적해를 증명하지 않습니다.**
- 처음 보고된 EPA 제품 배합과 최적화 배합을 각각 라울·UNIFAC·HANNA에 넣어 **모델 간 민감도**를 점검할 수 있습니다.

```powershell
# 연구자가 설정한 탐색 범위 변경 예시
python optimize_msds_gap.py --max-reallocation-pp 5 --relative-change 0.25 --starts 3 --maxiter 60
```

메탄올·아세톤 등 아직 물질별 법적 기재 한계농도가 검증되지 않은 경우에는 `--target-cas 67-56-1`만으로 실행하면 **명시적으로 오류**가 발생합니다. 정말로 비법적 탐색 시나리오만 수행하고 싶다면 아래처럼 별도 가정값을 지정하십시오. 이 결과는 MSDS 법적 평가 결과로 해석하면 안 됩니다.

```powershell
python optimize_msds_gap.py --target-cas 67-56-1 --exploratory-cutoff-wt-pct 1 --doses 0.039 0.099 --optimize-model raoult --check-models raoult
```

## 출력 파일

폴더: `D:\fire\chemexpo_paint_thinner_output\msds_optimized_screening_output`

| 파일 | 내용 |
|---|---|
| `01_scenario_status.csv` | 시나리오별 실행·실패·한계농도 불일치 및 사유 |
| `02_optimization_summary.csv` | 제품별 최적화 전후 국내 가중 기여율, 증기상 함량, 실제 배합 재분배량 |
| `03_baseline_vs_optimized_compositions.csv` | 성분별 원래 EPA 함량, 가상 첨가 후 함량, 최적화 함량 |
| `04_cross_model_checks.csv` | 최적화 모델 이외의 다른 두 모델에서 전후 결과 재계산 |
| `05_search_audit.csv` | SLSQP 실행 성공 여부와 선택한 최적 후보의 기여율 |
| `06_run_metadata.json` | 사용한 법규번호·시행일·알고리즘·제약조건·한계 및 실행 상태 |

## 학술적 검증이 필요한 사항

1. **실제 GHS 분류 및 성분별 기재 기준 확인:** 벤젠 이외 물질은 KOSHA 및 국내 법규와 매칭하여 CAS별 적용 조항을 검증해야 합니다. 0.1 wt%를 모든 물질에 적용하지 않습니다.
2. **실측 기액평형(VLE) 검증:** 독립 VLE 문헌자료에서 HANNA/UNIFAC 예측 오차를 확인해야 합니다. 세 모델의 합의는 실측 정확성을 보증하지 않습니다.
3. **제품 조성 근거 검증:** ChemExpo 제품 보고 함량의 출처와 신뢰도를 독립적으로 재검토해야 합니다.
4. **탐색 범위 민감도:** 현재의 ±50% 개별 변경·10%p 총 재배분은 **연구자가 정한 탐색 제약**이지 제조사 허용오차나 실험 사실이 아닙니다. 5·10·15%p 등으로 교차검증해 결과의 안정성을 제시해야 합니다.
5. **정보공백의 법적 결론 분리:** `0.099 wt%`라는 수치상 한계농도 미만은 실제 특정 성분 기재 생략이 적법함을 뜻하지 않습니다. 실제 한국 MSDS를 분석하거나 관계 조항을 완전히 검토하기 전에는 위반/부적절 판정 금지.
6. **통계 독립성:** 한 제품에서 혼합비만 바꾼 최적화 시나리오는 독립적인 실제 제품 표본으로 취급하지 않습니다.
