# 국내 MSDS 기반 도료 희석제 증기상 평가 및 대체후보 스크리닝

이 프로그램은 **국내 유통이 확인된 도료 희석제**의 원래 배합을 변경하지 않고, 공급사의 MSDS 성분정보를 25°C HANNA 평형 증기상 예측과 **국내 TWA (고용노동부고시 제2020-48호)**에 연결합니다. 국내 MSDS 작성 근거는 **고용노동부고시 제2026-26호**입니다.

**연구용 도구**입니다. 실제 작업자 흡입노출, 작업환경측정값, 법적 노출기준 초과, MSDS 기재 의무 위반, 특정 제품의 최종 대체 적합성을 판정하지 않습니다. 추가 실측과 도료/공정 호환성 확인이 필요합니다.

## 적용 기준

- 제품은 `market_country=KR`로 확인되고, `product_evidence_url`, `msds_revision`, `source_type`을 입력한 **국내 시판제품만** 분석 대상에 포함합니다. 온라인 URL의 진위는 연구자가 직접 확인해야 합니다.
- 한국 법제처 법규 JSON `law_data/행정규칙/화학물질 및 물리적 인자의 노출기준.json`에서 실제 별표 1(`별표키=000100`)의 **TWA ppm**만 사용합니다. STEL, mg/m³ 및 미국 NIOSH 수치로 대체하지 않습니다.
- `thinner_vle_reference_25C.csv`의 분자량, SMILES, 25°C 증기압은 물성 참고자료로만 사용합니다. 해당 CSV의 `OEL_NIOSH_ppm`은 **무시**합니다.
- `wt_pct_min` 및 `wt_pct_max` 값이 동일하고, 전체 성분의 함량 합이 100±0.1 wt%, 모든 CAS에 물성 및 국내 TWA가 있는 제품만 수치 예측합니다. **성분이 생략되거나 함량이 구간으로 표시된 MSDS**는 임의로 채우지 않고 `EXCLUDED`로 기록합니다.
- `use_group`은 제품의 사용 목적입니다. **대체후보 비교**는 HANNA 계산 성공, `compatibility_group` 일치, `compatibility_verified=yes`, `compatibility_evidence_url`이 모두 있고 **동일 그룹에 2개 이상 제품이 있을 때만** 출력합니다. 단순 명칭이 유사한 제품은 대체 가능한 것으로 가정하지 않습니다.
- 가격, 환기량, 작업시간, 공간크기는 입력하지 않아도 됩니다. 미국 ChemExpo는 이 프로그램의 **국내 제품 순위 입력값으로 사용하지 않습니다.**

## 입력 템플릿

`domestic_product_inputs/products.csv` 에 제품별 성분을 행 단위로 추가합니다. 헤더만 제공하므로, **실제 국내 제품을 아직 수집하지 않았다면 결과는 빈 목록**입니다. 가상의 국내 시판제품을 만들지 않습니다.

필수 컬럼: `product_id, product_name, manufacturer, market_country, product_evidence_url, msds_revision, source_type, use_group, compatibility_group, compatibility_verified, compatibility_evidence_url, CAS, component_name, wt_pct_min, wt_pct_max, msds_disclosure`.

`source_type`: `manufacturer_msds`, `supplier_msds`, `measured_verified_product` 중 하나. `msds_disclosure`: `yes/no/unknown` 기록용. **MSDS에 없는 성분은 HANNA만으로 발견할 수 없습니다.** 누락 검증은 제품과 동일성이 확보된 독립 실측 자료가 있는 경우에만 가능합니다.

## 실행

```powershell
cd D:\fire
python -m unittest -v test_screen_korean_paint_thinners.py
python screen_korean_paint_thinners.py --input domestic_product_inputs\products.csv --models hanna unifac raoult --hanna-repo D:\fire\HANNA
```

가용성이 없는 HANNA나 UNIFAC이 오류를 발생하면 해당 모델 상태를 `MODEL_FAILED`로 남깁니다. **라울 계산값을 AI 결과로 위장하지 않습니다.**

결과: `domestic_screening_output/01_domestic_data_audit.csv`, `02_korean_twa_crosswalk.csv`, `03_component_vapor_screen.csv`, `04_product_screen_summary.csv`, `05_eligible_selection_candidates.csv`, `06_run_metadata.json`.

## 지표 해석

`equilibrium_screening_priority = (predicted_partial_pressure_kPa / 101.325 × 1,000,000) / Korean_TWA_ppm`입니다. 이는 **평형 증기 분압을 표준기압의 ppm 상당량으로 환산하여 TWA로 나눈 우선검토 지표**이지, 작업장 농도/TWA의 비율이 아닙니다. 제품별 비교는 최대 개별성분 지표를 사용해 미확인 독성 가산성을 가정하지 않습니다. 여러 성분의 합계는 부가 설명용이지 혼합물 노출기준 초과를 뜻하지 않습니다.

한계: 25°C 증기압 고정; 제품 내 수지·고형분·비공개 물질을 무시할 수 없음; 평형값은 실제 증발속도나 실시간 환기·분무 공정 노출을 나타내지 않음; 제품 호환성 근거 없이는 대체 권고 금지. 데이터가 부족하면 **추천 대신 보류 및 부족 사유**를 표시하는 것이 올바른 결과입니다.

## 국내 규정 출처

- [국가법령정보센터: 화학물질의 분류·표시 및 물질안전보건자료에 관한 기준, 2026-26](https://www.law.go.kr/admRulLsInfoP.do?admRulSeq=2100000278320)
- [국가법령정보센터: 화학물질 및 물리적 인자의 노출기준, 2020-48](https://www.law.go.kr/admRulSc.do?menuId=5&query=%ED%99%94%ED%95%99%EB%AC%BC%EC%A7%88&subMenuId=41&tabMenuId=183)

코드가 읽는 법규 JSON 버전은 결과 `06_run_metadata.json`에 보존됩니다. 개정 시 자료를 갱신하고 근거를 다시 점검해야 합니다.