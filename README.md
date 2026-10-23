# fire — 국내 시판 도료 희석제의 기액평형·유해성 연구

이 저장소는 **이전 ChemExpo/EPA 기반 시범연구**와 **국내 제조사 MSDS에 기반한 연구**를 구분해 관리합니다.

## 이번 연구: 삼화·노루 MSDS 271건
- **[연구 데이터·실행방법](research/korean_thinners_271/README.md)**
- 제품 271개 **PDF 문서**, 성분 1,462행, 제품군 14종, 고유 CAS 62종
- 2026 다성분 HANNA 분석 후보: 46건(구조 기반 1차 필터), 기존 물성으로 Raoult 예비계산 가능: 6건
- 정확 배합비·실제 증발률·실측 작업장 노출·AI 최적 배합비는 **현재 확보되지 않음**
- 원본 MSDS PDF는 업로드하지 않음(로컬 보관, 파일명과 SHA-256으로 추적)

```powershell
python research/korean_thinners_271/src/validate_271.py
python research/korean_thinners_271/src/run_raoult_271.py --samples 2000
python -m unittest discover -s research/korean_thinners_271/src -p "test_*.py"
```

## 이전 코드 보존
- `screen_korean_paint_thinners.py`: 정확 조성·국내 TWA·동일 목적 확인 조건의 엄격 선별 v1
- `compare_thinner_raoult_unifac_hanna.py`: 이전 모델 비교 실험(국내 271건 결과 아님)
- `chemexpo_paint_thinner_output/`: 미국 ChemExpo 자료, **국내 시판제품 근거로 사용하지 않음**
- `domestic_product_inputs/`: 국내 MSDS 7건 원문 추출 기반 기존 시험 입력

**주의:** 기액평형 증기압은 시간당 실제 증발속도나 직업노출 농도를 의미하지 않습니다.
