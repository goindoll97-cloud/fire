> **2026-10-24 논문 본실험 확정: [RQ1 동일 용도 배수 · RQ2 Raoult/HANNA 순위·상대등급 역전 · RQ3 MSDS 하중상 불확실성](research/korean_thinners_271/docs/PAPER_RQ1_RQ2_RQ3_270.md).** 270개 MSDS 전체가 연구 모집단입니다. 오래된 271건 파일럿 안내 대신 이 문서와 `src/run_full_paper_experiment.py`, `src/run_paper_research_questions.py`를 우선 사용합니다. **법정 TWA·물성·TDS 미승인 항목을 임의 승인하지 않으며, 혼합물 지표는 성분별 합계 평형 잠재치입니다(작업장 실측 노출이 아님).**


> **권장 연구 실행 방법:** [원본 신나 MSDS 270건 번호별 00–09 전체 파이프라인](research/korean_thinners_271/docs/PIPELINE_00_09_270.md). 원문 검증·54 CAS 물성 승인·국내 TWA·동일 용도 공정 호환성 조건을 단계별로 관리합니다. `.env`의 실제 LAW_OC/KOSHA_MSDS_SERVICE_KEY는 GitHub에 공개하지 않습니다.

# 국내 시판 신나 MSDS · 기액평형 · 저유해성 후보 탐색

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
