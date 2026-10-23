# 국내 시판 신나 MSDS · 기액평형 · AI 배합비 연구

현재 연구 대상은 **사용자가 확보한 삼화 87건·노루 184건, 총 271개 MSDS PDF 문서**입니다. 이는 271개의 서로 다른 제품이나 국내 시장 대표표본을 의미하지 않습니다.

**[연구 데이터, 실행법 및 분석 한계](research/korean_thinners_271/README.md)**

## 저장소 구조

| 경로 | 역할 |
| --- | --- |
| `research/korean_thinners_271/data/` | 271개 MSDS 문서, 제3항 1,462개 성분기록, 62종 CAS, 14개 잠정 제품군, 원문 식별정보 |
| `research/korean_thinners_271/src/` | 조성 범위 검증, 25°C Raoult 예비계산, 검증조건을 요구하는 혼합비 탐색, 테스트 |
| `research/korean_thinners_271/docs/` | 후속 최적화에 필요한 검증자료·삭제 정리 기록 |
| `domestic_product_inputs/` | 기존 KCC 포함 국내 MSDS 7건/성분 39건 **별도 참고자료** |
| `audit_korean_msds_data.py`, `test_audit_korean_msds_data.py` | 위 7건 원문 추출자료 검증 |
| `thinner_vle_reference_25C.csv` | 7종 순물질의 25°C 참고 물성, 반드시 출처·단위 검증 |
| `kosha_msds/chem_ids.csv`, `kosha_msds/section09.csv` | 순물질 물성 관련 참고 수집표 |
| `law_data/행정규칙/화학물질 및 물리적 인자의 노출기준.json` | **과거 고용노동부 노출기준 스냅샷**. 현행 2026 고시 확인 전 사용 금지 |

## 실행

저장소 루트에서 실행합니다. 기본 데이터 검증 및 Raoult 예비계산에는 Python 표준 라이브러리만 필요합니다.

```powershell
python research/korean_thinners_271/src/validate_271.py
python research/korean_thinners_271/src/run_raoult_271.py --samples 2000
python -m unittest discover -s research/korean_thinners_271/src -p "test_*.py" -v
python -m unittest test_audit_korean_msds_data -v
python research/korean_thinners_271/src/optimize_verified_blend.py --write-template constraints_to_verify.json
```

## 현 단계 한계

- **Raoult 결과 6개 MSDS 문서만 계산**했습니다. 다성분 HANNA 모델/학습 체크포인트는 **미연동** 상태입니다.
- 평형증기압은 실제 시간당 증발률이나 작업장 공기 중 노출농도와 다릅니다.
- 조성 범위가 공개된 MSDS만으로 실제 제조사 배합비를 확정하지 않습니다.
- 자일렌 이성질체, 나프타 등 UVCB는 임의의 하나의 SMILES로 치환하지 않습니다.
- AI 최적 배합이나 제품군별 최악 순위는 아직 도출하지 않았습니다. 성분 구조·물성·국내 TWA·공정 호환성 검증이 필요합니다.

## 이전 연구파일 정리

미국 EPA ChemExpo 데이터·벤젠 가상첨가 예비실험·사용하지 않는 법령/서식 데이터 및 임시 산출물을 **연구 브랜치에서 삭제**했습니다. 복원 기준 커밋, 유지 파일 및 정리 범위는 [저장소 정리 기록](research/korean_thinners_271/docs/저장소_정리_기록.md)에 남겨두었습니다.
