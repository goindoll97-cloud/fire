# 270건 신나 연구 데이터 자동 수집: 02a–02d

> **대상**: 삼화·노루 신나 원본 270건(1,454개 성분행, 54 CAS). 기존 00–09 순차 파이프라인 **02 다음**에 수행하는 근거 수집 단계입니다. 
>
> **중요**: 자동 수집과 사람의 최종 검증은 구분합니다. PubChem 자료는 후보로 입력되며, 국내 법정 TWA·제조사 공정 호환성·안전성 순위를 자동 승인하지 않습니다.

## 첫 실행(Windows PowerShell)

실제 GitHub 저장소가 `F:\fire2`에 있고 `(chemprop)`이 활성화되어 있다면:

```powershell
cd F:\fire2
git status --short
git pull --ff-only origin main
python -m pip install pymupdf

# 기존 00–02 결과가 이미 있으므로 01/02를 다시 할 필요는 없습니다.
# 02a: PubChem CAS, SMILES, MolecularWeight, InChIKey 및 직접 25°C VP 원문 후보
python research/korean_thinners_271/src/02a_auto_prefill_pubchem.py

# 02b: LAW_OC 행정규칙 검색/본문/별표목록과 KOSHA_MSDS_SERVICE_KEY CAS 검색 근거
python research/korean_thinners_271/src/02b_auto_prefill_korean_twa.py

# 02c: 업로드했던 원본 PDF ZIP 두 개의 제1항 용도 증거를 자동 발췌
python research/korean_thinners_271/src/02c_auto_extract_msds_use.py `
  --noru-zip 'F:\fire2\노루_신나(1).zip' `
  --samhwa-zip 'F:\fire2\삼화_신나(1).zip'

# 02d: 처리된 내용 및 아직 원문 검수/승인해야 하는 항목을 1개 CSV로 정리
python research/korean_thinners_271/src/02d_review_queue.py
```

`git status --short`에 수정사항이 있거나 pull이 실패하면 강제로 초기화하지 마세요. 기존 `02_review/*.csv`의 수동 기입 내용은 보존합니다.

`02a`는 54개 CAS를 순서대로 처리합니다. 혼합 자일렌, UVCB 나프타 등 **단일 분자로 나타낼 수 없는 복합 CAS**는 구조 자동 조회를 차단합니다. PubChem이 CAS 동의어와 일치하는 **단일 CID**를 반환한 경우에만 빈 `SMILES`, `분자량_g_mol`, `구조공식출처`를 채웁니다. PUG View의 명시적 25°C 원문에서 유효한 압력 단위를 하나로 해석할 수 있는 경우에만 `25C_포화증기압_kPa`를 **미승인 후보값**으로 채웁니다. 조건이 20°C인 값은 25°C로 임의 변환하지 않습니다. 실험값/추정값과 원 출처의 타당성을 별도 확인해야 합니다.

`02b`는 `.env`의 `LAW_OC`, `KOSHA_MSDS_SERVICE_KEY`로 법제처 **현행 행정규칙 검색**, 일련번호가 유일하면 **본문 조회**, 그리고 **별표서식 목록**까지 수집하고 응답 SHA-256을 저장합니다. KOSHA MSDS는 `searchWrd=CAS`, `searchCnd=1`로 조회하며 참고자료입니다. 법제처 별표 목록과 본문 XML이 받아져도 실제 고시의 시행일·별표표·CAS별 TWA(ppm) 값을 정확하게 대조하기 전에는 법정 TWA 승인값을 입력하지 않습니다. 즉 **TWA 데이터 54행이 전부 자동 완성되는 기능은 아닙니다.**

`02c`는 원본 ZIP 파일을 다시 읽어 문서별 SHA-256과 일치하는 첫 1~2쪽에서 MSDS 제1항의 용도 관련 문맥을 추출합니다. TDS나 다른 도료/수지와의 공정 호환성을 추정해 승인하지 않습니다.

## 결과물

모든 출력은 Git에서 제외된 `research/korean_thinners_271/output/full270/` 아래에 저장됩니다.

| 출력 폴더 | 파일/내용 |
|---|---|
| `02_review/CAS_물성_검증.csv` | PubChem의 후보 SMILES·분자량과 직접 25°C VP 후보만 빈 칸에 입력. 승인 필드는 변경하지 않음 |
| `02a_pubchem/` | PubChem CID·InChIKey·구조·분자량·기상압력 원문후보 및 상태표 |
| `02b_twa_evidence/` | KOSHA 참고 MSDS CAS 검색 원문, 법제처 행정규칙 본문/별표목록 XML, 미검증 근거표 |
| `02c_product_evidence/` | 제품 270건의 원본 MSDS 제1항 사용용도 발췌 |
| `02d_review_queue/` | 54 CAS 물성, 54 CAS TWA, 270 제품 총 378개 검토행과 상세 보류 사유 |

한 번 조회해 저장된 파일에 연구자가 직접 기입한 검증 상태를 **다음 실행에서 임의로 덮어쓰지 않습니다.** PubChem 후보표는 동일 CAS별 갱신되므로 수동 검증기록은 `02_review/` 원본의 검증 필드에 기입하세요.

## 자동화 후 계속 실행

```powershell
# 실제 공식 출처와 원본 MSDS/TDS를 대조하여 검증할 항목만 검증 필드에서 '승인'
python research/korean_thinners_271/src/04_check_model_eligibility.py
python research/korean_thinners_271/src/05_predict_hanna_raoult.py --samples 200 --batch-size 8 --sampler gibbs
python research/korean_thinners_271/src/06_score_headspace_hazard.py
python research/korean_thinners_271/src/07_compare_same_use.py
python research/korean_thinners_271/src/08_prepare_blend_optimization.py
python research/korean_thinners_271/src/09_integrated_paper_report.py
```

`HANNA` 계산 및 **같은 용도 안에서 낮은 유해성 잠재지표를 갖는 신나 비교**는 기존의 모든 확인 조건(순물질·물성·성분 완결성·현행 국내 TWA·동일 용도/수지/공정/TDS·독립 제품 판정)이 충족된 경우에만 진행됩니다. 0건이 나와도 검증 누락을 조용히 0으로 취급하거나 안전한 제품으로 오인하지 않습니다.

## 간단 시험

```powershell
# 첫 실행 시 통신·API 키 문제가 없는지 CAS 3개만 빠르게 테스트
python research/korean_thinners_271/src/02a_auto_prefill_pubchem.py --limit 3
python research/korean_thinners_271/src/02b_auto_prefill_korean_twa.py --limit 3

# 인터넷 없이 API 파일·승인 가드 유닛테스트
python -m unittest discover -s research/korean_thinners_271/src -p 'test_auto_enrich_270.py' -v
```

- 공식 API 참고 문서: [PubChem PUG REST](https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest), [PUG View](https://pubchem.ncbi.nlm.nih.gov/docs/pug-view), [법제처 행정규칙 본문 API](https://open.law.go.kr/LSO/openApi/guideResult.do?htmlName=admrulInfoGuide), [행정규칙 별표 API](https://open.law.go.kr/LSO/openApi/guideResult.do?htmlName=admrulBylListGuide), [KOSHA MSDS 조회 서비스](https://www.data.go.kr/data/15157612/openapi.do).
- API 키는 실제 `.env`에만 존재하고 Git에 커밋되지 않습니다. 공식 XML은 로컬 생성 폴더에만 두고 인증키 값은 기록하지 않습니다.

## KOSHA HTTP 403과 로컬 실행기 누락 복구 (2026-10-24)

KOSHA 조회 54건이 모두 `REQUEST_FAILED_...HTTP_403`이라면
**검색된 CAS가 0건이 아니라 기관 API 게이트웨이 접근이 거부된 상황**입니다.
공공데이터포털 [물질안전보건자료 조회 서비스](https://www.data.go.kr/data/15157612/openapi.do)
(서비스 ID 15157612)의 *활용신청·승인 상태*와 발급받은 인증키를 확인해야 합니다.
`serviceKey`의 디코딩키/인코딩키 처리, 제공기관 권한 및 요청 IP 차단도 점검합니다.
인증키 값이나 전체 호출 URL을 채팅/GitHub/이슈에 붙여넣지 마세요.

```powershell
cd F:\fire2
git fetch origin main
$src = "research/korean_thinners_271/src"
# 아래 경로에 작성한 로컬 코드가 없다면 GitHub 승인본 복구
git restore --source=origin/main --worktree -- "$src/run_auto_enrichment.py" "$src/auto_enrich_270.py" "$src/02b_diagnose_kosha_access.py"
Test-Path "$src/run_auto_enrichment.py"
python "$src/02b_diagnose_kosha_access.py"
```

`02b_diagnose_kosha_access.py`는 벤젠 CAS 한 번만 요청합니다.
결과는 `output/full270/02b_twa_evidence/KOSHA_단일요청_접근진단.json`에
인증키 제외하여 저장하며 `http_status`, `provider_code`,
`provider_status`, `diagnosis`을 구분합니다.
포털에 등록된 키를 `.env`에 URL 인코딩된 형태로 저장한 경우,
이 진단기는 디코딩 후 요청 시 한 번만 URL 인코딩하도록 처리합니다.

공식 반환 오류코드 예:
- `20`: 키 누락 또는 서비스 이용 권한 거부
- `30`: 미등록 인증키
- `31`: 인증키 만료
- `22`/`23`: 호출량/호출속도 제한
- `29`: IP 차단

`403`만으로 구체적인 원인은 확정하지 않습니다. XML 기관 오류코드가
존재하면 함께 확인합니다. 또한 일괄 `02b`는 첫 403 이후 나머지 CAS 요청을
`SKIPPED_AFTER_HTTP_403_CHECK_API_PERMISSION`으로 표시해 불필요한 호출을
중단합니다. 데이터 없음 `NO_EXACT_CAS_IN_RESPONSE`과 인증 거부를 구분합니다.

KOSHA MSDS 조회는 **참고용**이고, 현행 한국 법정 TWA는
법제처 공식 고시 본문·별표에서 CAS·단위·시행일 확인이 별도로 필요합니다.


## 파일 누락 반복을 막는 로컬 Git 동기화 점검

사용자 PC가 `F:\\fire2`이고 GitHub `main`에 있는 Python 파일이 로컬에
누락되는 사례가 반복되었습니다. `git fetch`는 작업 폴더를 갱신하지 않습니다.
다음 명령은 **GitHub origin/main에 등록된 연구 Python 소스 파일 중 로컬에 없는
것만 복구**하며, **존재하는 로컬 수정 파일이나 .env, 실험 결과는 변경하지 않습니다.**

```powershell
cd F:\fire2
git fetch origin main

$src = "research/korean_thinners_271/src"
$missing = @(
  git ls-tree -r --name-only origin/main -- "$src" |
  Where-Object { -not (Test-Path -LiteralPath $_) }
)

Write-Host "누락 파일 수: $($missing.Count)"
$missing | ForEach-Object { Write-Host "복구 예정: $_" }
if ($missing.Count -gt 0) {
  git restore --source=origin/main --worktree -- $missing
}
```

동기화의 근본 원인은 다음 명령 출력으로 확인합니다.
```powershell
git rev-parse --show-toplevel
git branch --show-current
git status --short
git rev-parse --short HEAD
git rev-parse --short origin/main
git sparse-checkout list
```
`git sparse-checkout list`는 설정되지 않은 경우 오류 메시지를 출력할 수 있습니다.
이는 그 자체로 결함이 아닙니다. 로컬과 원격의 HEAD가 서로 다르다면 정상적으로
`git pull --ff-only origin main`이 가능한지 확인합니다. 로컬 수정 때문에
병합이 거부될 경우 데이터를 강제 초기화하지 않습니다.

단일 API 접근 진단기 `src/02b_diagnose_kosha_access.py`는 이제
`pipeline_270_core.py`를 import하지 않고 저장소 루트의 `.env`를 직접
안전하게 읽으므로 **독립적으로 실행**할 수 있습니다.
