도료 희석제 증기 조성 예측: 라울의 법칙 / UNIFAC / HANNA 비교
=============================================================

1. 포함된 파일
-------------------
compare_thinner_raoult_unifac_hanna.py    비교 분석 실행 코드
thinner_vle_reference_25C.csv           7개 물질의 25°C 증기압, 분자량, NIOSH 노출기준, SMILES
README_도료희석제_3모델비교.txt           본 안내서

2. 입력 파일 경로 (이미 코드에 적용)
--------------------------------
D:\fire\chemexpo_paint_thinner_output\04_model_input_candidates_REVIEW.csv

실제 EPA 제품 3종(Data Document ID): 1208563, 273605, 1183746

출력 폴더(자동 생성):
D:\fire\chemexpo_paint_thinner_output\thinner_vle_comparison_output\

3. 필요한 프로그램 및 설치
-------------------------
아래 명령은 Windows PowerShell에서 실행합니다.

A. 라울 + UNIFAC만 먼저 비교할 때 (Python 환경에 설치):
    python -m pip install numpy pandas thermo

    python compare_thinner_raoult_unifac_hanna.py --models raoult unifac

B. HANNA까지 비교할 때: 공식 학습 모델이 필요합니다.
   HANNA는 단순 pip install만으로 사용하지 않으며, 작성자가 공개한
   모델 파라미터와 ChemBERTa 의존성을 갖춘 환경이 필요합니다.

   (1) Git 설치 + Miniconda/Anaconda 설치 후:
    git clone https://github.com/marco-hoffmann/HANNA.git D:\fire\HANNA

   (2) HANNA 공식 환경 만들기 (conda 또는 Anaconda Prompt):
    conda env create -f D:\fire\HANNA\environment.yml
    conda activate HANNA_environment
    python -m pip install thermo

   (3) 비교 코드와 이 CSV 참조표를 같은 폴더에 두고:
    python compare_thinner_raoult_unifac_hanna.py --models raoult unifac hanna

   HANNA 코드는 공식 `utils.HANNA_predictor.HANNA_Predictor`의 `.predict`
   메서드(성분 SMILES, 몰분율 2차원 배열, 온도 Kelvin)를 직접 사용합니다.
   이 메서드는 ln(gamma)을 반환하므로 exp(ln(gamma))로 변환합니다.
   공식 훈련된 앙상블(`models/HANNA/ensemble/`)과 스케일러가 필요합니다.
   최초 실행 때 대용량 모델 파일 또는 사전학습 언어모델 가중치를 내려받아야
   할 수 있어 인터넷 연결이 필요합니다.

   설치하지 않은 모델은 UNAVAILABLE로 기록되며, 절대로 라울의 법칙
   gamma=1을 UNIFAC/HANNA 예측값처럼 대신 기록하지 않습니다.

   꼭 세 모델 모두 정상 실행되었는지 강제 확인하려면:
    python compare_thinner_raoult_unifac_hanna.py --strict-models

   선택사항: Dortmund 수정 UNIFAC도 같이 계산:
    python compare_thinner_raoult_unifac_hanna.py --models raoult unifac dortmund hanna

4. 출력파일
-----------
01_model_status.csv             모델별 성공/실패와 오류메시지
02_component_predictions.csv   개별성분 wt%, x, gamma, 포화증기압, y, 상대 지수
03_benzene_comparison_long.csv  벤젠 농도별 모델별 분석
04_benzene_comparison_wide.csv  제품/벤젠함량 행별로 모델 비교표
05_differences_vs_raoult.csv   라울 기준 모델간 차이
06_across_model_ranges.csv      성공한 모델간 최대최소 범위
07_run_metadata.json           모델 구성, 자료원, 연구상 한계

5. 예비분석 시나리오와 해석
------------------------
(a) EPA 미국 도료 희석제 3종의 공개문서상 'reported' 조성을 사용합니다.
    reported 값은 독립 GC 실측값을 의미하지 않습니다.
(b) 한국 2022년 논문에서 다른 제품의 벤젠 정량값 0.039, 0.040 wt%가
    보고되었다는 사실을 근거로, EPA 제품에 해당 농도의 벤젠이 있다고
    가정한 '반사실적' 시나리오를 만듭니다.
(c) 추가적으로 0.099, 0.100 wt% 경계값과 0.000 wt% 기준 조성을 계산합니다.
    벤젠을 추가한 만큼 EPA의 원성분 비율을 동일 비율로 줄여 합계를 100%로
    유지합니다. EPA 제품에서 벤젠이 실제 검출된 것은 아닙니다.
(d) 계산 온도는 정확히 25°C입니다. 첨부 참조표의 포화증기압이 25°C
    고정값이므로, 다른 온도로 계산하려면 포화증기압 데이터를 다시 구축해야 합니다.
(e) 혼합물 증기 몰분율 y_i = x_i * gamma_i * P_i(sat) / Σ(x_j * gamma_j * P_j(sat))
    Raoult: gamma=1; UNIFAC/HANNA: 해당 모델 계산 gamma.
(f) 상대 중요도 = (y_i / NIOSH REL_i) / Σ(y_j / NIOSH REL_j).
    이는 단순 상대 'OEL 가중 스크리닝 지표'일 뿐입니다.
    실제 근로자 노출, 암 발생위험, 검증된 혼합물 독성 가산지수로 해석하면 안 됩니다.
(g) 한국 산업안전보건법상 정확한 기재 생략 판단에는 각 CAS의 적용 유해성 분류,
    별표 6, 별도 의무 및 비공개 대체정보 등을 각각 검토해야 합니다.
    본 코드는 액상 0.1% 미만 여부만 표시하고 적법성 판정은 하지 않습니다.
(h) HANNA/UNIFAC에서 gamma가 크거나 작아지는 경우에도 모델 적용영역과
    실제 기액평형 실험 데이터로 별도의 검증을 하기 전에는 예측 결과입니다.
(i) 같은 계산조건에서의 모델 간 차이만으로 어떤 모델이 더 정확하다고 결론
    내릴 수 없습니다. 실험 VLE 원자료에 대한 독립 검증이 후속 과제입니다.

6. 참고 출처
-----------
HANNA (multicomponent): https://github.com/marco-hoffmann/HANNA
Thermo UNIFAC: https://thermo.readthedocs.io/activity_coefficients.html
NIST VP: https://webbook.nist.gov/
NIOSH NPG: https://www.cdc.gov/niosh/npg/
우성준·하권철(2022), 한국산업보건학회지 32(4): 340-349.

7. 오류 예시
------------
'No module named thermo' : 활성화한 Python 환경에서 pip install thermo
'HANNA repository absent': D:\fire\HANNA 폴더에 clone 필요
'HANNA checkpoints missing': 공식 pretrained model 폴더/파일 확인
'UNIFAC group assignment unavailable': 해당 CAS 매핑이 없어 계산 제외
'04 input not found': 실제 04 CSV 위치를 확인하거나 --input으로 경로 설정

연구용 파일이므로, 최종 논문 작성 시에는 모든 물성/노출기준의 출처,
모델 버전, 검증 오차와 법적 성분별 분류를 반드시 기록하십시오.
