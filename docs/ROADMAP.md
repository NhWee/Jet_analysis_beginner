# 보완 사항과 다음 단계 (Roadmap)

2026-11 기준 저장소 상태를 점검한 결과다.
**현재 위치 (2026-11-05 갱신): 파이프라인은 toy 와 Pythia truth 로 검증됐다. 실데이터 결과는 아직 없다.**

## 진행 상황

| # | 상태 | 비고 |
|---|---|---|
| A1 테스트 | ✅ | `tests/` 22개, 1초 내. `python -m pytest -q tests` |
| A2 CI | ✅ | `.github/workflows/tests.yml`: pytest + toy tier 분석 + Pythia 검증(소규모) |
| A3 lock 파일 | ⚠️ 부분 | pip 고정 목록 `requirements.txt` 만 있음. **conda `environment.lock.yml` 은 실제 사용하는 conda 환경에서 export 해야 한다** (이 작업 환경에는 conda 가 없었음) |
| A4 Pythia 검증 | ✅ | PASS, 결과는 [pythia_validation_results.md](pythia_validation_results.md). sample register P8-001/002 |
| A5 통계 오차 | ✅ | `eec_with_errors()`: jet 단위 SEM. 최종 결과에는 event 단위 bootstrap 권장 |
| A6 pT 구간 | ✅ | `eec_in_pt_bins()`. Pythia 에서 작은 각도 쪽이 R_L·pT 로 겹치는 것 확인 |
| A7 정규화 규약 | ✅ | `norm="sum_pt"`(기본) / `"jet_pt"`(CMS 방식) |
| A8 메모리 | ✅ 부분 | `chunk_size` 로 블록 처리, 메모리 상한 확보. numba 루프는 아직 필요 없어 보류 |
| B1 hadronisation on/off | ✅ | `--compare-hadronisation`. R_L < 0.01 에서 parton/hadron 모양 비 0.28 |

새로 발견한 함정: **EEC 를 서로 다른 constituent 수의 샘플끼리 비교할 때 self-pair 가중치(Σz²)
차이 때문에 보이는 bin 전체가 위아래로 밀린다.** 모양 비교는 공통 범위로 재정규화해서 한다.
자세한 건 [pythia_validation_results.md](pythia_validation_results.md) §3.

C3 (MOD 데이터셋): **존재·공개 확인, 다운로드는 아직** — 단계 3 경로 (a) 참고.

MOD 조건 Pythia 기준선 (7 TeV, R = 0.5, leading 2, pT > 375 GeV): ✅ PASS. pT 구간 4개 모두 작은 각도 쪽이
R_L·pT 로 한 곡선에 겹침, turnover ≈ 7–9 GeV. [pythia_validation_results.md](pythia_validation_results.md) §4.
pTHatMin 을 임계값 근처에 두면 jet 의 3–9% 를 잃는다는 것을 확인하고 pT-hat bias 가중 생성으로 해결.

다음 할 일: 네트워크에서 `zenodo.org` 허용 후 MOD 파일 1개 받아 구조 확인 → 첫 실데이터 EEC 를 이 기준선과 비교.

판단 근거와 진행 기록은 Notion 연구 허브에 남긴다 ([WORKFLOW.md](WORKFLOW.md) §4).
이 문서는 "무엇을, 어떤 순서로, 무엇이 되면 끝인가" 만 적는다.

---

## A. 먼저 보완할 부분 (실데이터 전에)

실데이터로 넘어가면 버그의 원인이 코드인지 데이터인지 구분하기 어려워진다.
아래는 지금 고치는 편이 훨씬 싸다.

| # | 문제 | 왜 중요한가 | 완료 기준 |
|---|---|---|---|
| A1 | **자동 테스트가 없다** (`tests/` 없음) | EEC 정규화 = 1, fastjet 필드 함정 가드 같은 핵심 불변량이 스크립트 출력으로만 확인된다. 깨져도 아무도 모른다 | `pytest` 로 정규화, `four()` 가드, ΔR wrap, 빈 입력 에러를 검사. toy 생성기 기반이라 수 초 안에 끝남 |
| A2 | **CI 가 없다** | A1 을 사람이 기억해서 돌려야 한다 | GitHub Actions 에서 push 마다 `pytest` + `tier_analyses.py` 소규모 실행 |
| A3 | **`environment.lock.yml` 이 없다** | 재현성 규칙(README)에 커밋하라고 적혀 있지만 실제로는 없다. fastjet/awkward 버전 차이는 결과를 조용히 바꿀 수 있다 | 실제 사용 환경에서 export 후 커밋 |
| A4 | **Pythia 검증 결과가 기록돼 있지 않다** | README 의 "Quick start 1번은 반드시 통과" 가 아직 증명되지 않았다. 모든 수치가 toy 출신 | `validate_eec_pythia.py` 실행 결과(정규화, ⟨n_const⟩, EEC 그림)를 seed 와 함께 `docs/` 에 기록, `data/README.md` sample register 에 등록 |
| A5 | **`eec()` 가 통계 오차를 반환하지 않는다** | 오차 없는 EEC 는 데이터·이론과 비교할 수 없다. 같은 jet 안의 쌍은 서로 상관되므로 단순 √N 은 틀린다 | jet 단위 bootstrap 또는 jet별 히스토그램의 분산으로 오차 계산 |
| A6 | **EEC 가 jet pT 구간별로 나뉘지 않는다** | 이 프로젝트의 핵심 물리(turnover 위치가 Λ_QCD/pT 로 움직임)는 pT 구간별 EEC 가 있어야 보인다. 지금은 모든 jet 평균 하나 | `eec()` 에 pT bin 인자 추가, `R_L·pT` 로 그렸을 때 turnover 가 겹치는지 확인 |
| A7 | **정규화 규약이 문헌과 다르다** | 여기서는 `Σ pT_const` 로 나누고 CMS ([arXiv:2402.13864](https://arxiv.org/abs/2402.13864)) 는 jet pT 기준이다. 숫자를 직접 비교하면 어긋난다 | 두 규약을 옵션으로 제공하고, 비교 시 어느 쪽을 썼는지 명시 |
| A8 | **쌍 계산이 O(N²) 메모리** | `ak.argcombinations` 는 모든 쌍을 메모리에 만든다. 고 pT jet(수십~100 개 constituent) × 수백만 jet 이면 터진다. `numba` 는 환경에 있지만 쓰이지 않는다 | 파일 단위 chunk 처리 + 필요 시 numba 루프로 히스토그램에 직접 누적 |

권장 순서: **A1 → A4 → A5/A6 → A2/A3 → A7/A8.**
A1 이 있어야 A5~A8 을 고칠 때 기존 결과가 깨지지 않았음을 보장할 수 있다.

---

## B. 다음 단계 (단계별)

각 단계는 이전 단계가 known-good 일 때만 시작한다.

### 단계 1 — Pythia 로 물리 기준선 만들기
- A4 완료. Pythia truth 에서 pT 구간별 EEC 를 그리고 세 영역(free hadron / turnover / scaling)이 보이는지 확인
- 강입자화 on/off (`HadronLevel:all = off`) 비교 → turnover 가 강입자화에서 온다는 것을 직접 확인
- 쿼크 jet vs 글루온 jet EEC 비교 (Pythia 파톤 정보로 태깅)
- **완료 기준:** turnover 위치가 대략 1/pT 로 움직이는 그림 1장

### 단계 2 — 첫 실데이터: NanoAOD (constituent 불필요)
- CMS 2016 NanoAOD 한 run 으로 inclusive jet 스펙트럼, R₃₂
- 여기서 배우는 것: golden JSON(lumi mask), trigger 선택과 turn-on, JEC 가 이미 적용됐는지 확인
- 첫 샘플을 `data/README.md` sample register 에 record ID 와 함께 등록
- **완료 기준:** toy 에서 돌린 NanoAOD 분석과 같은 코드가 실데이터에서 돌아간다

### 단계 3 — 첫 실데이터 EEC
두 경로 중 하나를 고른다 (결정 기록을 Notion 에 남길 것).

- **경로 (a) CMS 2011A MOD 데이터셋 — 존재와 공개 여부 확인됨 (2026-11-05).**
  CMS 2011A Jet primary dataset 의 AOD 를 MIT 그룹이 PFCandidate 포함 HDF5 로 재처리한 공개본.
  `energyflow` 로 CMSSW 없이 읽는다.

  | 항목 | 내용 | 확인 경로 |
  |---|---|---|
  | 데이터 | Zenodo [10.5281/zenodo.3340205](https://doi.org/10.5281/zenodo.3340205), CC BY 4.0, 2019-08 공개 | 웹 검색 |
  | 내용 | 1,785,625 jets, Jet300 trigger, 이벤트의 leading 2 jet, **anti-kT R = 0.5**, 보정 후 **pT > 375 GeV**, PFCandidate 포함 | 웹 검색 |
  | MC | Pythia 6 + CMS 검출기 시뮬레이션(SIM)과 generator-level(GEN), 생성 pT 구간 8개 (170 GeV ~ ∞) | `energyflow` 소스 |
  | 로더 | `energyflow` 1.4.0 `energyflow.datasets.mod.load(collection="CMS2011AJets", dataset="cms")` | `energyflow` 소스 |
  | 파일 | 로더 기준 18개, 파일당 약 111 MB → 전체 약 2 GB | 로더 소스 + 웹 검색 (파일 수 미검증) |

  ⚠️ **이 작업 환경에서는 다운로드하지 못했다** (`zenodo.org` 가 네트워크 정책에 막힘). 파일을 열어 본
  것은 아니므로 실제 branch 구조와 파일 수는 처음 내려받을 때 확인할 것.

  분석상 주의:
  - **R = 0.5, pT > 375 GeV 로 고정**이다. 지금 Pythia 기준선(R = 0.4, pT 100–500)과 겹치지 않는다.
    데이터와 비교할 Pythia 기준선을 같은 조건으로 다시 만들어야 한다.
  - Run 1 (7 TeV) 데이터다. CMS αs 측정(Run 2, 13 TeV)과 직접 비교 대상이 아니다.
  - SIM/GEN 쌍이 함께 있으므로 **검출기 보정(unfolding) 연습까지 이 데이터셋 하나로** 가능하다.
- **경로 (b) MiniAOD → PFNano (CMSSW Docker 1회).** Run 2, 높은 통계, CMS αs 측정과 같은 조건.
  작업량은 크지만 최종 목적지.

**완료 기준:** 실데이터 EEC 정규화 = 1, 그리고 공개된 결과(경로 a: 2201.07800, 경로 b: HEPData)와 형태 비교

### 단계 4 — 측정다운 측정
- 검출기 효과 보정: Pythia + 검출기 시뮬레이션(또는 공개 MC)으로 unfolding
- 체계 오차: JES/JER, 트래킹 효율(charged-only EEC), pileup
- E3C/E2C 비율로 αs 추출 시도 (CMS 결과 재현이 목표)

### 단계 5 — 딥러닝
- 먼저 NanoAOD jet-level feature 로 quark/gluon 분류 (워크플로우 연습)
- 그다음 constituent 기반 EFN/PFN (`energyflow`), 입력에 EEC 를 넣었을 때와 비교

---

## C. 바로 다음 한 주에 할 일

1. ~~`tests/test_eec.py`, `tests/test_jets.py` 작성 (A1)~~ ✅
2. ~~`validate_eec_pythia.py` 실행하고 결과 기록 (A4)~~ ✅
3. ~~CMS 2011A MOD 데이터셋 접근 가능 여부 확인 (단계 3a)~~ ✅ 공개 확인, 다운로드는 네트워크 허용 후
