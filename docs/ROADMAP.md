# 보완 사항과 다음 단계 (Roadmap)

2026-11 기준 저장소 상태를 점검한 결과다. **현재 위치: 파이프라인은 toy 이벤트로만 검증됐고,
Pythia 검증 결과와 실데이터 결과는 저장소에 아직 없다.**

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

- **경로 (a) CMS 2011A MOD 데이터셋 — 먼저 확인할 가치가 큼.**
  Komiske–Moult–Thaler–Zhu ([arXiv:2201.07800](https://arxiv.org/abs/2201.07800)) 가 쓴
  CMS 2011 Jet primary dataset 을 PFCandidate 포함 Python 형식으로 변환한 공개본이다.
  `energyflow` 패키지에서 읽을 수 있는 것으로 알려져 있다. 사실이라면 **CMSSW 없이 바로 EEC 를
  실데이터로 계산**할 수 있고, 같은 논문 결과와 직접 비교할 수 있다.
  ⚠️ 접근 방법·현재 가용성은 확인 필요. 확인되면 [tier_physics_programme.md](tier_physics_programme.md)
  §0 의 "EEC 는 CMSSW 가 최소 한 번 필요" 문장을 수정해야 한다.
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

1. `tests/test_eec.py`, `tests/test_jets.py` 작성 (A1)
2. `validate_eec_pythia.py` 실행하고 결과 기록 (A4)
3. CMS 2011A MOD 데이터셋 접근 가능 여부 확인 (단계 3a) — 가능하면 전체 로드맵이 크게 짧아진다
