# 작업 방식 (Workflow)

이 프로젝트는 **결과물뿐 아니라 과정을 남기는 것**이 목표의 일부다.
분석을 돌리고 끝내는 것이 아니라, 나중의 나와 다음 사람이 같은 경로를 다시 걸을 수
있어야 한다.

기록은 두 곳에 나눠 남긴다.

| 어디에 | 무엇을 | 왜 |
|---|---|---|
| **GitHub** (이 저장소) | 코드, 설정, 재현 방법, **결과 수치의 정본** | 실행 가능한 상태를 정확히 보존 |
| **Notion** | 연구 노트, **결과의 해석**, 판단 근거, 실패 기록 | 왜 그렇게 했는지는 코드에 안 남는다 |

원칙: **"무엇을 했다"는 GitHub, "왜 그렇게 했고 무엇이 안 됐다"는 Notion.**

### 연구 결과는 어디에 쓰나 — 양쪽에, 형태를 달리해서

헷갈리기 쉬운 지점이라 명시한다. **결과는 양쪽 모두에 들어간다.** 다만 담는 형태가 다르다.

| | GitHub | Notion |
|---|---|---|
| 성격 | **정본 (source of truth)** | 해석과 맥락 |
| 담는 것 | 전체 표, 모든 수치, 그림, 실행 명령, seed | 대표 숫자 1\~2개 + 그것이 뜻하는 바 |
| 갱신되면 | 코드를 다시 돌려 덮어쓴다 | 그날의 기록으로 남긴다 (덮어쓰지 않는다) |
| 예 | `docs/tier_analyses_results.md` 의 전체 표 | "n = 5.04 ± 0.07, 생성값 5.0 과 일치. 전에 적은 n = 4.85 와 '3% 는 threshold 효과' 라는 해석은 틀렸다 — 비가중 log 공간 fit 의 편향이었다. 추정량도 검증 대상이다" |

**판단 기준: 코드를 다시 돌려서 나오면 GitHub, 사람의 판단이 들어가면 Notion.**

수치 자체는 재현 가능하므로 저장소가 정본이다. 반면 "이 숫자가 왜 이런가", "예상과 달랐는가",
"그래서 다음에 무엇을 할 것인가" 는 재현되지 않으므로 Notion 에만 남는다.

**중복은 최소화한다.** 같은 표를 양쪽에 복사하면 한쪽만 갱신됐을 때 어긋난다. Notion 에는
대표 숫자와 **커밋 해시**만 적고 전체는 저장소를 보게 한다. 반대로 GitHub 문서에는 그날의
추측이나 미확정 판단을 적지 않는다 — 코드와 함께 갱신될 성질이 아니기 때문이다.

**Notion 기록은 사후 수정하지 않는다.** 나중에 틀린 것으로 밝혀져도 지우지 말고 새 항목에서
정정한다. 틀렸던 기록이 남아야 왜 그렇게 판단했는지 추적할 수 있다.
(실제 예: PFNano 를 공개 다운로드 형식으로 오판했던 건이 실패 보고 DB 에 그대로 남아 있다.)

---

## 1. 저장소 구조

```
Jet_analysis_beginner/
├── README.md                 셋업과 실행 순서. 신규 진입자의 첫 페이지
├── environment.yml           conda 환경 정의 (전체: Pythia, PyTorch 포함)
├── requirements.txt          pip 최소 환경 (toy 파이프라인, 노트북, 테스트)
├── requirements-lock.txt     문서 수치를 만든 정확한 버전
├── Dockerfile                재현용 컨테이너
├── docs/
│   ├── WORKFLOW.md           ← 지금 이 문서
│   ├── data_tiers.md         티어별로 무엇이 가능한지 (메커니즘)
│   ├── tier_physics_programme.md  티어별로 무엇을 측정할 수 있는지 (물리)
│   └── tier_analyses_results.md   실제 돌린 결과와 수치, 찾은 버그
├── src/
│   ├── jeteec/               재사용 라이브러리 (jets, eec, stats, kinematics, plot)
│   ├── toy_generator.py      의존성 없는 대체 샘플 생성기 (truth 기록 포함)
│   ├── validate_eec_pythia.py  실데이터 전 필수 검증
│   ├── make_sample.py        노트북용 연습 ROOT 파일
│   ├── inspect_file.py       ROOT 파일 내용 확인 — 새 파일은 항상 이것부터
│   ├── tier_demo.py          티어별 가능/불가능 비교
│   └── tier_analyses.py      티어별 실제 분석 4종
├── tests/                    단위·closure·문서 수치·노트북 테스트 (pytest)
├── data/                     내용은 커밋하지 않음. sample register 만 기록
├── results/figures|tables/   그림·표. 커밋하지 않음 — 단 문서가 참조하는 tier_*.png 는 예외
├── models/                   체크포인트. 커밋하지 않음
└── notebooks/                튜토리얼·탐색용. 커밋 전 output 제거
```

---

## 2. 재현 방법

```bash
# 환경 (또는 pip install -r requirements.txt)
mamba env create -f environment.yml
conda activate jet-eec

# 0) 테스트: 단위 검사, closure, 문서에 적힌 수치, 노트북 실행
pytest

# 1) 파이프라인이 옳은지 먼저 확인 (실데이터 전 필수)
python src/validate_eec_pythia.py --n-events 2000 --seed 12345

# 2) 어떤 데이터 티어가 필요한지 확인
python src/tier_demo.py --n-events 400 --seed 12345

# 3) 티어별 분석 실행
python src/tier_analyses.py --n-events 2000 --seed 12345
```

모든 스크립트는 `--seed` 를 받는다. **결과를 기록할 때는 seed 를 반드시 함께 적는다.**
seed 없는 숫자는 재현할 수 없으므로 기록으로서 가치가 없다.

**숫자에는 통계 오차를 붙인다.** 오차 없는 "3% 차이" 는 해석할 수 없다 — 이 저장소의 AOD
결과가 그랬다 (오차를 붙이니 모두 1σ 안이었다). toy 처럼 답을 아는 경우에는 그 답도 옆에 적는다.

`tests/test_tier_numbers.py` 가 `docs/tier_analyses_results.md` 의 수치를 고정한다. 코드를 고쳐
수치가 바뀌는 것이 의도된 것이면, 그림·문서·테스트를 **같은 커밋에서** 함께 갱신한다.

---

## 3. 커밋 규칙

한 커밋 = 한 가지 의미 단위. 제목은 영어 prefix + 한국어 요약, 본문은 **왜**를 적는다.

```
<type>: 한국어 요약

무엇을 왜 그렇게 했는지. 대안이 있었다면 왜 그것을 택하지 않았는지.
숫자가 있으면 seed 와 함께 적는다.
```

| type | 쓰는 경우 |
|---|---|
| `feat` | 새 기능·분석 추가 |
| `fix` | 버그 수정 (**무엇이 틀렸는지 본문에 반드시 기록**) |
| `docs` | 문서 |
| `test` | 검증·테스트 |
| `build` | 환경·의존성 |
| `chore` | 잡무 |

커밋 본문에 남긴 것은 `git log` 로 영원히 검색된다. 특히 `fix` 는 실패 기록의
일부이므로, 증상·원인·재발 방지책을 세 줄이라도 적는다.

**절대 커밋하지 않는 것:** 데이터 파일(`.root`, `.parquet`), 모델 체크포인트,
그림·표, notebook output. 모두 재생성 가능해야 하고, 재생성이 안 되면 그건
파이프라인이 불완전하다는 뜻이다.

예외 하나: README 와 `docs/` 가 보여주는 `results/figures/tier_*.png` 는 커밋한다 (없으면
GitHub 에서 그림이 깨진다). 고정 seed 로 `requirements-lock.txt` 환경에서 재생성한 것만 올린다.

---

## 4. Notion 기록 규칙

Notion에는 세 종류를 남긴다.

**작업 로그 (Work Log DB)** — 하루 또는 한 세션 단위로 한 항목.

- 무엇을 하려 했는가 (목표)
- 무엇을 했는가 (실행)
- 무엇을 알았는가 (결과, 숫자는 seed 와 함께)
- 다음에 무엇을 할 것인가
- 관련 커밋 해시

**실패 보고 (Failure Report DB)** — goal 에 명시된 항목. 잘 안 된 것마다 하나.

- 증상: 무엇이 어떻게 잘못됐는가
- 원인: 왜 그랬는가 (모르면 "미상"이라고 적는다)
- 어떻게 찾았는가 ← **가장 재사용 가치가 높은 항목**
- 재발 방지: 코드 가드, 테스트, 문서 중 무엇을 남겼는가

**결정 기록** — 갈림길에서 한쪽을 택했을 때. 나중에 "왜 이렇게 했지?" 를 막는다.

Notion 위치:

- [Jet Analysis — 연구 허브](https://app.notion.com/p/3c06974402a681278164ec1d2488baf6) (목표, 로드맵, 판단 근거)
  - 작업 로그 (Work Log) DB
  - 실패 보고 (Failure Report) DB

> 실패 보고를 남기는 기준: *같은 함정에 한 번 더 빠질 수 있다면 남긴다.*
> 사소해 보여도 30분 이상 태운 것은 전부 해당한다.

---

## 5. 지금까지의 판단 근거 (요약)

새로 합류한 사람이 "왜 이런 구조인가"를 묻지 않아도 되도록 정리한다.

**왜 conda 와 Docker 를 둘 다 쓰나.** 배타적 선택이 아니다. Docker 는 데이터 접근과
재현성(계층 1), conda 는 분석과 ML(계층 2). 표준 패턴은 Dockerfile 안에서 conda 환경을
구성하는 것이고, 이 저장소가 그렇게 돼 있다.

**왜 LHC 4개 실험만인가.** RHIC(STAR/sPHENIX)과 Fermilab(CDF/D0)은 공개 다운로드 포털이
없다. 데이터 보존은 돼 있으나 협업 멤버 접근이다. 실제로 손댈 수 있는 것은 CERN Open Data
의 CMS/ATLAS/ALICE/LHCb 뿐이다.

**왜 NanoAOD 부터 시작하나.** EEC 에는 못 쓰지만, uproot 만으로 실데이터를 다루며 selection,
trigger, JEC, unfolding 을 익힐 수 있다. 각 단계가 known-good 이어야 다음 단계의 이상을
싸게 진단할 수 있다.

**왜 PFNano 가 목적지인가.** constituent 가 있고 uproot 로 읽히며 AOD 보다 20배 작다.
단, PFNano 는 **다운로드하는 형식이 아니라** MiniAOD 위에 CMSSW 로 직접 생성하는 것이다.
CMSSW 를 최소 한 번은 거쳐야 한다.

**왜 toy generator 를 만들었나.** Pythia 컴파일 없이도 분석 코드를 돌려보기 위해서다.
물리가 아니므로 결론을 내는 데는 쓰지 않는다.

---

## 6. 알려진 함정

**fastjet 3.5.1.4 이하는 메모리에 빈틈 없이 이어지지 않은 배열을 잘못 읽는다.** fastjet 이
돌려준 constituent 를 골라내거나(acceptance cut) 순서를 바꾼 뒤 다시 `ClusterSequence` 에
넘기면 **조용히 틀린 jet** 이 나온다. 예외도 경고도 없다. 측정값: subjet 4-벡터 합이 부모 jet 과
pT 725 GeV, mass 1480 GeV 까지 어긋났다 (fastjet 3.5.1.4, 2000 이벤트, seed 12345. 처음 봤을 때는
168 GeV / 502 GeV). 첫 클러스터링은 멀쩡해 보였고 재클러스터링에서야 드러났다.

처음에는 `charge` 같은 여분 필드 탓으로 기록했지만 **틀린 진단이었다.** 필드를 지운 `four()`
가 통했던 것은 다시 묶으면서 배열이 연속으로 복사됐기 때문이다. `charge` 를 둔 채
`ak.to_packed` 만 해도 정확하다. fastjet 3.5.1.5 에서 고쳐졌다. 교훈: 고쳐졌다고 원인을 안
것은 아니다. **원인을 확인하는 대조 실험** (필드는 그대로, 배열 형태만 바꾸기) 을 했어야 했다.

가드는 `jeteec.jets.clusterable()` (`ak.to_packed`) 에 있고, `jeteec` 의 모든 클러스터링이
이것을 거친다. charge/PID 는 fastjet 이 constituent 에 그대로 달아서 돌려주므로 dR 로 다시
찾을 필요가 없다 (dR 매칭은 jet 의 3% 에서 틀렸다).

**fastjet 은 jet 을 pT 순서로 정렬해 주지 않는다.** `jets[:, 0]` 을 leading jet 으로 쓰면
3-jet 이벤트에서 엉뚱한 쌍을 고른다. `jeteec.jets.cluster` 가 정렬한다.

**fastjet 의 거리는 η 가 아니라 래피디티 y 로 잰다.** 질량에 비해 pT 가 크지 않은 부드러운
입자는 y 와 η 가 크게 달라서, η 로 잰 ΔR 과 클러스터링 결과가 어긋난다.

**fit 방법도 검증 대상이다.** log(count) 에 직선을 긋는 비가중 최소제곱 fit 은 입력 5.0 에서
4.85 를 냈고, 그 3% 를 물리 효과로 잘못 해석했다. 답을 아는 toy 에서 추정량의 편향부터
확인한다 (`tests/test_stats.py`).

**`.gitignore` 에서 `data/**` 만 쓰면 안 된다.** git 은 제외된 디렉터리 내부를 보지 않으므로
`!data/**/.gitkeep` 부정 패턴이 동작하지 않는다. `!data/**/` 를 함께 넣어야 한다.

**EEC 는 constituent 가 필요하다.** 샘플을 받기 전에 반드시 확인한다.

```python
import uproot
tree = uproot.open("your_file.root")["Events"]
print([k for k in tree.keys() if "PFCand" in k])   # 비어 있으면 그 샘플은 탈락
```
