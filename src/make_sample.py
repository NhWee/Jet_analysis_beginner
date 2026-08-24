"""연습용 ROOT 파일을 만든다.

진짜 Open Data 를 받기 전에, 같은 형식의 작은 파일로 코드를 먼저 익히기 위한 것이다.
파일 구조는 실제 실험 데이터와 같다 — 이벤트마다 입자 목록이 들어있는 ROOT TTree.

    python src/make_sample.py

결과: data/raw/sample_events.root  (이벤트 200개, 약 0.5 MB)

주의: 여기 들어있는 것은 실제 충돌 데이터가 아니라 흉내낸 것이다.
코드를 익히는 용도이고, 물리 결론을 내는 데 쓰면 안 된다.
"""

from __future__ import annotations

import pathlib
import sys

import awkward as ak
import uproot

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from toy_generator import ToyConfig, generate  # noqa: E402

REPO = pathlib.Path(__file__).resolve().parents[1]
OUT = REPO / "data" / "raw" / "sample_events.root"


def main() -> int:
    n_events = 200
    seed = 12345

    print(f"이벤트 {n_events}개 생성 중 (seed={seed}) ...")
    particles = generate(n_events, seed=seed, cfg=ToyConfig())

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with uproot.recreate(OUT) as f:
        # 실제 데이터처럼 이벤트당 입자 목록을 담는다.
        # 이름은 CMS NanoAOD 관례를 따랐다 (Particle_* / PFCand_*).
        f["Events"] = {
            "Particle_pt": particles.pt,      # 횡운동량 [GeV]
            "Particle_eta": particles.eta,    # 유사래피디티 (빔축 방향 각도)
            "Particle_phi": particles.phi,    # 방위각 [rad]
            "Particle_mass": particles.mass,  # 질량 [GeV]
        }

    size_kb = OUT.stat().st_size / 1024
    print(f"저장: {OUT.relative_to(REPO)}  ({size_kb:.0f} kB)")
    print(f"이벤트당 평균 입자 수: {ak.mean(ak.num(particles)):.0f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
