"""연습용 ROOT 파일을 만든다.

진짜 Open Data 를 받기 전에, 같은 형식의 작은 파일로 코드를 먼저 익히기 위한 것이다.
파일 구조는 실제 실험 데이터와 같다 — 이벤트마다 입자 목록이 들어있는 ROOT TTree.

    python src/make_sample.py                         # 기본: 이벤트 200개, seed 12345
    python src/make_sample.py --n-events 1000 --seed 7

결과: data/raw/sample_events.root  (기본 설정에서 약 0.35 MB)

주의: 여기 들어있는 것은 실제 충돌 데이터가 아니라 흉내낸 것이다.
코드를 익히는 용도이고, 물리 결론을 내는 데 쓰면 안 된다.
"""

from __future__ import annotations

import argparse
import pathlib
import sys

import awkward as ak
import uproot

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from toy_generator import ToyConfig, generate  # noqa: E402

REPO = pathlib.Path(__file__).resolve().parents[1]
OUT = REPO / "data" / "raw" / "sample_events.root"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n-events", type=int, default=200)
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--out", type=pathlib.Path, default=OUT, help="출력 파일 경로")
    args = ap.parse_args(argv)

    print(f"이벤트 {args.n_events}개 생성 중 (seed={args.seed}) ...")
    particles = generate(args.n_events, seed=args.seed, cfg=ToyConfig())

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with uproot.recreate(args.out) as f:
        # 실제 데이터처럼 이벤트당 입자 목록을 담는다.
        # 이름은 <collection>_<variable> 이라는 NanoAOD 식 규칙을 따랐다
        # (실제 CMS 파일에서는 GenPart_pt, PFCands_pt 같은 이름이 된다).
        f["Events"] = {
            "Particle_pt": particles.pt,      # 횡운동량 [GeV]
            "Particle_eta": particles.eta,    # 유사래피디티 (빔축 방향 각도)
            "Particle_phi": particles.phi,    # 방위각 [rad], -π ~ π
            "Particle_mass": particles.mass,  # 질량 [GeV]
        }

    size_kb = args.out.stat().st_size / 1024
    try:
        shown = args.out.resolve().relative_to(REPO)
    except ValueError:
        shown = args.out
    print(f"저장: {shown}  ({size_kb:.0f} kB)")
    print(f"이벤트당 평균 입자 수: {ak.mean(ak.num(particles)):.0f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
