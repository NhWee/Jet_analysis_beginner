"""ROOT 파일 안에 무엇이 들어있는지 보여준다.

**새 데이터를 받으면 분석 전에 항상 이것부터 실행한다.**
실험마다, 릴리스마다 변수 이름이 다르기 때문에 이름을 확인하지 않고 짠 코드는 반드시 깨진다.

사용법
------
로컬 파일:
    python src/inspect_file.py data/raw/sample_events.root

인터넷 URL (uproot 는 URL 도 직접 연다 — 다운로드 불필요):
    python src/inspect_file.py https://.../파일.root

특정 단어가 들어간 변수만 보기:
    python src/inspect_file.py 파일.root --grep jet
"""

from __future__ import annotations

import argparse

import uproot


def show_tree(tree, grep: str | None, max_branches: int) -> None:
    print(f"  이벤트 수: {tree.num_entries:,}")

    names = list(tree.keys())
    if grep:
        low = grep.lower()
        names = [n for n in names if low in n.lower()]
        print(f"  '{grep}' 를 포함하는 변수: {len(names)}개")
    else:
        print(f"  변수(branch) 수: {len(names)}개")

    if not names:
        print("  (해당하는 변수 없음)")
        return

    print()
    width = min(max(len(n) for n in names), 45)
    for n in names[:max_branches]:
        b = tree[n]
        # typename 끝의 [] 는 '이벤트마다 길이가 다른 목록' 이라는 뜻이다.
        # 그런 변수가 입자/jet 처럼 개수가 변하는 대상을 담는다.
        print(f"    {n:<{width}}  {b.typename}")
    if len(names) > max_branches:
        print(f"    ... 그 외 {len(names) - max_branches}개 "
              f"(--max 로 더 보기)")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("path", help="ROOT 파일 경로 또는 URL")
    ap.add_argument("--grep", help="이 단어가 들어간 변수만 표시 (예: jet)")
    ap.add_argument("--max", type=int, default=60, help="표시할 변수 최대 개수")
    args = ap.parse_args()

    print(f"열기: {args.path}")
    try:
        f = uproot.open(args.path)
    except Exception as e:                                   # noqa: BLE001
        print(f"\n실패: {type(e).__name__}")
        print(f"  {e}")
        print("\n확인할 것:")
        print("  - 경로/URL 이 맞는가")
        print("  - URL 이면 인터넷이 되는가 (회사/학교 방화벽 확인)")
        return 1

    print(f"\n파일 안의 객체: {f.keys()}")

    trees = {k: v for k, v in f.items() if hasattr(v, "num_entries")}
    if not trees:
        print("\nTTree 가 없다. 이 파일은 히스토그램만 담고 있을 수 있다.")
        return 1

    for name, tree in trees.items():
        print(f"\n=== TTree: {name} ===")
        show_tree(tree, args.grep, args.max)

    print("\n" + "-" * 60)
    print("다음 단계: 위 변수 이름을 노트북의 이름과 맞춘다.")
    print("  jet 이 이미 재구성돼 있으면  -> 클러스터링 단계를 건너뛴다")
    print("  입자(constituent)가 있으면   -> 노트북 그대로 쓸 수 있다")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
