"""Download CMS 2011A MOD files (data and SIM) from Zenodo, md5-checked, into data/raw/mod/.

Uses energyflow's own file table (Zenodo record per sub-dataset, md5 per file) and its fetcher,
so the files land exactly where ``analyze_mod.py``, ``analyze_mod_sim.py`` and ``unfold_mod.py``
look for them. Already present files with the right md5 are not downloaded again.

Usage:
  python src/download_mod.py CMS:3 SIM170:1 SIM300:6 SIM470:6 SIM600:5 SIM800:3 SIM1000:2 \
      SIM1400:2 SIM1800:2            # the inputs of docs/mod_sim_results.md (about 6 GB)

``NAME:N`` downloads files 0 .. N-1 of a sub-dataset (CMS = the data; SIM170 ... SIM1800 = the
pT-hat slices of the simulation). Three files at a time took 22 min for the above (2026-11-24).
"""

from __future__ import annotations

import argparse
import pathlib
import time
from concurrent.futures import ThreadPoolExecutor

REPO = pathlib.Path(__file__).resolve().parents[1]
RAW = REPO / "data" / "raw" / "mod"


def jobs(specs):
    import energyflow.datasets.mod as m

    md5 = m._get_dataset_info("CMS2011AJets")["md5_hashes"]
    record = {s[0]: s[2] for grp in m.COLLECTIONS["CMS2011AJets"].values()
              if isinstance(grp, dict) and "subdatasets" in grp for s in grp["subdatasets"]}
    out = []
    for spec in specs:
        name, n = spec.split(":")
        full = "CMS_Jet300_pT375-infGeV" if name == "CMS" else f"{name}_Jet300_pT375-infGeV"
        if full not in record:
            raise SystemExit(f"unknown sub-dataset {name}")
        for i in range(int(n)):
            fn = f"{full}_{i}_compressed.h5"
            out.append((fn, m.ZENODO_URL_PATTERN.format(record[full], fn),
                        f"datasets/CMS2011AJets/{full}", md5[fn]))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("specs", nargs="+", help="NAME:N_FILES, e.g. CMS:3 SIM300:6")
    ap.add_argument("--parallel", type=int, default=3)
    args = ap.parse_args()
    from energyflow.utils.data_utils import _get_filepath

    todo = jobs(args.specs)
    for _, _, sub, _ in todo:          # create directories up front (threads race otherwise)
        (RAW / sub).mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    def get(job):
        fn, url, sub, h = job
        for attempt in range(4):
            try:
                _get_filepath(fn, url, str(RAW), cache_subdir=sub, file_hash=h)
                print(f"ok {fn} ({time.time() - t0:.0f} s)", flush=True)
                return True
            except Exception as e:  # network errors: retry with backoff
                print(f"retry {fn}: {e}", flush=True)
                time.sleep(10 * 2**attempt)
        print(f"FAILED {fn}", flush=True)
        return False

    with ThreadPoolExecutor(args.parallel) as ex:
        ok = list(ex.map(get, todo))
    print(f"{sum(ok)}/{len(ok)} files in {time.time() - t0:.0f} s")
    return 0 if all(ok) else 1


if __name__ == "__main__":
    raise SystemExit(main())
