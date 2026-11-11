"""The beginner notebook must keep running against the current code and sample format.

Executes the code cells in order, without Jupyter, on a freshly made sample file.
"""

import json
import pathlib
import warnings

import matplotlib.pyplot as plt

import make_sample

NOTEBOOK = (pathlib.Path(__file__).resolve().parents[1]
            / "notebooks" / "01_what_is_a_jet.ipynb")


def test_notebook_runs_top_to_bottom(tmp_path):
    sample = tmp_path / "sample_events.root"
    assert make_sample.main(["--out", str(sample)]) == 0

    cells = [c for c in json.loads(NOTEBOOK.read_text(encoding="utf-8"))["cells"]
             if c["cell_type"] == "code"]
    ns = {"__name__": "__notebook__"}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")          # plt.show() on a headless backend
        for k, cell in enumerate(cells):
            code = compile("".join(cell["source"]), f"{NOTEBOOK.name}[code cell {k}]", "exec")
            exec(code, ns)
            if "DATA" in ns:
                ns["DATA"] = sample              # point the notebook at the fresh sample
    plt.close("all")

    # the closure check the notebook teaches must actually hold
    assert abs(ns["total"].pt - ns["jet"].pt) < 1e-9
    # and the jets it labels "jet 0, jet 1, ..." must be in pT order
    jets = ns["jets"]
    assert all(list(j) == sorted(j, reverse=True) for j in jets.pt.to_list())
