"""Conversion of MOD-format PF candidates into the pipeline's constituent arrays.

Uses tiny hand-written arrays in the documented MOD column order
(pt, y, phi, m, pid, vertex); no data file is needed."""

from __future__ import annotations

import awkward as ak
import numpy as np
import pytest

pytest.importorskip("energyflow")

from analyze_mod import pfcs_to_constituents  # noqa: E402
from jeteec import eec  # noqa: E402

JET = np.array([
    # pt     y     phi   m      pid  vertex
    [100.0, 0.10, 0.20, 0.000,   22, 0],   # photon
    [50.0, 0.12, 0.25, 0.1396, 211, 0],   # pi+ from leading vertex
    [5.0, 0.30, 0.10, 0.1396, -211, 2],   # pi- from a pileup vertex
    [3.0, 0.05, 0.30, 0.4977, 130, 0],    # K0L (neutral)
])


def test_chs_drops_only_charged_pileup():
    c = pfcs_to_constituents([JET], chs=True)
    np.testing.assert_allclose(ak.to_numpy(c.pt[0]), [100.0, 50.0, 3.0])
    c = pfcs_to_constituents([JET], chs=False)
    assert len(c[0]) == 4


def test_charged_only():
    c = pfcs_to_constituents([JET], chs=False, charged_only=True)
    np.testing.assert_allclose(ak.to_numpy(c.pt[0]), [50.0, 5.0])


def test_drop_pids_uses_abs_pdg_id():
    c = pfcs_to_constituents([JET], chs=False, drop_pids=(130,))
    np.testing.assert_allclose(ak.to_numpy(c.pt[0]), [100.0, 50.0, 5.0])
    c = pfcs_to_constituents([JET], chs=False, drop_pids=(211,))
    np.testing.assert_allclose(ak.to_numpy(c.pt[0]), [100.0, 3.0])


def test_rapidity_converted_to_eta():
    c = pfcs_to_constituents([JET], chs=False)
    eta = ak.to_numpy(c.eta[0])
    assert eta[0] == pytest.approx(0.10)          # massless: eta == y
    assert abs(eta[3]) > abs(JET[3, 1])           # massive: |eta| > |y|
    np.testing.assert_allclose(ak.to_numpy(c.rapidity[0]), JET[:, 1], atol=1e-9)


def test_jagged_layout_and_normalisation():
    c = pfcs_to_constituents([JET, JET[:2]], chs=False)
    assert ak.to_list(ak.num(c)) == [4, 2]
    assert eec.normalisation_check(c) == pytest.approx(1.0, abs=1e-12)
