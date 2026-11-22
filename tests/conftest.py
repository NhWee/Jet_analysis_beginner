"""Shared pytest setup: headless plotting and small toy samples. Tests run on the
dependency-free toy generator, so they need no Pythia and no data and finish in seconds."""

from __future__ import annotations

import pathlib
import sys

import matplotlib

matplotlib.use("Agg")

SRC = pathlib.Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(SRC))  # also set in pyproject.toml; kept for running without it

import awkward as ak  # noqa: E402
import pytest  # noqa: E402

from jeteec import jets as jets_mod  # noqa: E402
from toy_generator import ToyConfig, generate  # noqa: E402


@pytest.fixture(scope="session")
def toy_particles():
    """300 toy events with a ``charge`` field — small, but every feature is present."""
    return generate(300, seed=12345, cfg=ToyConfig())


@pytest.fixture(scope="session")
def toy_jets(toy_particles):
    """(jets, constituents) flattened over events: one entry per jet."""
    jets, const = jets_mod.cluster(toy_particles)
    return ak.flatten(jets), ak.flatten(const, axis=1)
