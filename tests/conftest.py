"""Shared fixtures. Tests run on the dependency-free toy generator, so they need no
Pythia and no data and finish in seconds."""

from __future__ import annotations

import pathlib
import sys

import awkward as ak
import pytest

SRC = pathlib.Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(SRC))

from jeteec import jets as jets_mod  # noqa: E402
from toy_generator import generate  # noqa: E402


@pytest.fixture(scope="session")
def toy_particles():
    """Toy events with an extra ``charge`` field, as real PF candidates would carry."""
    return generate(200, seed=2024)


@pytest.fixture(scope="session")
def toy_jets(toy_particles):
    """(jets, constituents) flattened over events: one entry per jet."""
    jets, const = jets_mod.cluster(toy_particles)
    return ak.flatten(jets), ak.flatten(const, axis=1)
