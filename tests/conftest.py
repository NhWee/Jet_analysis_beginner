"""Shared pytest setup: headless plotting and small toy samples."""

import matplotlib

matplotlib.use("Agg")

import pytest  # noqa: E402

from toy_generator import ToyConfig, generate  # noqa: E402


@pytest.fixture(scope="session")
def toy_particles():
    """300 toy events with a ``charge`` field — small, but every feature is present."""
    return generate(300, seed=12345, cfg=ToyConfig())
