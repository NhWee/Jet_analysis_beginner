"""Jet clustering.

Thin wrapper over scikit-hep ``fastjet`` (awkward/array-oriented interface) so the rest of
the analysis never touches fastjet internals directly. Works identically on Pythia truth
particles and on experimental constituents (PFCandidates / tracks / topoclusters).
"""

from __future__ import annotations

from dataclasses import dataclass

import awkward as ak
import fastjet
import vector

vector.register_awkward()

# Algorithm name -> fastjet enum
_ALGOS = {
    "antikt": fastjet.antikt_algorithm,
    "kt": fastjet.kt_algorithm,
    "ca": fastjet.cambridge_algorithm,
}


@dataclass(frozen=True)
class JetConfig:
    """Jet definition and selection cuts.

    Parameters
    ----------
    algorithm : {"antikt", "kt", "ca"}
        Sequential recombination algorithm. Standard LHC jets use anti-kt.
    R : float
        Jet radius parameter.
    min_pt : float
        Minimum jet pT in GeV.
    max_abs_eta : float
        Maximum |eta| for accepted jets (detector acceptance).
    """

    algorithm: str = "antikt"
    R: float = 0.4
    min_pt: float = 100.0
    max_abs_eta: float = 2.0

    def definition(self) -> "fastjet.JetDefinition":
        try:
            algo = _ALGOS[self.algorithm]
        except KeyError:
            raise ValueError(
                f"unknown algorithm {self.algorithm!r}; choose from {sorted(_ALGOS)}"
            ) from None
        return fastjet.JetDefinition(algo, self.R)


def as_particles(pt, eta, phi, mass=None, energy=None) -> ak.Array:
    """Build an awkward array of Lorentz vectors from per-event jagged kinematics.

    Provide either ``mass`` or ``energy``. Arrays must share the same jagged layout
    (one sublist per event).
    """
    if mass is None and energy is None:
        raise ValueError("provide either mass or energy")
    fields = {"pt": pt, "eta": eta, "phi": phi}
    fields["mass" if mass is not None else "E"] = mass if mass is not None else energy
    return ak.zip(fields, with_name="Momentum4D")


def four(v: ak.Array) -> ak.Array:
    """Strip a record down to exactly four momentum fields.

    ⚠️ Always do this before calling fastjet. A record carrying extra fields (charge,
    PID, vertex index, ...) can silently yield wrong jets, and the failure is data
    dependent so it often looks fine on a first test. Carry the extra fields separately
    and re-associate them afterwards by dR matching.
    """
    return ak.zip(
        {"pt": v.pt, "eta": v.eta, "phi": v.phi, "mass": v.mass},
        with_name="Momentum4D",
    )


def cluster(particles: ak.Array, config: JetConfig | None = None):
    """Cluster particles into jets.

    Parameters
    ----------
    particles : ak.Array
        Jagged array of ``Momentum4D`` records, one sublist per event.
    config : JetConfig, optional
        Jet definition and cuts. Defaults to anti-kt R=0.4, pT > 100 GeV, |eta| < 2.

    Returns
    -------
    jets : ak.Array
        Selected jets per event.
    constituents : ak.Array
        Constituents of each selected jet (needed for EEC).
    """
    config = config or JetConfig()
    seq = fastjet.ClusterSequence(four(particles), config.definition())

    jets = seq.inclusive_jets(min_pt=config.min_pt)
    constituents = seq.constituents(min_pt=config.min_pt)

    # Acceptance cut, applied identically to jets and their constituents
    keep = abs(jets.eta) < config.max_abs_eta
    return jets[keep], constituents[keep]
