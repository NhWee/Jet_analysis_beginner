"""Jet clustering.

Thin wrapper over scikit-hep ``fastjet`` (awkward/array-oriented interface) so the rest of
the analysis never touches fastjet internals directly. Works identically on Pythia truth
particles and on experimental constituents (PFCandidates / tracks / topoclusters).

Two fastjet behaviours are handled here, so callers cannot trip over them:

* **Non-contiguous input gives silently wrong jets in fastjet <= 3.5.1.4.** Those versions
  misread an awkward array whose sublists are not stored back to back in memory -- what
  you get after selecting or reordering the constituents ``ClusterSequence.constituents()``
  returns (a ``ListArray`` over an ``IndexedArray``). Nothing is raised: reclustering such
  constituents into subjets missed the jet by hundreds of GeV in this project (README,
  "Known trap"). Fixed in fastjet 3.5.1.5. Every call into fastjet below goes through
  :func:`clusterable`, which makes a packed copy, so the result is right on any version.
* **Jets come back in no particular order.** :func:`cluster` sorts them by descending pT,
  so ``jets[:, 0]`` is the leading jet.

Extra per-particle fields (charge, PID, vertex index, ...) are harmless and are carried
through to the constituents, so nothing has to be stripped beforehand or re-matched by dR
afterwards.
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


def clusterable(particles: ak.Array) -> ak.Array:
    """Return a packed copy of ``particles`` that every fastjet version reads correctly.

    All fields are kept. Use this before any direct call to ``fastjet.ClusterSequence``;
    the functions in this module already do.
    """
    return ak.to_packed(particles)


def four(v: ak.Array) -> ak.Array:
    """Strip a record down to its four momentum fields ``pt, eta, phi, mass`` (packed).

    Clustering does not need this -- extra fields are harmless and :func:`cluster` keeps
    them on the constituents -- but it is handy when a plain four-vector array is wanted.
    """
    return ak.to_packed(ak.zip(
        {"pt": v.pt, "eta": v.eta, "phi": v.phi, "mass": v.mass},
        with_name="Momentum4D",
    ))


def cluster(particles: ak.Array, config: JetConfig | None = None):
    """Cluster particles into jets.

    Parameters
    ----------
    particles : ak.Array
        Jagged array of ``Momentum4D`` records, one sublist per event. Extra fields
        (charge, PID, ...) are allowed and are carried through to the constituents.
    config : JetConfig, optional
        Jet definition and cuts. Defaults to anti-kt R=0.4, pT > 100 GeV, |eta| < 2.

    Returns
    -------
    jets : ak.Array
        Selected jets per event, sorted by descending pT.
    constituents : ak.Array
        Constituents of each selected jet, ``[event][jet][particle]``, with every field
        of the input records (needed for EEC, and for charged-only observables).
    """
    config = config or JetConfig()
    seq = fastjet.ClusterSequence(clusterable(particles), config.definition())

    jets = seq.inclusive_jets(min_pt=config.min_pt)
    constituents = seq.constituents(min_pt=config.min_pt)

    # fastjet returns jets in no particular order; make jets[:, 0] the leading jet
    order = ak.argsort(jets.pt, axis=-1, ascending=False)
    jets, constituents = jets[order], constituents[order]

    # Acceptance cut, applied identically to jets and their constituents
    keep = abs(jets.eta) < config.max_abs_eta
    return ak.to_packed(jets[keep]), ak.to_packed(constituents[keep])


def subjets(constituents: ak.Array, algorithm: str = "kt", R: float = 0.2) -> ak.Array:
    """Recluster each jet's constituents into inclusive subjets, ``[jet][subjet]``.

    The building block for trimming, filtering and Lund-plane studies. ``constituents``
    must be a ``[jet][particle]`` array, e.g. ``ak.flatten(constituents, axis=1)`` of the
    output of :func:`cluster`. This is the call that went wrong on fastjet <= 3.5.1.4
    without :func:`clusterable`.
    """
    jetdef = JetConfig(algorithm=algorithm, R=R).definition()
    seq = fastjet.ClusterSequence(clusterable(constituents), jetdef)
    return seq.inclusive_jets()
