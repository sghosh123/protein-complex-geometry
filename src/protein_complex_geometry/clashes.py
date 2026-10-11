"""Chunked inter-structure distance checks to avoid allocating an N-by-M matrix."""
from __future__ import annotations
from dataclasses import dataclass
from typing import Iterable
import numpy as np
from .pdbio import AtomRecord, PDBStructure


@dataclass(frozen=True)
class AtomPairDistance:
    atom_a: str
    atom_b: str
    distance_angstrom: float
    chain_a: str
    residue_a: str
    chain_b: str
    residue_b: str

    def as_dict(self) -> dict:
        return {"atom_a": self.atom_a, "atom_b": self.atom_b,
                "distance_angstrom": self.distance_angstrom,
                "chain_a": self.chain_a, "residue_a": self.residue_a,
                "chain_b": self.chain_b, "residue_b": self.residue_b}


def interstructure_pairs(structure_a: PDBStructure, structure_b: PDBStructure, *,
                          cutoff_angstrom: float = 2.0,
                          chain_a: str | None = None, chain_b: str | None = None,
                          heavy_only: bool = True, chunk_size: int = 512,
                          max_pairs: int | None = None) -> list[AtomPairDistance]:
    """List pairs at or below a distance cutoff, using Angstroms.

    The two structures are treated as separate partners. No bond exclusion is
    applied; use this for inter-partner clashes/contacts, not intramolecular
    clash validation. ``max_pairs`` can cap output size for dense overlaps.
    """
    if not np.isfinite(cutoff_angstrom) or cutoff_angstrom <= 0:
        raise ValueError("cutoff_angstrom must be a positive finite number")
    if chunk_size < 1:
        raise ValueError("chunk_size must be >= 1")
    if max_pairs is not None and max_pairs < 1:
        raise ValueError("max_pairs must be >= 1 when provided")
    aa = structure_a.selected_atoms(chain=chain_a, heavy_only=heavy_only)
    bb = structure_b.selected_atoms(chain=chain_b, heavy_only=heavy_only)
    if not aa or not bb:
        raise ValueError("atom selection is empty; check structures, chains, and heavy-atom filtering")
    xyz_a = np.asarray([a.xyz for a in aa])
    xyz_b = np.asarray([a.xyz for a in bb])
    cutoff2 = float(cutoff_angstrom) ** 2
    results: list[AtomPairDistance] = []
    # Block both dimensions so the temporary distance matrix is bounded by
    # chunk_size x chunk_size even when the second structure is very large.
    for start_a in range(0, len(aa), chunk_size):
        block_a = xyz_a[start_a:start_a + chunk_size]
        for start_b in range(0, len(bb), chunk_size):
            block_b = xyz_b[start_b:start_b + chunk_size]
            d2 = np.sum((block_a[:, None, :] - block_b[None, :, :]) ** 2, axis=2)
            hit_i, hit_j = np.where(d2 <= cutoff2)
            for i, j in zip(hit_i.tolist(), hit_j.tolist()):
                a = aa[start_a + i]; b = bb[start_b + j]
                results.append(AtomPairDistance(
                    a.label, b.label, float(np.sqrt(d2[i, j])), a.chain_id,
                    f"{a.residue_name}{a.residue_id}{a.insertion_code}", b.chain_id,
                    f"{b.residue_name}{b.residue_id}{b.insertion_code}"))
                if max_pairs is not None and len(results) >= max_pairs:
                    return sorted(results, key=lambda p: p.distance_angstrom)
    return sorted(results, key=lambda p: p.distance_angstrom)
