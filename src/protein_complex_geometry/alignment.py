"""PDB atom matching and structure superposition."""
from __future__ import annotations
from dataclasses import dataclass
from typing import Iterable
import numpy as np
from .geometry import RigidTransform, kabsch_transform, rmsd
from .pdbio import AtomRecord, PDBStructure


@dataclass(frozen=True)
class AlignmentResult:
    transform: RigidTransform
    rmsd_before_angstrom: float
    rmsd_after_angstrom: float
    n_atoms: int
    matched_labels: tuple[str, ...]

    def as_dict(self) -> dict:
        return {"n_atoms": self.n_atoms,
                "rmsd_before_angstrom": self.rmsd_before_angstrom,
                "rmsd_after_angstrom": self.rmsd_after_angstrom,
                "transform": self.transform.as_dict(),
                "matched_atoms_reference_to_mobile": list(self.matched_labels)}


def _selected_atom_map(structure: PDBStructure, *, chain: str | None,
                        residue_ids: set[int] | None, atom_name: str) -> dict[tuple, AtomRecord]:
    """Select one atom per residue/insertion/atom key, preferring blank altloc then A."""
    found: dict[tuple, AtomRecord] = {}
    for atom in structure.atoms:
        if chain is not None and atom.chain_id != chain:
            continue
        if residue_ids is not None and atom.residue_id not in residue_ids:
            continue
        if atom.atom_name.upper() != atom_name.upper() or atom.altloc not in {"", "A"}:
            continue
        # Residue name is deliberately not part of the key: equivalent backbone
        # atoms can be aligned across mutations or homologous protein chains.
        key = (atom.residue_id, atom.insertion_code, atom.atom_name.upper())
        previous = found.get(key)
        if previous is None:
            found[key] = atom
        elif previous.altloc != "" and atom.altloc == "":
            found[key] = atom
        elif previous.altloc == "" and atom.altloc != "":
            continue
        else:
            raise ValueError(f"duplicate atom mapping for residue/atom key {key} in {structure.path}")
    return found


def matching_atom_pairs(reference: PDBStructure, mobile: PDBStructure, *,
                        reference_chain: str | None = None,
                        mobile_chain: str | None = None,
                        residue_ids: Iterable[int] | None = None,
                        reference_residue_ids: Iterable[int] | None = None,
                        mobile_residue_ids: Iterable[int] | None = None,
                        atom_name: str = "CA") -> tuple[list[AtomRecord], list[AtomRecord]]:
    """Match alignment atoms between two structures.

    By default, atoms are matched by residue number, insertion code, and atom
    name. The residue name is not required to match, which permits C-alpha fits
    across mutations/homologs. If both ``reference_residue_ids`` and
    ``mobile_residue_ids`` are provided, selected atoms are paired in ascending
    residue order; this supports corresponding but differently numbered regions.
    The researcher is responsible for ensuring that these ordered selections are
    biologically corresponding. The legacy ``residue_ids`` option applies the
    same residue-number filter to both structures and cannot be combined with
    the separate selectors.
    """
    if not atom_name.strip():
        raise ValueError("atom_name cannot be empty")
    has_ref_specific = reference_residue_ids is not None
    has_mob_specific = mobile_residue_ids is not None
    if has_ref_specific != has_mob_specific:
        raise ValueError("provide both reference_residue_ids and mobile_residue_ids for ordered pairing")
    if residue_ids is not None and (has_ref_specific or has_mob_specific):
        raise ValueError("residue_ids cannot be combined with separate reference/mobile residue selectors")

    if has_ref_specific:
        ref_map = _selected_atom_map(reference, chain=reference_chain,
                                     residue_ids=set(reference_residue_ids), atom_name=atom_name)
        mob_map = _selected_atom_map(mobile, chain=mobile_chain,
                                     residue_ids=set(mobile_residue_ids), atom_name=atom_name)
        sort_key = lambda key: (key[0], key[1], key[2])
        ref_keys = sorted(ref_map, key=sort_key)
        mob_keys = sorted(mob_map, key=sort_key)
        if len(ref_keys) != len(mob_keys):
            raise ValueError(f"ordered selections contain different atom counts: reference={len(ref_keys)}, mobile={len(mob_keys)}")
        if not ref_keys:
            raise ValueError("no corresponding atoms found; check chains, residue selections, and atom name")
        return [ref_map[k] for k in ref_keys], [mob_map[k] for k in mob_keys]

    residue_filter = set(residue_ids) if residue_ids is not None else None
    ref_map = _selected_atom_map(reference, chain=reference_chain,
                                 residue_ids=residue_filter, atom_name=atom_name)
    mob_map = _selected_atom_map(mobile, chain=mobile_chain,
                                 residue_ids=residue_filter, atom_name=atom_name)
    common = sorted(set(ref_map) & set(mob_map))
    if not common:
        raise ValueError("no corresponding atoms found; check chain IDs, residue numbering, and atom name")
    return [ref_map[k] for k in common], [mob_map[k] for k in common]


def align_structures(reference: PDBStructure, mobile: PDBStructure, *,
                     reference_chain: str | None = None,
                     mobile_chain: str | None = None,
                     residue_ids: Iterable[int] | None = None,
                     reference_residue_ids: Iterable[int] | None = None,
                     mobile_residue_ids: Iterable[int] | None = None,
                     atom_name: str = "CA") -> tuple[PDBStructure, AlignmentResult]:
    ref_atoms, mob_atoms = matching_atom_pairs(
        reference, mobile, reference_chain=reference_chain, mobile_chain=mobile_chain,
        residue_ids=residue_ids, reference_residue_ids=reference_residue_ids,
        mobile_residue_ids=mobile_residue_ids, atom_name=atom_name)
    ref_xyz = np.asarray([a.xyz for a in ref_atoms])
    mob_xyz = np.asarray([a.xyz for a in mob_atoms])
    before = rmsd(mob_xyz, ref_xyz)
    transform, after = kabsch_transform(mob_xyz, ref_xyz)
    aligned = mobile.transformed(transform)
    mapped_labels = tuple(f"{ref.label} -> {mob.label}" for ref, mob in zip(ref_atoms, mob_atoms))
    result = AlignmentResult(transform, before, after, len(ref_atoms), mapped_labels)
    return aligned, result
