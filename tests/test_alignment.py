import numpy as np
import pytest
from protein_complex_geometry.alignment import align_structures, matching_atom_pairs
from protein_complex_geometry.geometry import euler_rotation_matrix
from protein_complex_geometry.pdbio import PDBStructure

def _atoms(points, chain, residue_names=None):
    if residue_names is None: residue_names = ["ALA", "GLY", "SER", "THR"][:len(points)]
    return [(i+1, "CA", residue_names[i], chain, i+1, *map(float, p), "C") for i, p in enumerate(points)]

def test_align_pdbs_by_atom_residue_key(tmp_path, write_pdb, noncollinear_points):
    rotation = euler_rotation_matrix(33, -12, 70); shift = np.array([8., -2., 4.])
    reference_path = write_pdb(tmp_path / "ref.pdb", _atoms(noncollinear_points, "B"))
    mobile_path = write_pdb(tmp_path / "mob.pdb", _atoms(noncollinear_points @ rotation + shift, "C"))
    aligned, result = align_structures(PDBStructure.read(reference_path), PDBStructure.read(mobile_path),
                                      reference_chain="B", mobile_chain="C")
    assert result.n_atoms == 4
    assert result.rmsd_after_angstrom < 0.002
    np.testing.assert_allclose(aligned.coordinates, noncollinear_points, atol=1e-3)

def test_alignment_residue_filter(tmp_path, write_pdb, noncollinear_points):
    path1 = write_pdb(tmp_path / "ref.pdb", _atoms(noncollinear_points, "A"))
    path2 = write_pdb(tmp_path / "mob.pdb", _atoms(noncollinear_points, "B"))
    r, m = matching_atom_pairs(PDBStructure.read(path1), PDBStructure.read(path2),
                               reference_chain="A", mobile_chain="B", residue_ids={1, 3, 4})
    assert len(r) == len(m) == 3

def test_alignment_fails_when_no_correspondence(tmp_path, write_pdb):
    p1 = write_pdb(tmp_path / "ref.pdb", [(1, "CA", "ALA", "A", 1, 0, 0, 0, "C")])
    p2 = write_pdb(tmp_path / "mob.pdb", [(1, "CA", "GLY", "B", 2, 1, 0, 0, "C")])
    with pytest.raises(ValueError, match="no corresponding"):
        matching_atom_pairs(PDBStructure.read(p1), PDBStructure.read(p2))

def test_alignment_requires_at_least_three_matching_atoms(tmp_path, write_pdb):
    atoms = [(1, "CA", "ALA", "A", 1, 0, 0, 0, "C"),
             (2, "CA", "GLY", "A", 2, 1, 0, 0, "C")]
    ref = write_pdb(tmp_path / "ref.pdb", atoms)
    mob = write_pdb(tmp_path / "mob.pdb", [(a[0], a[1], a[2], "B", a[4], a[5]+1, a[6], a[7], a[8]) for a in atoms])
    with pytest.raises(ValueError, match="at least three"):
        align_structures(PDBStructure.read(ref), PDBStructure.read(mob))

def test_align_different_residue_ranges_and_residue_names(tmp_path, write_pdb, noncollinear_points):
    from protein_complex_geometry.geometry import euler_rotation_matrix
    rotation = euler_rotation_matrix(-22, 16, 118)
    shift = np.array([-5., 8., 2.])
    names_ref = ["ALA", "GLY", "SER", "THR"]
    names_mob = ["VAL", "LEU", "ASP", "PHE"]
    ref_atoms = [(i+1, "CA", names_ref[i], "B", 41+i, *map(float, p), "C")
                 for i, p in enumerate(noncollinear_points)]
    mob_xyz = noncollinear_points @ rotation + shift
    mob_atoms = [(i+1, "CA", names_mob[i], "C", 44+i, *map(float, p), "C")
                 for i, p in enumerate(mob_xyz)]
    ref_path = write_pdb(tmp_path / "reference.pdb", ref_atoms)
    mob_path = write_pdb(tmp_path / "mobile.pdb", mob_atoms)
    aligned, result = align_structures(
        PDBStructure.read(ref_path), PDBStructure.read(mob_path),
        reference_chain="B", mobile_chain="C",
        reference_residue_ids={41, 42, 43, 44}, mobile_residue_ids={44, 45, 46, 47},
        atom_name="CA")
    assert result.n_atoms == 4
    assert result.rmsd_after_angstrom < 0.002
    assert "ALA41" in result.matched_labels[0] and "VAL44" in result.matched_labels[0]
    np.testing.assert_allclose(aligned.coordinates, noncollinear_points, atol=1e-3)

def test_separate_residue_selectors_must_be_provided_in_pairs(tmp_path, write_pdb):
    path1 = write_pdb(tmp_path / "ref.pdb", [(1, "CA", "ALA", "A", 1, 0, 0, 0, "C")])
    path2 = write_pdb(tmp_path / "mob.pdb", [(1, "CA", "GLY", "B", 2, 1, 0, 0, "C")])
    with pytest.raises(ValueError, match="provide both"):
        matching_atom_pairs(PDBStructure.read(path1), PDBStructure.read(path2),
                            reference_residue_ids={1})
