import pytest
from protein_complex_geometry.clashes import interstructure_pairs
from protein_complex_geometry.pdbio import PDBStructure

def test_pair_distances_cutoff_and_sort(tmp_path, write_pdb):
    a = write_pdb(tmp_path / "a.pdb", [(1, "CA", "ALA", "A", 1, 0, 0, 0, "C"),
                                       (2, "CB", "ALA", "A", 1, 10, 0, 0, "C")])
    b = write_pdb(tmp_path / "b.pdb", [(1, "CA", "GLY", "B", 8, 1.2, 0, 0, "C"),
                                       (2, "CB", "GLY", "B", 8, 12.2, 0, 0, "C")])
    pairs = interstructure_pairs(PDBStructure.read(a), PDBStructure.read(b), cutoff_angstrom=2.0)
    assert len(pairs) == 1
    assert pairs[0].distance_angstrom == pytest.approx(1.2)

def test_pair_distance_empty_selection_raises(tmp_path, write_pdb):
    a = write_pdb(tmp_path / "a.pdb", [(1, "CA", "ALA", "A", 1, 0, 0, 0, "C")])
    b = write_pdb(tmp_path / "b.pdb", [(1, "CA", "GLY", "B", 8, 1.2, 0, 0, "C")])
    with pytest.raises(ValueError, match="empty"):
        interstructure_pairs(PDBStructure.read(a), PDBStructure.read(b), chain_a="Z")

def test_pair_distance_invalid_cutoff(tmp_path, write_pdb):
    a = write_pdb(tmp_path / "a.pdb", [(1, "CA", "ALA", "A", 1, 0, 0, 0, "C")])
    b = write_pdb(tmp_path / "b.pdb", [(1, "CA", "GLY", "B", 8, 1.2, 0, 0, "C")])
    with pytest.raises(ValueError, match="positive finite"):
        interstructure_pairs(PDBStructure.read(a), PDBStructure.read(b), cutoff_angstrom=0)

def test_chunked_search_is_independent_of_chunk_size(tmp_path, write_pdb):
    a = write_pdb(tmp_path / "a.pdb", [
        (1, "CA", "ALA", "A", 1, 0, 0, 0, "C"),
        (2, "CB", "ALA", "A", 1, 10, 0, 0, "C"),
        (3, "N", "GLY", "A", 2, 15, 0, 0, "N"),
    ])
    b = write_pdb(tmp_path / "b.pdb", [
        (1, "CA", "SER", "B", 8, 1, 0, 0, "C"),
        (2, "CB", "SER", "B", 8, 11, 0, 0, "C"),
        (3, "O", "THR", "B", 9, 16, 0, 0, "O"),
    ])
    sa, sb = PDBStructure.read(a), PDBStructure.read(b)
    large_block = interstructure_pairs(sa, sb, cutoff_angstrom=2.0, chunk_size=512)
    small_block = interstructure_pairs(sa, sb, cutoff_angstrom=2.0, chunk_size=1)
    assert [(p.atom_a, p.atom_b, p.distance_angstrom) for p in large_block] == [
        (p.atom_a, p.atom_b, p.distance_angstrom) for p in small_block]
