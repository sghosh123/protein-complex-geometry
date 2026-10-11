import numpy as np
import pytest
from protein_complex_geometry.geometry import transform_from_euler
from protein_complex_geometry.pdbio import PDBStructure

def test_read_write_preserves_noncoordinate_lines(tmp_path, write_pdb):
    path = write_pdb(tmp_path / "a.pdb", [(1, "CA", "ALA", "A", 1, 1, 2, 3, "C")])
    structure = PDBStructure.read(path)
    output = tmp_path / "out.pdb"
    structure.transformed(transform_from_euler(translation=(3, 0, 0))).write(output)
    text = output.read_text()
    assert text.startswith("HEADER    SYNTHETIC TEST STRUCTURE")
    assert "TER\nEND\n" in text
    reread = PDBStructure.read(output)
    np.testing.assert_allclose(reread.coordinates, [[4, 2, 3]], atol=1e-3)

def test_rejects_empty_and_bad_pdb(tmp_path):
    empty = tmp_path / "empty.pdb"; empty.write_text("END\n")
    with pytest.raises(ValueError, match="no ATOM"):
        PDBStructure.read(empty)
    malformed = tmp_path / "malformed.pdb"; malformed.write_text("ATOM      1  CA  ALA A   1\n")
    with pytest.raises(ValueError, match="fewer than 54"):
        PDBStructure.read(malformed)

def test_rejects_multimodel_pdb(tmp_path, write_pdb):
    path = write_pdb(tmp_path / "one.pdb", [(1, "CA", "ALA", "A", 1, 0, 0, 0, "C")])
    text = path.read_text()
    path.write_text("MODEL        1\n" + text + "MODEL        2\n" + text)
    with pytest.raises(ValueError, match="multi-model"):
        PDBStructure.read(path)

def test_coordinate_count_and_nonfinite_validation(tmp_path, write_pdb):
    path = write_pdb(tmp_path / "a.pdb", [(1, "CA", "ALA", "A", 1, 0, 0, 0, "C")])
    structure = PDBStructure.read(path)
    with pytest.raises(ValueError, match="shape"):
        structure.with_coordinates([[0, 0, 0], [1, 1, 1]])
    with pytest.raises(ValueError, match="non-finite"):
        structure.with_coordinates([[float("nan"), 0, 0]])

def test_element_fallback_identifies_hydrogen_and_copper(tmp_path, write_pdb):
    path = write_pdb(tmp_path / "elements.pdb", [
        (1, "1H", "CYS", "A", 1, 0, 0, 0, "H"),
        (2, "CU1", "CU1", "B", 2, 2, 0, 0, "CU"),
        (3, "CA", "ALA", "A", 3, 4, 0, 0, "C"),
    ])
    lines = path.read_text().splitlines(keepends=True)
    # Clear element columns and use standard left justification for two-letter Cu atom names.
    lines = [(line[:12] + "CU1 " + line[16:76] + "  " + line[78:])
             if line.startswith("ATOM") and line[17:20].strip() == "CU1" else
             (line[:76] + "  " + line[78:] if line.startswith(("ATOM", "HETATM")) else line)
             for line in lines]
    path.write_text("".join(lines))
    structure = PDBStructure.read(path)
    assert [atom.element for atom in structure.atoms] == ["H", "CU", "C"]
    assert [atom.atom_name for atom in structure.selected_atoms(heavy_only=True)] == ["CU1", "CA"]


def test_writer_rejects_coordinates_that_overflow_pdb_fields(tmp_path, write_pdb):
    path = write_pdb(tmp_path / "a.pdb", [(1, "CA", "ALA", "A", 1, 0, 0, 0, "C")])
    structure = PDBStructure.read(path)
    with pytest.raises(ValueError, match="fixed-width PDB field"):
        structure.with_coordinates([[10000.0, 0, 0]])

def test_transform_rejects_stale_anisou_records(tmp_path, write_pdb):
    path = write_pdb(tmp_path / "anisou.pdb", [(1, "CA", "ALA", "A", 1, 0, 0, 0, "C")])
    text = path.read_text()
    atom_line = next(line for line in text.splitlines(keepends=True) if line.startswith("ATOM"))
    path.write_text(text.replace(atom_line, atom_line + "ANISOU    1  CA  ALA A   1     1000   1000   1000      0      0      0       C\n"))
    structure = PDBStructure.read(path)
    with pytest.raises(ValueError, match="ANISOU"):
        structure.with_coordinates(structure.coordinates + 1.0)

def test_cryst1_frame_warning_is_explicit(tmp_path, write_pdb):
    path = write_pdb(tmp_path / "crystal.pdb", [(1, "CA", "ALA", "A", 1, 0, 0, 0, "C")])
    path.write_text("CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1           1\n" + path.read_text())
    structure = PDBStructure.read(path)
    assert any("CRYST1" in warning for warning in structure.coordinate_frame_warnings())
