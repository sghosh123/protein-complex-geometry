import json
from protein_complex_geometry.cli import main

def _atoms(points, chain):
    names = ["ALA", "GLY", "SER", "THR"]
    return [(i+1, "CA", names[i], chain, i+1, *map(float, p), "C") for i, p in enumerate(points)]

def test_transform_cli_writes_pdb_and_report(tmp_path, write_pdb):
    inp = write_pdb(tmp_path / "in.pdb", [(1, "CA", "ALA", "A", 1, 1, 0, 0, "C")])
    out, report = tmp_path / "out.pdb", tmp_path / "transform.json"
    assert main(["transform", "--input", str(inp), "--output", str(out), "--rotate-z", "90",
                 "--origin", "0,0,0", "--translate", "1,2,3", "--report", str(report)]) == 0
    data = report.read_text(); doc = json.loads(data)
    assert doc["translation_angstrom"] == [1., 2., 3.]
    assert out.exists()

def test_clashes_cli_writes_header_even_without_hits(tmp_path, write_pdb):
    a = write_pdb(tmp_path / "a.pdb", [(1, "CA", "ALA", "A", 1, 0, 0, 0, "C")])
    b = write_pdb(tmp_path / "b.pdb", [(1, "CA", "GLY", "B", 8, 100, 0, 0, "C")])
    out = tmp_path / "pairs.csv"
    assert main(["clashes", "--structure-a", str(a), "--structure-b", str(b),
                 "--cutoff", "2", "--output", str(out)]) == 0
    assert "distance_angstrom" in out.read_text()

def test_align_cli_supports_different_residue_ranges_and_logs_mapping(tmp_path, write_pdb, noncollinear_points):
    import numpy as np
    from protein_complex_geometry.geometry import euler_rotation_matrix
    reference_names = ["ALA", "GLY", "SER", "THR"]
    mobile_names = ["VAL", "LEU", "ASP", "PHE"]
    rotation = euler_rotation_matrix(12, -37, 81)
    mobile_xyz = noncollinear_points @ rotation + np.array([3., -6., 5.])
    ref_rows = [(i+1, "CA", reference_names[i], "B", 41+i, *map(float, point), "C")
                for i, point in enumerate(noncollinear_points)]
    mob_rows = [(i+1, "CA", mobile_names[i], "C", 44+i, *map(float, point), "C")
                for i, point in enumerate(mobile_xyz)]
    ref = write_pdb(tmp_path / "ref.pdb", ref_rows)
    mob = write_pdb(tmp_path / "mob.pdb", mob_rows)
    output, report = tmp_path / "aligned.pdb", tmp_path / "alignment.json"
    assert main(["align", "--reference", str(ref), "--mobile", str(mob),
                 "--reference-chain", "B", "--mobile-chain", "C", "--atom-name", "CA",
                 "--reference-residues", "41-44", "--mobile-residues", "44-47",
                 "--output", str(output), "--report", str(report)]) == 0
    data = json.loads(report.read_text())
    assert data["n_atoms"] == 4
    assert data["rmsd_after_angstrom"] < 0.002
    assert "ALA41" in data["matched_atoms_reference_to_mobile"][0]
    assert "VAL44" in data["matched_atoms_reference_to_mobile"][0]
