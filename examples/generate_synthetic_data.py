"""Regenerate the tiny, public-domain-style synthetic PDBs used in examples."""
from pathlib import Path
import numpy as np
from protein_complex_geometry.geometry import euler_rotation_matrix

ROOT = Path(__file__).parent / "data"

def write_pdb(path: Path, rows):
    lines = ["HEADER    SYNTHETIC GEOMETRY DEMO; NOT A BIOLOGICAL STRUCTURE\n"]
    for serial, atom, resname, chain, resid, xyz, element in rows:
        atom_field = (f" {atom:<3s}" if len(element) == 1 else f"{atom:<4s}")
        x, y, z = xyz
        lines.append(f"ATOM  {serial:5d} {atom_field} {resname:>3s} {chain:1s}{resid:4d}    "
                     f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{20.0:6.2f}          {element:>2s}  \n")
    lines.extend(["TER\n", "END\n"])
    path.write_text("".join(lines), encoding="utf-8")

def residues(chain, xyz, *, first_residue=30, names=None):
    names = names or ["ALA", "GLY", "SER", "THR"]
    return [(i + 1, "CA", names[i], chain, first_residue + i, p, "C") for i, p in enumerate(xyz)]

def main():
    ROOT.mkdir(parents=True, exist_ok=True)
    reference = np.array([[0., 0., 0.], [1., 0., 0.], [0., 2., 0.], [0., 0., 3.]])
    mobile = reference @ euler_rotation_matrix(z_degrees=90) + np.array([8., -2., 4.])
    write_pdb(ROOT / "reference.pdb", residues("B", reference, first_residue=30))
    write_pdb(ROOT / "mobile.pdb", residues("C", mobile, first_residue=44,
                                               names=["VAL", "LEU", "ASP", "PHE"]))
    write_pdb(ROOT / "partner_a.pdb", [
        (1, "CA", "ALA", "A", 1, (0., 0., 0.), "C"),
        (2, "CB", "ALA", "A", 1, (0., 1.5, 0.), "C"),
        (3, "SG", "CYS", "A", 2, (0., 4.0, 0.), "S"),
    ])
    write_pdb(ROOT / "partner_b.pdb", [
        (1, "CA", "GLY", "B", 8, (1.5, 0., 0.), "C"),
        (2, "CB", "GLY", "B", 8, (0., 3.4, 0.), "C"),
        (3, "CU1", "CU1", "B", 9, (3.1, 4.0, 0.), "CU"),
    ])
    print(f"Wrote synthetic PDB examples into {ROOT}")

if __name__ == "__main__":
    main()
