from pathlib import Path
import numpy as np
import pytest

@pytest.fixture
def write_pdb(tmp_path):
    def _write(path, atoms):
        """atoms: iterable of (serial, name, resname, chain, resid, x, y, z, element)."""
        lines = ["HEADER    SYNTHETIC TEST STRUCTURE\n"]
        for serial, name, resname, chain, resid, x, y, z, element in atoms:
            line = (f"ATOM  {serial:5d} {name:>4s} {resname:>3s} {chain:1s}{resid:4d}    "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{20.00:6.2f}          {element:>2s}  \n")
            lines.append(line)
        lines.append("TER\nEND\n")
        path = Path(path); path.write_text("".join(lines), encoding="utf-8")
        return path
    return _write

@pytest.fixture
def noncollinear_points():
    return np.array([[0., 0., 0.], [1., 0., 0.], [0., 2., 0.], [0., 0., 3.]])
