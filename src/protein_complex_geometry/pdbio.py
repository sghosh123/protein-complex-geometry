"""Conservative PDB coordinate I/O that preserves non-coordinate records."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable
import numpy as np

_ELEMENT_SYMBOLS = {
    "H", "HE", "LI", "BE", "B", "C", "N", "O", "F", "NE", "NA", "MG", "AL", "SI", "P", "S", "CL", "AR",
    "K", "CA", "SC", "TI", "V", "CR", "MN", "FE", "CO", "NI", "CU", "ZN", "GA", "GE", "AS", "SE", "BR", "KR",
    "RB", "SR", "Y", "ZR", "NB", "MO", "TC", "RU", "RH", "PD", "AG", "CD", "IN", "SN", "SB", "TE", "I", "XE",
    "CS", "BA", "LA", "CE", "PR", "ND", "PM", "SM", "EU", "GD", "TB", "DY", "HO", "ER", "TM", "YB", "LU",
    "HF", "TA", "W", "RE", "OS", "IR", "PT", "AU", "HG", "TL", "PB", "BI", "PO", "AT", "RN", "FR", "RA",
    "AC", "TH", "PA", "U", "NP", "PU", "AM", "CM", "BK", "CF", "ES", "FM", "MD", "NO", "LR", "RF", "DB",
    "SG", "BH", "HS", "MT", "DS", "RG", "CN", "NH", "FL", "MC", "LV", "TS", "OG"
}
_STANDARD_AA = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"}

def _infer_element(atom_field: str, atom_name: str, residue_name: str, record_type: str) -> str:
    """Best-effort PDB element fallback when columns 77-78 are blank."""
    name = atom_name.strip().upper()
    if not name:
        return "?"
    # CA in a standard amino-acid ATOM record means alpha carbon, not calcium.
    if record_type == "ATOM" and residue_name.upper() in _STANDARD_AA and name == "CA":
        return "C"
    # Leading-digit atom names (e.g. 1HG) begin with the element after the digit.
    letters = name[1:] if name[0].isdigit() else name
    if not letters:
        return "?"
    # PDB convention: two-letter elements are generally left-justified in the atom field.
    if atom_field and not atom_field[0].isspace() and len(letters) >= 2 and letters[:2] in _ELEMENT_SYMBOLS:
        return letters[:2]
    return letters[0]


@dataclass(frozen=True)
class AtomRecord:
    line_index: int
    record_type: str
    serial: int | None
    atom_name: str
    altloc: str
    residue_name: str
    chain_id: str
    residue_id: int
    insertion_code: str
    x: float
    y: float
    z: float
    element: str
    occupancy: float | None = None

    @property
    def xyz(self) -> np.ndarray:
        return np.array([self.x, self.y, self.z], dtype=float)

    @property
    def residue_key(self) -> tuple[str, int, str, str]:
        return (self.chain_id, self.residue_id, self.insertion_code, self.residue_name)

    @property
    def atom_key(self) -> tuple[str, int, str, str, str]:
        return (*self.residue_key, self.atom_name)

    @property
    def label(self) -> str:
        return (f"{self.chain_id or '_'}:{self.residue_name}{self.residue_id}"
                f"{self.insertion_code}:{self.atom_name}#{self.serial if self.serial is not None else '?'}")


class PDBStructure:
    """Single-model PDB with fixed-width coordinate updates and preserved records."""
    def __init__(self, path: str | Path, lines: list[str], atoms: list[AtomRecord]):
        self.path = Path(path)
        self.lines = lines
        self.atoms = atoms

    @classmethod
    def read(cls, path: str | Path) -> "PDBStructure":
        path = Path(path)
        try:
            lines = path.read_text(encoding="utf-8", errors="replace").splitlines(keepends=True)
        except OSError as exc:
            raise OSError(f"could not read PDB file {path}: {exc}") from exc
        if not lines:
            raise ValueError(f"PDB file is empty: {path}")
        model_records = sum(line[:6].strip() == "MODEL" for line in lines)
        if model_records > 1:
            raise ValueError("multi-model PDB files are not supported; split the model first")
        atoms: list[AtomRecord] = []
        for idx, line in enumerate(lines):
            rec = line[:6].strip()
            if rec not in {"ATOM", "HETATM"}:
                continue
            if len(line.rstrip("\r\n")) < 54:
                raise ValueError(f"malformed coordinate record at {path}:{idx+1}: fewer than 54 columns")
            try:
                serial = int(line[6:11].strip()) if line[6:11].strip() else None
                atom_name = line[12:16].strip()
                residue_name = line[17:20].strip()
                chain_id = line[21:22].strip()
                residue_id = int(line[22:26].strip())
                insertion_code = line[26:27].strip()
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                occupancy = float(line[54:60]) if len(line) >= 60 and line[54:60].strip() else None
            except (ValueError, IndexError) as exc:
                raise ValueError(f"could not parse coordinate record at {path}:{idx+1}") from exc
            if not np.isfinite([x, y, z]).all():
                raise ValueError(f"non-finite coordinate at {path}:{idx+1}")
            element = line[76:78].strip().upper() if len(line) >= 78 else ""
            if not element:
                element = _infer_element(line[12:16], atom_name, residue_name, rec)
            atoms.append(AtomRecord(idx, rec, serial, atom_name, line[16:17].strip(),
                                    residue_name, chain_id, residue_id, insertion_code,
                                    x, y, z, element, occupancy))
        if not atoms:
            raise ValueError(f"no ATOM/HETATM coordinate records found in {path}")
        return cls(path, lines, atoms)

    @property
    def coordinates(self) -> np.ndarray:
        return np.asarray([[a.x, a.y, a.z] for a in self.atoms], dtype=float)

    def coordinate_frame_warnings(self) -> list[str]:
        """Describe preserved PDB frame records that may need review after transforms."""
        records = {line[:6].strip() for line in self.lines}
        warnings = []
        if "CRYST1" in records:
            warnings.append("CRYST1 unit-cell metadata is preserved unchanged; review/remove/update it for the transformed coordinate frame")
        if records & {"SCALE1", "SCALE2", "SCALE3", "ORIGX1", "ORIGX2", "ORIGX3", "MTRIX1", "MTRIX2", "MTRIX3"}:
            warnings.append("SCALE/ORIGX/MTRIX frame records are preserved unchanged and may not match transformed coordinates")
        return warnings

    def selected_atoms(self, *, chain: str | None = None,
                       atom_names: Iterable[str] | None = None,
                       residue_ids: set[int] | None = None,
                       heavy_only: bool = False,
                       include_alternates: bool = False) -> list[AtomRecord]:
        names = {n.strip().upper() for n in atom_names} if atom_names is not None else None
        out = []
        for atom in self.atoms:
            if chain is not None and atom.chain_id != chain:
                continue
            if names is not None and atom.atom_name.upper() not in names:
                continue
            if residue_ids is not None and atom.residue_id not in residue_ids:
                continue
            if heavy_only and atom.element in {"H", "D"}:
                continue
            if not include_alternates and atom.altloc not in {"", "A"}:
                continue
            out.append(atom)
        return out

    def with_coordinates(self, coordinates: Iterable[Iterable[float]]) -> "PDBStructure":
        tensor_records = {line[:6].strip() for line in self.lines} & {"ANISOU", "SIGATM", "SIGUIJ"}
        if tensor_records:
            names = ", ".join(sorted(tensor_records))
            raise ValueError(f"cannot transform PDB with coordinate-associated {names} records without updating their tensors")
        xyz = np.asarray(coordinates, dtype=float)
        if xyz.shape != (len(self.atoms), 3):
            raise ValueError(f"coordinate array must have shape ({len(self.atoms)}, 3); got {xyz.shape}")
        if not np.isfinite(xyz).all():
            raise ValueError("coordinates contain non-finite values")
        formatted = [[f"{value:8.3f}" for value in row] for row in xyz]
        if any(len(field) != 8 for row in formatted for field in row):
            raise ValueError("coordinates exceed the fixed-width PDB field; use a format supporting larger coordinates")
        new_lines = list(self.lines)
        for atom, fields in zip(self.atoms, formatted):
            old = new_lines[atom.line_index]
            body = old.rstrip("\r\n")
            ending = old[len(body):]
            body = body.ljust(54)
            new_lines[atom.line_index] = body[:30] + "".join(fields) + body[54:] + ending
        return PDBStructure(self.path, new_lines, [
            AtomRecord(a.line_index, a.record_type, a.serial, a.atom_name, a.altloc,
                       a.residue_name, a.chain_id, a.residue_id, a.insertion_code,
                       float(x), float(y), float(z), a.element, a.occupancy)
            for a, (x, y, z) in zip(self.atoms, xyz)
        ])

    def transformed(self, transform) -> "PDBStructure":
        return self.with_coordinates(transform.apply(self.coordinates))

    def write(self, path: str | Path) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("".join(self.lines), encoding="utf-8")
        return path
