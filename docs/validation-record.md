# Validation record — protein-complex-geometry 0.1.0a1

Validation date: 2026-10-11

## Automated checks

- `pytest -q`: **30 passed**.
- Source compilation check: passed.
- `pcg --help`: passed.
- Python wheel built successfully with the locally available build tools.
- The wheel installed into an isolated target directory and its CLI was exercised using the bundled synthetic PDB files.

## End-to-end synthetic demonstration

- Alignment: 4 C-alpha atom pairs; RMSD changed from 8.6891 Å before fitting to below displayed precision after fitting (synthetic coordinates were generated from a known rigid transform).
- Transform: the aligned PDB was rotated and written with a JSON transform report.
- Inter-structure close-pair scan: 3 pairs at or below 2.0 Å; the minimum reported distance was 0.600 Å.
- The demo intentionally pairs reference residues 30–33 with differently numbered mobile residues 44–47.

## What this does not establish

These checks use synthetic small PDBs and known rigid transforms. They do not validate the default selections, residue correspondences, copper coordination, biological plausibility, crystallographic metadata handling for any specific input, or downstream simulation behavior for a research structure. Review the transform and fit reports, inspect resulting structures visually, and validate selections against your original systems.

The project remains an alpha release. License choice and author/citation metadata must be resolved before publishing for external reuse.
