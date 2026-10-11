# protein-complex-geometry

**Reusable, testable geometry utilities for preparing and checking protein–protein complex structures.**

This package generalizes recurring structure-preparation tasks from the research workspace into a small Python toolkit: superpose a mobile protein using selected Cα atoms, apply rigid-body rotations/translations to PDB coordinates, and detect close inter-structure atom pairs before simulation or visualization.

## Capabilities

- **Rigid-body alignment:** Kabsch fit, with reflection disallowed, using explicitly matched residue/atom keys.
- **Coordinate transforms:** rotations around x/y/z axes in degrees and translations in Å; rotation about a user-specified pivot or the coordinate centroid.
- **Inter-structure contacts/clashes:** chunked heavy-atom pair distances with residue and atom labels, CSV output, configurable cutoff, chain filters, and optional output cap.
- **Structure provenance:** optional JSON reports contain fit RMSD, number of matched atoms, transform matrix, input/output names, and selected fit settings.
- **Conservative file handling:** fixed-width PDB coordinate fields are updated while other PDB lines are preserved.

## Install

```bash
python -m pip install -e '.[test]'
pytest -q
pcg --help
```

Requires Python 3.10+ and NumPy. The core package does not require MDAnalysis, OpenMM, VMD, or Biopython.

## Quick start

### 1. Align one structure onto another

```bash
pcg align \
  --reference reference.pdb \
  --mobile mobile.pdb \
  --reference-chain B \
  --mobile-chain C \
  --atom-name CA \
  --residues '30-33,41-46' \
  --output mobile_aligned.pdb \
  --report alignment.json
```

The residue numbers are assumed to correspond between the two structures. The fit uses only the matched atoms; every coordinate record in the mobile PDB is transformed. Inspect the reported matched-atom count and post-fit RMSD before proceeding.

For different numbering or equivalent segments in different proteins, select each range independently and pair them in ascending residue order (you are responsible for confirming the correspondence):

```bash
pcg align --reference reference.pdb --mobile mobile.pdb \
  --reference-chain B --mobile-chain C --atom-name CA \
  --reference-residues '41-46' --mobile-residues '44-49' \
  --output mobile_aligned.pdb --report alignment.json
```

### 2. Rotate and translate a structure

```bash
pcg transform \
  --input partner.pdb \
  --rotate-x 20 --rotate-y 20 --rotate-z 0 \
  --origin centroid \
  --translate '-41,5,-25' \
  --output partner_placed.pdb \
  --report transform.json
```

Rotations are applied x, then y, then z. `--origin centroid` uses the centroid of all parsed ATOM/HETATM records. To rotate around the Cartesian origin, use `--origin '0,0,0'`. The pivot remains fixed by the rotation itself, then the explicit translation is added.

### 3. Check two placed partners for close contacts

```bash
pcg clashes \
  --structure-a partner_a.pdb \
  --structure-b partner_b.pdb \
  --cutoff 2.0 \
  --output close_pairs.csv
```

This calculates **inter-structure** pairs only. The default cutoff is 2.0 Å, and hydrogens/deuteriums are excluded where their element is identified. Increase the cutoff (for example, 4.5 Å) to report broader contact pairs. This tool is not a full intramolecular clash validator and does not apply force-field radii.

## Geometry conventions

- Coordinates in PDB files are in **Ångström (Å)**.
- Rotation angles are in **degrees**.
- Matrices use row-vector application: `xyz_new = xyz @ rotation + translation`.
- The alignment report contains the post-fit RMSD in Å.
- Alignment requires at least three non-collinear, corresponding points; it will not silently fit a reflection.
- Atom matching uses residue number, insertion code, residue name, and atom name. Chain IDs are specified independently to allow equivalent chains with different IDs.

Read [`docs/methods-and-limitations.md`](docs/methods-and-limitations.md) before production use. The current synthetic test record is in [`docs/validation-record.md`](docs/validation-record.md).

## Reproducible demo

After installation, run `bash examples/run_demo.sh`. The small files in `examples/data/` are synthetic, non-biological coordinates and exist only to demonstrate the commands end to end. Output is written to `demo_output/`.

## Repository map

```text
src/protein_complex_geometry/  geometry, PDB I/O, fitting, clash routines, CLI
 tests/                         synthetic numerical and CLI tests
 examples/                      reproducible command examples
 docs/                          assumptions, limitations, provenance
```

## Citation and licensing

A recommended citation format is provided in `CITATION.cff`. The repository license is intentionally marked unresolved; choose and review a license before public release. See `LICENSE-STATUS.md`.
