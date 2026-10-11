# Methods, assumptions, and limitations

## Kabsch superposition

Given corresponding coordinate rows `P` (mobile) and `Q` (reference), each set is centered and the covariance matrix is decomposed with SVD. A proper rotation is selected by correcting the determinant so that reflections are not allowed. The transform is applied to all ATOM/HETATM records in the mobile PDB, not just the fitting atoms.

The default CLI matches atoms by `(residue number, insertion code, atom name)` after applying the requested chain filters; residue names need not match, which permits C-alpha fits across mutations. For different numbering (for example, reference residues 41–46 against mobile residues 44–49), supply both `--reference-residues` and `--mobile-residues`. These selected atoms are paired in ascending residue order, so the researcher must verify the correspondence. The implementation does **not** infer sequence correspondence, chain permutations, missing-residue mapping, or structural equivalence. Inspect the reported matched-atom pairs, count, and fit RMSD. The bundled synthetic demo intentionally uses reference residues 30–33 and mobile residues 44–47 to exercise ordered pairing.

## Coordinate transforms

Rotations use active right-handed rotations. Euler rotations apply x, then y, then z. A transform with pivot `o`, rotation `R`, and user shift `t` is represented under the row-vector convention as `x_new = x @ R + (o - o @ R + t)`. Thus the pivot is invariant under rotation and is moved only by the user-specified shift.

All PDB coordinates and shifts are in Å; angles are in degrees. Avoid mixing a VMD view transform with a coordinate transform: this package edits atom coordinates.

## Pairwise contact/clash calculations

The `clashes` command computes Euclidean distances between atoms from two separately supplied files. It uses a chunked N-by-M calculation to reduce peak memory use. A pair is included when `distance <= cutoff`. The output is sorted by distance. If `--max-pairs` is specified, the scan returns after that number has been encountered; do not interpret a capped report as the globally nearest subset for large files.

This is a geometric screen, not an energetic clash score. No van der Waals radii, bonded exclusions, periodic boundaries, waters policy beyond the heavy-atom filter, or force-field-specific criteria are applied. The two coordinate files must already be in the same reference frame and use the same units.

## PDB I/O

The current release supports single-model, fixed-width PDB coordinate files. It updates coordinate columns for ATOM/HETATM records and preserves other records as text. Multi-model PDBs are rejected to avoid silently transforming/aligning an ensemble as if it were one structure. mmCIF, compressed files, crystallographic symmetry, unit-cell transforms, and trajectory formats are out of scope. Writes fail explicitly if a coordinate would overflow the fixed-width PDB coordinate fields; use mmCIF or a different format for such cases. Coordinate-associated `ANISOU`, `SIGATM`, and `SIGUIJ` records cause coordinate edits to fail rather than leaving their tensors stale. `CRYST1`, `SCALE`, `ORIGX`, and `MTRIX` records are preserved with an explicit warning; review or remove/update these records before using transformed structures as crystallographic models.

Alternate locations other than blank/A are not considered in fit matching or omitted from default pairwise contact selections. All coordinate records are transformed together. Element detection depends on the element column when present; malformed or non-standard PDB files should be inspected manually.

## Validation status

Automated tests use synthetic coordinates where expected transforms are known. They do not validate any particular MBD/chaperone structure, atom-selection string, protonation state, copper coordination geometry, or simulation force field. Before using an output in a scientific simulation:

1. inspect before/after structures visually;
2. confirm that the fitting atoms are biologically corresponding;
3. compare reported RMSD and atom counts with an independent method or VMD `measure fit`;
4. inspect inter-partner contacts/clashes with an appropriate system-specific threshold;
5. preserve the input files and transform JSON report alongside downstream simulation metadata.
