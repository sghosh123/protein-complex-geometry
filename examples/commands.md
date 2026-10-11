# Example commands

Assume `reference.pdb` contains the desired MBD reference orientation, while `chaperone.pdb` is the mobile partner.

## Align homologous backbone atoms

```bash
pcg align \
  --reference reference.pdb \
  --mobile chaperone.pdb \
  --reference-chain B \
  --mobile-chain A \
  --atom-name CA \
  --residues '41-46' \
  --output chaperone_aligned.pdb \
  --report chaperone_alignment.json
```

Use the biological mapping, not these example chain/residue values, for your actual structures. If the equivalent residues are numbered differently, the built-in number-based matcher will not infer that mapping.

## Place a partner using an explicit rigid transform

```bash
pcg transform \
  --input chaperone_aligned.pdb \
  --rotate-z 180 \
  --origin centroid \
  --translate '9,9,9' \
  --output chaperone_encounter_pose.pdb \
  --report placement.json
```

## Assess inter-partner clashes

```bash
pcg clashes \
  --structure-a mbd.pdb \
  --structure-b chaperone_encounter_pose.pdb \
  --cutoff 2.0 \
  --output clash_pairs.csv
```

Check visual orientation and the full contact report. A single cutoff does not replace chemistry-aware validation.


## Run the fully synthetic end-to-end example

```bash
bash examples/run_demo.sh
```

The included `examples/data/*.pdb` coordinates are synthetic and are only intended to smoke-test the package.


## Fit different residue-number ranges

For a homologous backbone segment, the reference and mobile residues may have different numbers or residue identities. Make the mapping explicit:

```bash
pcg align \
  --reference reference.pdb --mobile mobile.pdb \
  --reference-chain B --mobile-chain C --atom-name CA \
  --reference-residues '41-46' --mobile-residues '44-49' \
  --output mobile_aligned.pdb --report alignment.json
```

The selected atoms are paired by ascending residue order. This is safe only when the two ranges represent corresponding atoms; the program cannot infer homology for you.
