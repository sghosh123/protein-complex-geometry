#!/usr/bin/env bash
set -euo pipefail

mkdir -p demo_output
pcg align \
  --reference examples/data/reference.pdb \
  --mobile examples/data/mobile.pdb \
  --reference-chain B --mobile-chain C \
  --atom-name CA --reference-residues '30-33' --mobile-residues '44-47' \
  --output demo_output/mobile_aligned.pdb \
  --report demo_output/alignment.json

pcg transform \
  --input demo_output/mobile_aligned.pdb \
  --rotate-z 180 --origin centroid --translate '0,0,0' \
  --output demo_output/mobile_rotated.pdb \
  --report demo_output/transform.json

pcg clashes \
  --structure-a examples/data/partner_a.pdb \
  --structure-b examples/data/partner_b.pdb \
  --cutoff 2.0 --output demo_output/close_pairs.csv
