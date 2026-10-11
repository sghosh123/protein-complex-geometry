# Source-provenance notes

This repository was designed after reviewing user-provided structure-preparation scripts in the copper-trafficking workspace, including these workflow patterns:

- `align_yeast_mbd.tcl`: Cα-based fit between apo/holo MBD structures using a residue selection and VMD `measure fit`.
- `align_yeast_chap_ref_mbd_holo_final.tcl`: fit chaperone structures against a reference MBD using selected backbone Cα atoms.
- `change_holo_chap_rotate.tcl`: rigid rotation of a partner around a Cartesian axis.
- `change_holo_chap_translate.tcl`: translate a partner after using reference anchor atoms and an additional offset.
- `change_apo_mbd_rotate.tcl` and `change_apo_mbd_translate.tcl`: sequential rotations and fixed coordinate translations of the MBD.
- `show_y_apo_mbd_holo_chap_geom.tcl`: visualize multiple placed encounter geometries.

These scripts were treated as scientific use cases, not as production-ready modules. The Python code in this repository is a clean generalization and does not copy those Tcl files verbatim.

## Corrective design choices

- Removed workstation-specific absolute paths.
- Replaced brittle atom indices with named atom/residue/chain matching.
- Made rotation order, pivot, angle units, and translation units explicit.
- Added guards for missing correspondence, too few points, collinear fits, and improper reflection transforms.
- Added fit and transform JSON reports.
- Added inter-structure geometric close-contact reporting.
- Added synthetic tests for numerical invariants and PDB text preservation.

The biological residue selection and mapping must still be supplied and checked by the researcher for each system. No assumption is made that a single residue selection is valid across all MBD/chaperone homologs.
