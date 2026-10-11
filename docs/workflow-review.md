# Workflow review: improvements over ad hoc coordinate scripts

| Previous workflow risk | Current handling |
|---|---|
| Absolute `/home/...` paths | Input and output paths are command-line arguments |
| Atom indices tied to a particular loaded PDB | Named atoms, chain filters, residue IDs, insertion codes, residue names |
| Hidden choice of fit atoms | `--atom-name`, chain arguments, and `--residues` are explicit |
| No automated check that fitting is well-defined | Requires at least three non-collinear pairs and reports RMSD |
| Rotation origin/order not self-documenting | Pivot and Euler order are documented and written to JSON |
| Manual translation offsets hard to audit later | Translation is in arguments and preserved in the report |
| No quantified inter-partner overlap screen | Chunked inter-structure atom-pair output with distance cutoff |

## Not solved automatically

This toolkit does not know whether a given encounter geometry is biologically plausible. Rigidly superposing domains can produce clashes, distortions at unresolved loops, or misaligned homologous residues when sequence/structure mapping is wrong. It also does not optimize an interface, relax a model, or validate copper coordination chemistry. Those require separate analysis and, where suitable, molecular mechanics or QM/MM calculations.
