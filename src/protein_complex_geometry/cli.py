"""Command-line interface for protein-complex-geometry."""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import sys
import numpy as np
from .alignment import align_structures
from .clashes import interstructure_pairs
from .geometry import transform_from_euler
from .pdbio import PDBStructure


def _triple(text: str, name: str) -> tuple[float, float, float]:
    try:
        values = tuple(float(v.strip()) for v in text.split(","))
    except ValueError as exc:
        raise ValueError(f"{name} must be three comma-separated numbers") from exc
    if len(values) != 3 or not np.isfinite(values).all():
        raise ValueError(f"{name} must be three finite comma-separated numbers")
    return values  # type: ignore[return-value]


def parse_residues(text: str | None) -> set[int] | None:
    if not text:
        return None
    result: set[int] = set()
    try:
        for item in text.split(","):
            item = item.strip()
            if not item:
                raise ValueError
            if "-" in item:
                lo_s, hi_s = item.split("-", 1)
                lo, hi = int(lo_s), int(hi_s)
                if hi < lo:
                    raise ValueError
                result.update(range(lo, hi + 1))
            else:
                result.add(int(item))
    except ValueError as exc:
        raise argparse.ArgumentTypeError("residues must look like '30-33,41-46,52'") from exc
    return result


def _parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="pcg", description="Rigid-body protein structure geometry tools")
    sub = p.add_subparsers(dest="command", required=True)
    a = sub.add_parser("align", help="fit a mobile PDB onto a reference PDB")
    a.add_argument("--reference", required=True); a.add_argument("--mobile", required=True)
    a.add_argument("--output", required=True); a.add_argument("--report")
    a.add_argument("--reference-chain"); a.add_argument("--mobile-chain")
    a.add_argument("--atom-name", default="CA")
    a.add_argument("--residues", type=parse_residues,
                   help="same residue-number selection on both structures")
    a.add_argument("--reference-residues", type=parse_residues,
                   help="reference fit residues; pair by ascending order with --mobile-residues")
    a.add_argument("--mobile-residues", type=parse_residues,
                   help="mobile fit residues; pair by ascending order with --reference-residues")
    t = sub.add_parser("transform", help="rotate and translate all coordinates in a PDB")
    t.add_argument("--input", required=True); t.add_argument("--output", required=True)
    t.add_argument("--rotate-x", type=float, default=0.0); t.add_argument("--rotate-y", type=float, default=0.0)
    t.add_argument("--rotate-z", type=float, default=0.0)
    t.add_argument("--origin", default="0,0,0", help="rotation origin: x,y,z or 'centroid'")
    t.add_argument("--translate", default="0,0,0", help="translation in Angstroms: x,y,z")
    t.add_argument("--report")
    c = sub.add_parser("clashes", help="find inter-structure heavy-atom pairs at/below a cutoff")
    c.add_argument("--structure-a", required=True); c.add_argument("--structure-b", required=True)
    c.add_argument("--cutoff", type=float, default=2.0, help="distance cutoff in Angstroms (default: 2.0)")
    c.add_argument("--chain-a"); c.add_argument("--chain-b")
    c.add_argument("--include-hydrogen", action="store_true")
    c.add_argument("--chunk-size", type=int, default=512); c.add_argument("--max-pairs", type=int)
    c.add_argument("--output", required=True, help="CSV output path")
    return p


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        if args.command == "align":
            ref = PDBStructure.read(args.reference); mob = PDBStructure.read(args.mobile)
            metadata_warnings = mob.coordinate_frame_warnings()
            for warning in metadata_warnings:
                print(f"pcg: warning: {warning}", file=sys.stderr)
            aligned, result = align_structures(ref, mob,
                reference_chain=args.reference_chain, mobile_chain=args.mobile_chain,
                residue_ids=args.residues, reference_residue_ids=args.reference_residues,
                mobile_residue_ids=args.mobile_residues, atom_name=args.atom_name)
            aligned.write(args.output)
            report = result.as_dict()
            report.update({"reference": str(Path(args.reference)), "mobile": str(Path(args.mobile)),
                           "output": str(Path(args.output)), "atom_name": args.atom_name,
                           "reference_chain": args.reference_chain, "mobile_chain": args.mobile_chain,
                           "residue_filter_same_number": sorted(args.residues) if args.residues is not None else None,
                           "reference_residue_filter": sorted(args.reference_residues) if args.reference_residues is not None else None,
                           "mobile_residue_filter": sorted(args.mobile_residues) if args.mobile_residues is not None else None,
                           "coordinate_frame_warnings": metadata_warnings})
            if args.report:
                Path(args.report).parent.mkdir(parents=True, exist_ok=True)
                Path(args.report).write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
            print(f"Matched atoms: {result.n_atoms}\nRMSD before: {result.rmsd_before_angstrom:.4f} Å"
                  f"\nRMSD after:  {result.rmsd_after_angstrom:.4f} Å\nWrote: {args.output}")
            return 0
        if args.command == "transform":
            structure = PDBStructure.read(args.input)
            metadata_warnings = structure.coordinate_frame_warnings()
            for warning in metadata_warnings:
                print(f"pcg: warning: {warning}", file=sys.stderr)
            translation = _triple(args.translate, "translation")
            if args.origin.strip().lower() == "centroid":
                origin = tuple(structure.coordinates.mean(axis=0))
            else:
                origin = _triple(args.origin, "origin")
            transform = transform_from_euler(args.rotate_x, args.rotate_y, args.rotate_z,
                                             origin=origin, translation=translation)
            structure.transformed(transform).write(args.output)
            report = {"input": args.input, "output": args.output,
                      "rotation_degrees_xyz": [args.rotate_x, args.rotate_y, args.rotate_z],
                      "rotation_origin_angstrom": list(origin), "translation_angstrom": list(translation),
                      "transform": transform.as_dict(), "coordinate_frame_warnings": metadata_warnings}
            if args.report:
                Path(args.report).parent.mkdir(parents=True, exist_ok=True)
                Path(args.report).write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
            print(f"Transformed {len(structure.atoms)} coordinate records. Wrote: {args.output}")
            return 0
        if args.command == "clashes":
            a = PDBStructure.read(args.structure_a); b = PDBStructure.read(args.structure_b)
            pairs = interstructure_pairs(a, b, cutoff_angstrom=args.cutoff,
                chain_a=args.chain_a, chain_b=args.chain_b, heavy_only=not args.include_hydrogen,
                chunk_size=args.chunk_size, max_pairs=args.max_pairs)
            output = Path(args.output); output.parent.mkdir(parents=True, exist_ok=True)
            fieldnames = list(pairs[0].as_dict()) if pairs else ["atom_a", "atom_b", "distance_angstrom", "chain_a", "residue_a", "chain_b", "residue_b"]
            with output.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=fieldnames); writer.writeheader()
                writer.writerows(pair.as_dict() for pair in pairs)
            print(f"Pairs within {args.cutoff:g} Å: {len(pairs)}\nWrote: {output}")
            if args.max_pairs is not None and len(pairs) >= args.max_pairs:
                print("pcg: warning: output cap reached; additional qualifying pairs may be omitted", file=sys.stderr)
            if pairs: print(f"Minimum inter-structure distance in reported pairs: {pairs[0].distance_angstrom:.3f} Å")
            return 0
    except (OSError, ValueError, np.linalg.LinAlgError) as exc:
        print(f"pcg: error: {exc}", file=sys.stderr)
        return 2
    return 2
