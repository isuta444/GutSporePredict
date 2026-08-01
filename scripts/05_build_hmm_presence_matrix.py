#!/usr/bin/env python3
"""Build best-hit and presence/absence tables from HMMER tblout files."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class HMMHit:
    """One parsed HMMER tblout hit."""

    target_id: str
    query_id: str
    full_evalue: float
    full_score: float
    best_domain_evalue: float
    best_domain_score: float


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(
        description=(
            "Select the best HMMER hit for each genome/gene pair and "
            "build presence/absence tables."
        )
    )
    parser.add_argument(
        "--hmmsearch-dir",
        type=Path,
        required=True,
        help="Directory containing one subdirectory per genome.",
    )
    parser.add_argument(
        "--targets",
        type=Path,
        required=True,
        help="TSV file containing a gene_id column.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory for output TSV files.",
    )
    parser.add_argument(
        "--max-evalue",
        type=float,
        default=1e-5,
        help="Maximum full-sequence E-value for presence.",
    )
    parser.add_argument(
        "--max-domain-evalue",
        type=float,
        default=1e-4,
        help="Maximum best-domain E-value for presence.",
    )
    return parser.parse_args()


def load_target_ids(path: Path) -> list[str]:
    """Load target gene identifiers from a TSV file."""

    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        if reader.fieldnames is None or "gene_id" not in reader.fieldnames:
            raise ValueError(
                f"Target TSV must contain a gene_id column: {path}"
            )

        target_ids = [
            row["gene_id"].strip()
            for row in reader
            if row["gene_id"].strip()
        ]

    if not target_ids:
        raise ValueError(f"No target genes found in: {path}")

    return target_ids


def parse_tblout(path: Path) -> list[HMMHit]:
    """Parse HMMER --tblout output."""

    hits: list[HMMHit] = []

    with path.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()

            if not line or line.startswith("#"):
                continue

            fields = line.split()

            if len(fields) < 18:
                raise ValueError(
                    f"Invalid HMMER tblout row at {path}:{line_number}"
                )

            try:
                hits.append(
                    HMMHit(
                        target_id=fields[0],
                        query_id=fields[2],
                        full_evalue=float(fields[4]),
                        full_score=float(fields[5]),
                        best_domain_evalue=float(fields[7]),
                        best_domain_score=float(fields[8]),
                    )
                )
            except ValueError as exc:
                raise ValueError(
                    f"Non-numeric HMMER value at {path}:{line_number}"
                ) from exc

    return hits


def best_hit(hits: list[HMMHit]) -> HMMHit | None:
    """Return the strongest hit by score and E-value."""

    if not hits:
        return None

    return max(
        hits,
        key=lambda hit: (
            hit.full_score,
            hit.best_domain_score,
            -hit.full_evalue,
            -hit.best_domain_evalue,
        ),
    )


def main() -> None:
    """Build output tables."""

    args = parse_args()
    target_ids = load_target_ids(args.targets)

    genome_dirs = sorted(
        path
        for path in args.hmmsearch_dir.iterdir()
        if path.is_dir()
    )

    if not genome_dirs:
        raise ValueError(
            f"No genome directories found in: {args.hmmsearch_dir}"
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)

    best_hits_path = args.output_dir / "best_hits.tsv"
    presence_path = args.output_dir / "presence_absence.tsv"

    best_hit_rows: list[dict[str, str]] = []
    presence_by_genome: dict[str, dict[str, str]] = {}

    for genome_dir in genome_dirs:
        genome_id = genome_dir.name
        presence_by_genome[genome_id] = {}

        for gene_id in target_ids:
            tblout = genome_dir / f"{gene_id}.tblout"

            if not tblout.is_file():
                raise FileNotFoundError(
                    f"Missing HMMER tblout file: {tblout}"
                )

            hit = best_hit(parse_tblout(tblout))

            if hit is None:
                present = False
                row = {
                    "genome_id": genome_id,
                    "gene_id": gene_id,
                    "target_id": "",
                    "full_evalue": "",
                    "full_score": "",
                    "best_domain_evalue": "",
                    "best_domain_score": "",
                    "present": "0",
                }
            else:
                present = (
                    hit.full_evalue <= args.max_evalue
                    and hit.best_domain_evalue <= args.max_domain_evalue
                )
                row = {
                    "genome_id": genome_id,
                    "gene_id": gene_id,
                    "target_id": hit.target_id,
                    "full_evalue": f"{hit.full_evalue:.6g}",
                    "full_score": f"{hit.full_score:.3f}",
                    "best_domain_evalue": (
                        f"{hit.best_domain_evalue:.6g}"
                    ),
                    "best_domain_score": (
                        f"{hit.best_domain_score:.3f}"
                    ),
                    "present": "1" if present else "0",
                }

            best_hit_rows.append(row)
            presence_by_genome[genome_id][gene_id] = (
                "1" if present else "0"
            )

    best_hit_columns = [
        "genome_id",
        "gene_id",
        "target_id",
        "full_evalue",
        "full_score",
        "best_domain_evalue",
        "best_domain_score",
        "present",
    ]

    with best_hits_path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=best_hit_columns,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(best_hit_rows)

    with presence_path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["genome_id", *target_ids])

        for genome_id in sorted(presence_by_genome):
            writer.writerow(
                [
                    genome_id,
                    *[
                        presence_by_genome[genome_id][gene_id]
                        for gene_id in target_ids
                    ],
                ]
            )

    print(f"Wrote: {best_hits_path}")
    print(f"Wrote: {presence_path}")
    print(
        f"Genomes: {len(genome_dirs)}, "
        f"targets: {len(target_ids)}, "
        f"comparisons: {len(best_hit_rows)}"
    )


if __name__ == "__main__":
    main()
