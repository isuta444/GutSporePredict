#!/usr/bin/env python3
"""Merge gene-call matrices sharing the same genomes."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        type=Path,
        action="append",
        required=True,
        help="Input gene_calls.tsv; specify more than once.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
    )
    return parser.parse_args()


def load_matrix(
    path: Path,
) -> tuple[list[str], dict[str, dict[str, str]]]:
    """Load one three-state gene-call matrix."""

    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        if reader.fieldnames is None or "genome_id" not in reader.fieldnames:
            raise ValueError(
                f"Matrix must contain genome_id: {path}"
            )

        gene_ids = [
            field
            for field in reader.fieldnames
            if field != "genome_id"
        ]

        rows: dict[str, dict[str, str]] = {}

        for row in reader:
            genome_id = row["genome_id"].strip()

            if not genome_id:
                raise ValueError(
                    f"Empty genome_id in: {path}"
                )

            if genome_id in rows:
                raise ValueError(
                    f"Duplicate genome_id in {path}: {genome_id}"
                )

            rows[genome_id] = {
                gene_id: row[gene_id].strip()
                for gene_id in gene_ids
            }

    return gene_ids, rows


def main() -> None:
    """Merge matrices and write a combined TSV."""

    args = parse_args()

    all_gene_ids: list[str] = []
    merged: dict[str, dict[str, str]] = {}
    expected_genomes: set[str] | None = None

    for path in args.input:
        gene_ids, rows = load_matrix(path)
        genomes = set(rows)

        if expected_genomes is None:
            expected_genomes = genomes
        elif genomes != expected_genomes:
            raise ValueError(
                f"Genome sets differ between matrices: {path}"
            )

        duplicated = set(all_gene_ids).intersection(gene_ids)

        if duplicated:
            raise ValueError(
                "Duplicate gene columns across matrices: "
                + ", ".join(sorted(duplicated))
            )

        all_gene_ids.extend(gene_ids)

        for genome_id, calls in rows.items():
            merged.setdefault(genome_id, {}).update(calls)

    args.output.parent.mkdir(parents=True, exist_ok=True)

    with args.output.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["genome_id", *all_gene_ids])

        for genome_id in sorted(merged):
            writer.writerow(
                [
                    genome_id,
                    *[
                        merged[genome_id][gene_id]
                        for gene_id in all_gene_ids
                    ],
                ]
            )

    print(f"Wrote: {args.output}")
    print(f"Genomes: {len(merged)}")
    print(f"Genes: {len(all_gene_ids)}")


if __name__ == "__main__":
    main()
