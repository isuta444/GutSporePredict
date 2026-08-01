#!/usr/bin/env python3
"""Build a three-state HMMER gene-call matrix."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class DomainHit:
    """One parsed HMMER domain hit."""

    target_id: str
    full_evalue: float
    full_score: float
    domain_evalue: float
    domain_score: float
    profile_coverage: float


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser()
    parser.add_argument("--hmmsearch-dir", type=Path, required=True)
    parser.add_argument("--targets", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)

    parser.add_argument(
        "--present-max-evalue",
        type=float,
        default=1e-5,
    )
    parser.add_argument(
        "--present-max-domain-evalue",
        type=float,
        default=1e-4,
    )
    parser.add_argument(
        "--present-min-coverage",
        type=float,
        default=0.45,
    )

    parser.add_argument(
        "--uncertain-max-evalue",
        type=float,
        default=1e-3,
    )
    parser.add_argument(
        "--uncertain-max-domain-evalue",
        type=float,
        default=1e-3,
    )
    parser.add_argument(
        "--uncertain-min-coverage",
        type=float,
        default=0.25,
    )

    return parser.parse_args()


def load_gene_ids(path: Path) -> list[str]:
    """Load target gene IDs."""

    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        if reader.fieldnames is None or "gene_id" not in reader.fieldnames:
            raise ValueError(
                f"Target TSV must contain gene_id: {path}"
            )

        gene_ids = []
        seen = set()

        for row in reader:
            gene_id = row["gene_id"].strip()

            if gene_id and gene_id not in seen:
                seen.add(gene_id)
                gene_ids.append(gene_id)

    if not gene_ids:
        raise ValueError(f"No target genes found: {path}")

    return gene_ids


def parse_domtblout(path: Path) -> list[DomainHit]:
    """Parse one HMMER domtblout file."""

    hits = []

    with path.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()

            if not line or line.startswith("#"):
                continue

            fields = line.split(maxsplit=22)

            if len(fields) < 22:
                raise ValueError(
                    f"Invalid domtblout row: {path}:{line_number}"
                )

            profile_length = int(fields[5])
            hmm_from = int(fields[15])
            hmm_to = int(fields[16])

            coverage = (
                hmm_to - hmm_from + 1
            ) / profile_length

            hits.append(
                DomainHit(
                    target_id=fields[0],
                    full_evalue=float(fields[6]),
                    full_score=float(fields[7]),
                    domain_evalue=float(fields[12]),
                    domain_score=float(fields[13]),
                    profile_coverage=coverage,
                )
            )

    return hits


def best_hit(hits: list[DomainHit]) -> DomainHit | None:
    """Return the strongest domain hit."""

    if not hits:
        return None

    return max(
        hits,
        key=lambda hit: (
            hit.domain_score,
            hit.full_score,
            -hit.domain_evalue,
            hit.profile_coverage,
        ),
    )


def classify_hit(
    hit: DomainHit | None,
    args: argparse.Namespace,
) -> str:
    """Classify a hit as present, uncertain or absent."""

    if hit is None:
        return "0"

    if (
        hit.full_evalue <= args.present_max_evalue
        and hit.domain_evalue <= args.present_max_domain_evalue
        and hit.profile_coverage >= args.present_min_coverage
    ):
        return "1"

    if (
        hit.full_evalue <= args.uncertain_max_evalue
        and hit.domain_evalue <= args.uncertain_max_domain_evalue
        and hit.profile_coverage >= args.uncertain_min_coverage
    ):
        return "?"

    return "0"


def main() -> None:
    """Build three-state gene calls."""

    args = parse_args()
    gene_ids = load_gene_ids(args.targets)

    genome_dirs = sorted(
        path
        for path in args.hmmsearch_dir.iterdir()
        if path.is_dir()
    )

    if not genome_dirs:
        raise ValueError(
            f"No genome directories found: {args.hmmsearch_dir}"
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)

    calls_path = args.output_dir / "gene_calls.tsv"
    details_path = args.output_dir / "gene_call_details.tsv"

    call_rows = []
    detail_rows = []

    for genome_dir in genome_dirs:
        genome_calls = {"genome_id": genome_dir.name}

        for gene_id in gene_ids:
            path = genome_dir / f"{gene_id}.domtblout"

            if not path.is_file():
                raise FileNotFoundError(path)

            hit = best_hit(parse_domtblout(path))
            call = classify_hit(hit, args)
            genome_calls[gene_id] = call

            detail_rows.append(
                {
                    "genome_id": genome_dir.name,
                    "gene_id": gene_id,
                    "call": call,
                    "target_id": "" if hit is None else hit.target_id,
                    "full_evalue": (
                        ""
                        if hit is None
                        else f"{hit.full_evalue:.6g}"
                    ),
                    "full_score": (
                        ""
                        if hit is None
                        else f"{hit.full_score:.3f}"
                    ),
                    "domain_evalue": (
                        ""
                        if hit is None
                        else f"{hit.domain_evalue:.6g}"
                    ),
                    "domain_score": (
                        ""
                        if hit is None
                        else f"{hit.domain_score:.3f}"
                    ),
                    "profile_coverage": (
                        ""
                        if hit is None
                        else f"{hit.profile_coverage:.3f}"
                    ),
                }
            )

        call_rows.append(genome_calls)

    with calls_path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["genome_id", *gene_ids],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(call_rows)

    detail_columns = [
        "genome_id",
        "gene_id",
        "call",
        "target_id",
        "full_evalue",
        "full_score",
        "domain_evalue",
        "domain_score",
        "profile_coverage",
    ]

    with details_path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=detail_columns,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(detail_rows)

    print(f"Wrote: {calls_path}")
    print(f"Wrote: {details_path}")
    print(
        f"Genomes: {len(genome_dirs)}, "
        f"genes: {len(gene_ids)}"
    )


if __name__ == "__main__":
    main()
