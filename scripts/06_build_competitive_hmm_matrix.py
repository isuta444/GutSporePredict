#!/usr/bin/env python3
"""Build a competitive HMMER presence/absence matrix."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class DomainHit:
    """One HMMER domain hit."""

    genome_id: str
    gene_id: str
    target_id: str
    target_length: int
    profile_length: int
    full_evalue: float
    full_score: float
    domain_evalue: float
    domain_score: float
    hmm_from: int
    hmm_to: int

    @property
    def profile_coverage(self) -> float:
        """Return the fraction of the HMM profile covered."""

        aligned = self.hmm_to - self.hmm_from + 1
        return aligned / self.profile_length


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser()
    parser.add_argument("--hmmsearch-dir", type=Path, required=True)
    parser.add_argument("--targets", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--max-evalue", type=float, default=1e-5)
    parser.add_argument("--max-domain-evalue", type=float, default=1e-4)
    parser.add_argument("--minimum-profile-coverage", type=float, default=0.45)
    return parser.parse_args()


def load_targets(path: Path) -> list[str]:
    """Load gene IDs from the target TSV."""

    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        if reader.fieldnames is None or "gene_id" not in reader.fieldnames:
            raise ValueError("Target TSV must contain gene_id.")

        return [
            row["gene_id"].strip()
            for row in reader
            if row["gene_id"].strip()
        ]


def parse_domtblout(
    path: Path,
    genome_id: str,
    gene_id: str,
) -> list[DomainHit]:
    """Parse one HMMER domtblout file."""

    hits: list[DomainHit] = []

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

            hits.append(
                DomainHit(
                    genome_id=genome_id,
                    gene_id=gene_id,
                    target_id=fields[0],
                    target_length=int(fields[2]),
                    profile_length=int(fields[5]),
                    full_evalue=float(fields[6]),
                    full_score=float(fields[7]),
                    domain_evalue=float(fields[12]),
                    domain_score=float(fields[13]),
                    hmm_from=int(fields[15]),
                    hmm_to=int(fields[16]),
                )
            )

    return hits


def main() -> None:
    """Build competitive assignments and presence matrix."""

    args = parse_args()
    target_ids = load_targets(args.targets)

    passing_hits: list[DomainHit] = []

    genome_dirs = sorted(
        path for path in args.hmmsearch_dir.iterdir() if path.is_dir()
    )

    for genome_dir in genome_dirs:
        for gene_id in target_ids:
            path = genome_dir / f"{gene_id}.domtblout"

            if not path.is_file():
                raise FileNotFoundError(path)

            for hit in parse_domtblout(path, genome_dir.name, gene_id):
                if (
                    hit.full_evalue <= args.max_evalue
                    and hit.domain_evalue <= args.max_domain_evalue
                    and hit.profile_coverage
                    >= args.minimum_profile_coverage
                ):
                    passing_hits.append(hit)

    # One protein may match several related HMMs.
    # Assign it only to the profile with the highest domain score.
    best_by_protein: dict[tuple[str, str], DomainHit] = {}

    for hit in passing_hits:
        key = (hit.genome_id, hit.target_id)
        current = best_by_protein.get(key)

        if current is None or (
            hit.domain_score,
            -hit.domain_evalue,
            hit.profile_coverage,
        ) > (
            current.domain_score,
            -current.domain_evalue,
            current.profile_coverage,
        ):
            best_by_protein[key] = hit

    assignments = sorted(
        best_by_protein.values(),
        key=lambda hit: (
            hit.genome_id,
            hit.gene_id,
            -hit.domain_score,
        ),
    )

    args.output_dir.mkdir(parents=True, exist_ok=True)

    assignments_path = args.output_dir / "competitive_assignments.tsv"
    matrix_path = args.output_dir / "presence_absence.tsv"

    with assignments_path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "genome_id",
                "gene_id",
                "target_id",
                "full_evalue",
                "full_score",
                "domain_evalue",
                "domain_score",
                "profile_coverage",
            ]
        )

        for hit in assignments:
            writer.writerow(
                [
                    hit.genome_id,
                    hit.gene_id,
                    hit.target_id,
                    f"{hit.full_evalue:.6g}",
                    f"{hit.full_score:.3f}",
                    f"{hit.domain_evalue:.6g}",
                    f"{hit.domain_score:.3f}",
                    f"{hit.profile_coverage:.3f}",
                ]
            )

    assigned_pairs = {
        (hit.genome_id, hit.gene_id)
        for hit in assignments
    }

    with matrix_path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["genome_id", *target_ids])

        for genome_dir in genome_dirs:
            writer.writerow(
                [
                    genome_dir.name,
                    *[
                        "1"
                        if (genome_dir.name, gene_id) in assigned_pairs
                        else "0"
                        for gene_id in target_ids
                    ],
                ]
            )

    print(f"Wrote: {assignments_path}")
    print(f"Wrote: {matrix_path}")
    print(f"Passing domain hits: {len(passing_hits)}")
    print(f"Competitive assignments: {len(assignments)}")


if __name__ == "__main__":
    main()
