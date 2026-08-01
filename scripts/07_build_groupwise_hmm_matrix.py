#!/usr/bin/env python3
"""Build HMM presence matrix with groupwise one-to-one assignment."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class DomainHit:
    """One passing HMMER domain hit."""

    genome_id: str
    gene_id: str
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
    parser.add_argument("--competition-groups", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--max-evalue", type=float, default=1e-5)
    parser.add_argument("--max-domain-evalue", type=float, default=1e-4)
    parser.add_argument(
        "--minimum-profile-coverage",
        type=float,
        default=0.45,
    )
    return parser.parse_args()


def load_target_ids(path: Path) -> list[str]:
    """Load target gene IDs."""

    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        if reader.fieldnames is None or "gene_id" not in reader.fieldnames:
            raise ValueError("Target TSV must contain gene_id.")

        targets = [
            row["gene_id"].strip()
            for row in reader
            if row["gene_id"].strip()
        ]

    if not targets:
        raise ValueError(f"No targets found in: {path}")

    return targets


def load_competition_groups(path: Path) -> dict[str, tuple[str, ...]]:
    """Load profile competition groups."""

    groups: dict[str, list[str]] = defaultdict(list)

    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        required = {"group_id", "gene_id"}

        if reader.fieldnames is None or not required.issubset(
            reader.fieldnames
        ):
            raise ValueError(
                "Competition TSV must contain group_id and gene_id."
            )

        for row in reader:
            group_id = row["group_id"].strip()
            gene_id = row["gene_id"].strip()

            if group_id and gene_id:
                groups[group_id].append(gene_id)

    return {
        group_id: tuple(gene_ids)
        for group_id, gene_ids in groups.items()
    }


def parse_domtblout(
    path: Path,
    genome_id: str,
    gene_id: str,
    *,
    max_evalue: float,
    max_domain_evalue: float,
    minimum_profile_coverage: float,
) -> list[DomainHit]:
    """Parse and filter HMMER domtblout hits."""

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

            profile_length = int(fields[5])
            hmm_from = int(fields[15])
            hmm_to = int(fields[16])

            profile_coverage = (
                hmm_to - hmm_from + 1
            ) / profile_length

            hit = DomainHit(
                genome_id=genome_id,
                gene_id=gene_id,
                target_id=fields[0],
                full_evalue=float(fields[6]),
                full_score=float(fields[7]),
                domain_evalue=float(fields[12]),
                domain_score=float(fields[13]),
                profile_coverage=profile_coverage,
            )

            if (
                hit.full_evalue <= max_evalue
                and hit.domain_evalue <= max_domain_evalue
                and hit.profile_coverage >= minimum_profile_coverage
            ):
                hits.append(hit)

    return hits


def hit_key(hit: DomainHit) -> tuple[float, float, float]:
    """Return the ordering key for individual hits."""

    return (
        hit.domain_score,
        hit.full_score,
        hit.profile_coverage,
    )


def best_group_assignment(
    genes: tuple[str, ...],
    hits_by_gene: dict[str, list[DomainHit]],
) -> list[DomainHit]:
    """Find the maximum-score one-to-one profile assignment."""

    best_hits: list[DomainHit] = []
    best_score = float("-inf")

    def search(
        index: int,
        used_targets: set[str],
        selected: list[DomainHit],
        score: float,
    ) -> None:
        nonlocal best_hits, best_score

        if index == len(genes):
            if score > best_score:
                best_score = score
                best_hits = list(selected)
            return

        gene_id = genes[index]

        # A profile may remain unassigned.
        search(
            index + 1,
            used_targets,
            selected,
            score,
        )

        for hit in hits_by_gene.get(gene_id, []):
            if hit.target_id in used_targets:
                continue

            used_targets.add(hit.target_id)
            selected.append(hit)

            search(
                index + 1,
                used_targets,
                selected,
                score + hit.domain_score,
            )

            selected.pop()
            used_targets.remove(hit.target_id)

    search(0, set(), [], 0.0)
    return best_hits


def main() -> None:
    """Build assignments and presence/absence matrix."""

    args = parse_args()
    target_ids = load_target_ids(args.targets)
    groups = load_competition_groups(args.competition_groups)

    grouped_gene_ids = {
        gene_id
        for genes in groups.values()
        for gene_id in genes
    }

    unknown_grouped = grouped_gene_ids.difference(target_ids)

    if unknown_grouped:
        raise ValueError(
            "Competition groups contain unknown targets: "
            + ", ".join(sorted(unknown_grouped))
        )

    genome_dirs = sorted(
        path
        for path in args.hmmsearch_dir.iterdir()
        if path.is_dir()
    )

    assignments: list[DomainHit] = []

    for genome_dir in genome_dirs:
        genome_id = genome_dir.name
        hits_by_gene: dict[str, list[DomainHit]] = {}

        for gene_id in target_ids:
            path = genome_dir / f"{gene_id}.domtblout"

            if not path.is_file():
                raise FileNotFoundError(path)

            hits_by_gene[gene_id] = parse_domtblout(
                path,
                genome_id,
                gene_id,
                max_evalue=args.max_evalue,
                max_domain_evalue=args.max_domain_evalue,
                minimum_profile_coverage=(
                    args.minimum_profile_coverage
                ),
            )

        # Non-competing profiles: select the best hit independently.
        for gene_id in target_ids:
            if gene_id in grouped_gene_ids:
                continue

            hits = hits_by_gene[gene_id]

            if hits:
                assignments.append(max(hits, key=hit_key))

        # Related profiles: optimize all assignments together.
        for genes in groups.values():
            assignments.extend(
                best_group_assignment(
                    genes,
                    hits_by_gene,
                )
            )

    assignments.sort(
        key=lambda hit: (
            hit.genome_id,
            target_ids.index(hit.gene_id),
        )
    )

    args.output_dir.mkdir(parents=True, exist_ok=True)

    assignments_path = args.output_dir / "assignments.tsv"
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
            genome_id = genome_dir.name

            writer.writerow(
                [
                    genome_id,
                    *[
                        "1"
                        if (genome_id, gene_id) in assigned_pairs
                        else "0"
                        for gene_id in target_ids
                    ],
                ]
            )

    print(f"Wrote: {assignments_path}")
    print(f"Wrote: {matrix_path}")
    print(f"Assignments: {len(assignments)}")


if __name__ == "__main__":
    main()
