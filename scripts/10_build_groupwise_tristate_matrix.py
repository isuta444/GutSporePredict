#!/usr/bin/env python3
"""Build groupwise competitive three-state HMM gene calls."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class DomainHit:
    """One HMMER domain hit."""

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
    """Load unique target gene IDs."""

    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        if reader.fieldnames is None or "gene_id" not in reader.fieldnames:
            raise ValueError(
                f"Target TSV must contain gene_id: {path}"
            )

        gene_ids: list[str] = []
        seen: set[str] = set()

        for row in reader:
            gene_id = row["gene_id"].strip()

            if gene_id and gene_id not in seen:
                seen.add(gene_id)
                gene_ids.append(gene_id)

    if not gene_ids:
        raise ValueError(f"No target genes found: {path}")

    return gene_ids


def load_competition_groups(
    path: Path,
) -> dict[str, tuple[str, ...]]:
    """Load related HMM profiles that must compete one-to-one."""

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

            profile_length = int(fields[5])
            hmm_from = int(fields[15])
            hmm_to = int(fields[16])

            profile_coverage = (
                hmm_to - hmm_from + 1
            ) / profile_length

            hits.append(
                DomainHit(
                    genome_id=genome_id,
                    gene_id=gene_id,
                    target_id=fields[0],
                    full_evalue=float(fields[6]),
                    full_score=float(fields[7]),
                    domain_evalue=float(fields[12]),
                    domain_score=float(fields[13]),
                    profile_coverage=profile_coverage,
                )
            )

    return hits


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


def candidate_hits(
    hits: list[DomainHit],
    args: argparse.Namespace,
) -> list[DomainHit]:
    """Return hits that pass at least the uncertain threshold."""

    return [
        hit
        for hit in hits
        if classify_hit(hit, args) != "0"
    ]


def hit_key(hit: DomainHit) -> tuple[float, float, float, float]:
    """Return the ordering key for an individual hit."""

    return (
        hit.domain_score,
        hit.full_score,
        -hit.domain_evalue,
        hit.profile_coverage,
    )


def best_group_assignment(
    genes: tuple[str, ...],
    hits_by_gene: dict[str, list[DomainHit]],
) -> list[DomainHit]:
    """Find the maximum-score one-to-one group assignment."""

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
    """Build competitive three-state gene calls."""

    args = parse_args()
    gene_ids = load_gene_ids(args.targets)
    groups = load_competition_groups(args.competition_groups)

    grouped_gene_ids = {
        gene_id
        for group_gene_ids in groups.values()
        for gene_id in group_gene_ids
    }

    unknown_gene_ids = grouped_gene_ids.difference(gene_ids)

    if unknown_gene_ids:
        raise ValueError(
            "Competition groups contain unknown genes: "
            + ", ".join(sorted(unknown_gene_ids))
        )

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

    all_call_rows: list[dict[str, str]] = []
    all_detail_rows: list[dict[str, str]] = []

    for genome_dir in genome_dirs:
        genome_id = genome_dir.name
        hits_by_gene: dict[str, list[DomainHit]] = {}

        for gene_id in gene_ids:
            domtblout = genome_dir / f"{gene_id}.domtblout"

            if not domtblout.is_file():
                raise FileNotFoundError(domtblout)

            hits_by_gene[gene_id] = parse_domtblout(
                domtblout,
                genome_id,
                gene_id,
            )

        selected_by_gene: dict[str, DomainHit] = {}

        for gene_id in gene_ids:
            if gene_id in grouped_gene_ids:
                continue

            candidates = candidate_hits(
                hits_by_gene[gene_id],
                args,
            )

            if candidates:
                selected_by_gene[gene_id] = max(
                    candidates,
                    key=hit_key,
                )

        for group_gene_ids in groups.values():
            group_candidates = {
                gene_id: candidate_hits(
                    hits_by_gene[gene_id],
                    args,
                )
                for gene_id in group_gene_ids
            }

            for hit in best_group_assignment(
                group_gene_ids,
                group_candidates,
            ):
                selected_by_gene[hit.gene_id] = hit

        genome_calls: dict[str, str] = {
            "genome_id": genome_id,
        }

        for gene_id in gene_ids:
            hit = selected_by_gene.get(gene_id)
            call = classify_hit(hit, args)
            genome_calls[gene_id] = call

            all_detail_rows.append(
                {
                    "genome_id": genome_id,
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

        all_call_rows.append(genome_calls)

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
        writer.writerows(all_call_rows)

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
        writer.writerows(all_detail_rows)

    print(f"Wrote: {calls_path}")
    print(f"Wrote: {details_path}")
    print(
        f"Genomes: {len(genome_dirs)}, "
        f"genes: {len(gene_ids)}"
    )


if __name__ == "__main__":
    main()
