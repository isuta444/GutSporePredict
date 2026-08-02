#!/usr/bin/env python3
"""Merge lineage-aware Clostridia SpoIIQ evidence into gene calls."""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class DomainHit:
    """One HMMER domain hit."""

    target_id: str
    target_length: int
    full_evalue: float
    full_score: float
    domain_evalue: float
    domain_score: float
    profile_coverage: float


@dataclass(frozen=True)
class EvaluatedCandidate:
    """Clostridia-type SpoIIQ candidate with context evidence."""

    hit: DomainHit
    spoiid_target: str
    spoiid_distance: int | None
    sequence_strength: str
    context_strength: str
    evidence_class: str


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(
        description=(
            "Merge Clostridia-specific SpoIIQ evidence into an "
            "existing sporulation gene-call matrix."
        )
    )

    parser.add_argument(
        "--input-calls",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "--input-details",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "--clostridia-search-dir",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "--sporulation-search-dir",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "--output-calls",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "--output-details",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "--audit-output",
        type=Path,
        required=True,
    )

    return parser.parse_args()


def split_protein_id(
    protein_id: str,
) -> tuple[str, int] | None:
    """Split a Prodigal-style protein ID."""

    match = re.match(r"(.+)_([0-9]+)$", protein_id)

    if match is None:
        return None

    return match.group(1), int(match.group(2))


def parse_domtblout(path: Path) -> list[DomainHit]:
    """Parse one HMMER domtblout file."""

    if not path.is_file():
        return []

    hits: list[DomainHit] = []

    with path.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()

            if not line or line.startswith("#"):
                continue

            fields = line.split(maxsplit=22)

            if len(fields) < 22:
                raise ValueError(f"Invalid domtblout row: {path}:{line_number}")

            profile_length = int(fields[5])
            hmm_from = int(fields[15])
            hmm_to = int(fields[16])

            hits.append(
                DomainHit(
                    target_id=fields[0],
                    target_length=int(fields[2]),
                    full_evalue=float(fields[6]),
                    full_score=float(fields[7]),
                    domain_evalue=float(fields[12]),
                    domain_score=float(fields[13]),
                    profile_coverage=(hmm_to - hmm_from + 1) / profile_length,
                )
            )

    return hits


def deduplicate_hits(
    hits: list[DomainHit],
) -> list[DomainHit]:
    """Keep the strongest domain row for each target protein."""

    best: dict[str, DomainHit] = {}

    for hit in hits:
        current = best.get(hit.target_id)

        if current is None or hit_key(hit) > hit_key(current):
            best[hit.target_id] = hit

    return list(best.values())


def hit_key(hit: DomainHit) -> tuple[float, float, float, float]:
    """Rank individual HMM hits."""

    return (
        hit.domain_score,
        hit.full_score,
        -hit.domain_evalue,
        hit.profile_coverage,
    )


def nearest_spoiid(
    candidate_id: str,
    spoiid_hits: list[DomainHit],
) -> tuple[str, int | None]:
    """Find the nearest SpoIID candidate on the same contig."""

    candidate = split_protein_id(candidate_id)

    if candidate is None:
        return "", None

    candidate_contig, candidate_number = candidate
    neighbors: list[tuple[DomainHit, int]] = []

    for hit in spoiid_hits:
        parsed = split_protein_id(hit.target_id)

        if parsed is None:
            continue

        contig, gene_number = parsed

        if contig != candidate_contig:
            continue

        neighbors.append(
            (
                hit,
                gene_number - candidate_number,
            )
        )

    if not neighbors:
        return "", None

    best_hit, distance = min(
        neighbors,
        key=lambda item: (
            abs(item[1]),
            -item[0].domain_score,
        ),
    )

    return best_hit.target_id, distance


def classify_sequence(hit: DomainHit) -> str:
    """Classify sequence evidence."""

    length_ok = 180 <= hit.target_length <= 330

    if (
        hit.full_evalue <= 1e-20
        and hit.domain_evalue <= 1e-20
        and hit.profile_coverage >= 0.80
        and length_ok
    ):
        return "strong"

    if (
        hit.full_evalue <= 1e-10
        and hit.domain_evalue <= 1e-10
        and hit.profile_coverage >= 0.60
        and 150 <= hit.target_length <= 400
    ):
        return "moderate"

    if (
        hit.full_evalue <= 1e-5
        and hit.domain_evalue <= 1e-4
        and hit.profile_coverage >= 0.40
    ):
        return "weak"

    return "insufficient"


def classify_context(distance: int | None) -> str:
    """Classify SpoIID neighborhood evidence."""

    if distance is None:
        return "none"

    absolute_distance = abs(distance)

    if absolute_distance <= 2:
        return "adjacent"

    if absolute_distance <= 10:
        return "local"

    return "distant"


def combine_evidence(
    sequence_strength: str,
    context_strength: str,
) -> str:
    """Combine sequence and genomic-context evidence."""

    if sequence_strength == "strong" and context_strength == "adjacent":
        return "high"

    if sequence_strength == "strong" and context_strength in {"local", "none"}:
        return "moderate"

    if sequence_strength == "moderate" and context_strength == "adjacent":
        return "moderate"

    if sequence_strength in {"strong", "moderate"}:
        return "provisional"

    if sequence_strength == "weak" and context_strength == "adjacent":
        return "provisional"

    return "uncertain"


def evidence_rank(candidate: EvaluatedCandidate) -> tuple[float, ...]:
    """Rank evaluated candidates."""

    class_rank = {
        "uncertain": 0,
        "provisional": 1,
        "moderate": 2,
        "high": 3,
    }

    context_rank = {
        "none": 0,
        "distant": 0,
        "local": 1,
        "adjacent": 2,
    }

    sequence_rank = {
        "insufficient": 0,
        "weak": 1,
        "moderate": 2,
        "strong": 3,
    }

    return (
        float(class_rank[candidate.evidence_class]),
        float(context_rank[candidate.context_strength]),
        float(sequence_rank[candidate.sequence_strength]),
        candidate.hit.domain_score,
        candidate.hit.full_score,
        candidate.hit.profile_coverage,
    )


def evaluate_genome(
    genome_id: str,
    clostridia_dir: Path,
    sporulation_dir: Path,
) -> EvaluatedCandidate | None:
    """Evaluate Clostridia-type SpoIIQ candidates in one genome."""

    clostridia_hits = deduplicate_hits(
        parse_domtblout(clostridia_dir / genome_id / "spoIIQ_Clostridia.domtblout")
    )

    if not clostridia_hits:
        return None

    spoiid_hits = deduplicate_hits(
        parse_domtblout(sporulation_dir / genome_id / "spoIID.domtblout")
    )

    evaluated: list[EvaluatedCandidate] = []

    for hit in clostridia_hits:
        spoiid_target, distance = nearest_spoiid(
            hit.target_id,
            spoiid_hits,
        )

        sequence_strength = classify_sequence(hit)
        context_strength = classify_context(distance)

        evaluated.append(
            EvaluatedCandidate(
                hit=hit,
                spoiid_target=spoiid_target,
                spoiid_distance=distance,
                sequence_strength=sequence_strength,
                context_strength=context_strength,
                evidence_class=combine_evidence(
                    sequence_strength,
                    context_strength,
                ),
            )
        )

    return max(evaluated, key=evidence_rank)


def merged_call(
    original_call: str,
    candidate: EvaluatedCandidate | None,
) -> str:
    """Merge the original and Clostridia-specific calls."""

    if original_call == "1":
        return "1"

    if candidate is None:
        return original_call

    if candidate.evidence_class == "high":
        return "1"

    if candidate.evidence_class == "moderate":
        return "?"

    # Provisional evidence is retained in the audit table but does not
    # promote an otherwise absent call.
    return original_call


def main() -> None:
    """Merge SpoIIQ lineage evidence."""

    args = parse_args()

    with args.input_calls.open(
        encoding="utf-8",
        newline="",
    ) as handle:
        call_reader = csv.DictReader(handle, delimiter="\t")

        if call_reader.fieldnames is None:
            raise ValueError("Input call matrix has no header.")

        call_columns = list(call_reader.fieldnames)
        call_rows = list(call_reader)

    if "spoIIQ" not in call_columns:
        raise ValueError("Input call matrix has no spoIIQ column.")

    original_by_genome = {row["genome_id"]: row["spoIIQ"] for row in call_rows}
    candidates: dict[str, EvaluatedCandidate | None] = {}

    for row in call_rows:
        genome_id = row["genome_id"]
        candidates[genome_id] = evaluate_genome(
            genome_id,
            args.clostridia_search_dir,
            args.sporulation_search_dir,
        )

        row["spoIIQ"] = merged_call(
            row["spoIIQ"],
            candidates[genome_id],
        )

    args.output_calls.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with args.output_calls.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=call_columns,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(call_rows)

    with args.input_details.open(
        encoding="utf-8",
        newline="",
    ) as handle:
        detail_reader = csv.DictReader(handle, delimiter="\t")

        if detail_reader.fieldnames is None:
            raise ValueError("Input detail table has no header.")

        detail_columns = list(detail_reader.fieldnames)
        detail_rows = list(detail_reader)

    for row in detail_rows:
        if row["gene_id"] != "spoIIQ":
            continue

        candidate = candidates.get(row["genome_id"])
        original_call = row["call"]
        final_call = merged_call(original_call, candidate)
        row["call"] = final_call

        if candidate is not None and final_call != original_call:
            hit = candidate.hit

            row["target_id"] = hit.target_id
            row["full_evalue"] = f"{hit.full_evalue:.6g}"
            row["full_score"] = f"{hit.full_score:.3f}"
            row["domain_evalue"] = f"{hit.domain_evalue:.6g}"
            row["domain_score"] = f"{hit.domain_score:.3f}"
            row["profile_coverage"] = f"{hit.profile_coverage:.3f}"

    with args.output_details.open(
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

    audit_columns = [
        "genome_id",
        "original_call",
        "final_call",
        "protein_id",
        "target_length",
        "full_evalue",
        "full_score",
        "domain_evalue",
        "domain_score",
        "profile_coverage",
        "spoIID_target",
        "spoIID_distance",
        "sequence_strength",
        "context_strength",
        "evidence_class",
    ]

    final_by_genome = {row["genome_id"]: row["spoIIQ"] for row in call_rows}

    with args.audit_output.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=audit_columns,
            delimiter="\t",
        )
        writer.writeheader()

        for genome_id in sorted(original_by_genome):
            candidate = candidates.get(genome_id)

            if candidate is None:
                writer.writerow(
                    {
                        "genome_id": genome_id,
                        "original_call": original_by_genome[genome_id],
                        "final_call": final_by_genome[genome_id],
                    }
                )
                continue

            hit = candidate.hit

            writer.writerow(
                {
                    "genome_id": genome_id,
                    "original_call": original_by_genome[genome_id],
                    "final_call": final_by_genome[genome_id],
                    "protein_id": hit.target_id,
                    "target_length": hit.target_length,
                    "full_evalue": f"{hit.full_evalue:.6g}",
                    "full_score": f"{hit.full_score:.3f}",
                    "domain_evalue": (f"{hit.domain_evalue:.6g}"),
                    "domain_score": (f"{hit.domain_score:.3f}"),
                    "profile_coverage": (f"{hit.profile_coverage:.3f}"),
                    "spoIID_target": candidate.spoiid_target,
                    "spoIID_distance": (
                        ""
                        if candidate.spoiid_distance is None
                        else candidate.spoiid_distance
                    ),
                    "sequence_strength": (candidate.sequence_strength),
                    "context_strength": (candidate.context_strength),
                    "evidence_class": candidate.evidence_class,
                }
            )

    changes = sum(
        original_by_genome[genome_id] != final_by_genome[genome_id]
        for genome_id in original_by_genome
    )

    counts: dict[str, int] = defaultdict(int)

    for candidate in candidates.values():
        if candidate is not None:
            counts[candidate.evidence_class] += 1

    print(f"Wrote: {args.output_calls}")
    print(f"Wrote: {args.output_details}")
    print(f"Wrote: {args.audit_output}")
    print(f"SpoIIQ calls changed: {changes}")

    for label in (
        "high",
        "moderate",
        "provisional",
        "uncertain",
    ):
        print(f"{label:12s}: {counts[label]}")


if __name__ == "__main__":
    main()
