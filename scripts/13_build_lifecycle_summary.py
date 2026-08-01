#!/usr/bin/env python3
"""Build sporulation, germination and lifecycle summaries."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

STATUS_RANK = {
    "absent": 0,
    "uncertain": 1,
    "partial": 2,
    "complete": 3,
}


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--stage-evaluations",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "--minimum-assessment",
        type=float,
        default=0.5,
    )
    return parser.parse_args()


def classify_sporulation(
    stage_rows: dict[str, dict[str, str]],
) -> str:
    """Classify the ST001-ST008 sporulation programme.

    Sporulation requires evidence of developmental continuity.
    Isolated matches to broadly conserved regulators, sigma factors,
    peptidoglycan proteins or hydrolases are not sufficient.
    """

    statuses = {
        stage_id: stage_rows[stage_id]["status"]
        for stage_id in (
            "ST001",
            "ST002",
            "ST003",
            "ST004",
            "ST005",
            "ST006",
            "ST007",
            "ST008",
        )
    }

    positive_states = {"complete", "partial"}

    entry_supported = (
        statuses["ST001"] in positive_states
        and statuses["ST002"] in positive_states
    )

    developmental_core_supported = (
        statuses["ST003"] in positive_states
        or statuses["ST004"] in positive_states
    )

    later_stage_count = sum(
        statuses[stage_id] in positive_states
        for stage_id in (
            "ST005",
            "ST006",
            "ST007",
            "ST008",
        )
    )

    if not entry_supported:
        return "absent"

    if not developmental_core_supported:
        return "absent"

    if later_stage_count < 2:
        return "partial"

    if all(
        status == "complete"
        for status in statuses.values()
    ):
        return "complete"

    if any(
        status == "partial"
        for status in statuses.values()
    ):
        return "partial"

    if any(
        status == "uncertain"
        for status in statuses.values()
    ):
        return "uncertain"

    return "partial"


def confidence_label(assessment_fraction: float) -> str:
    """Convert assessment coverage into a confidence label."""

    if assessment_fraction >= 0.8:
        return "high"

    if assessment_fraction >= 0.5:
        return "moderate"

    return "low"


def lifecycle_label(
    sporulation_status: str,
    germination_status: str,
) -> str:
    """Combine sporulation and germination calls."""

    if (
        sporulation_status == "complete"
        and germination_status == "complete"
    ):
        return "complete_lifecycle"

    if sporulation_status == "absent":
        return "no_sporulation_evidence"

    if (
        sporulation_status == "complete"
        and germination_status in {"partial", "uncertain", "absent"}
    ):
        return "sporulation_with_reduced_germination"

    if sporulation_status == "partial":
        return "incomplete_sporulation_programme"

    return "uncertain_lifecycle"


def main() -> None:
    """Build lifecycle summary table."""

    args = parse_args()

    required_columns = {
        "genome_id",
        "stage_id",
        "status",
        "score",
        "assessment_fraction",
    }

    by_genome: dict[str, dict[str, dict[str, str]]] = {}

    with args.stage_evaluations.open(
        encoding="utf-8",
        newline="",
    ) as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        if reader.fieldnames is None:
            raise ValueError("Empty stage evaluation table.")

        missing = required_columns.difference(reader.fieldnames)

        if missing:
            raise ValueError(
                "Missing columns: " + ", ".join(sorted(missing))
            )

        for row in reader:
            genome_id = row["genome_id"].strip()
            stage_id = row["stage_id"].strip()

            by_genome.setdefault(genome_id, {})[stage_id] = row

    required_stage_ids = {
        f"ST{index:03d}"
        for index in range(1, 10)
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)

    with args.output.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "genome_id",
                "sporulation_status",
                "sporulation_score",
                "sporulation_assessment",
                "sporulation_confidence",
                "germination_status",
                "germination_score",
                "germination_assessment",
                "germination_confidence",
                "lifecycle_status",
                "overall_score",
            ]
        )

        for genome_id in sorted(by_genome):
            stage_rows = by_genome[genome_id]
            missing_stages = required_stage_ids.difference(stage_rows)

            if missing_stages:
                raise ValueError(
                    f"{genome_id} is missing stages: "
                    + ", ".join(sorted(missing_stages))
                )

            sporulation_status = classify_sporulation(stage_rows)

            sporulation_scores = [
                float(stage_rows[f"ST{index:03d}"]["score"])
                for index in range(1, 9)
            ]
            sporulation_assessments = [
                float(
                    stage_rows[
                        f"ST{index:03d}"
                    ]["assessment_fraction"]
                )
                for index in range(1, 9)
            ]

            sporulation_score = (
                sum(sporulation_scores) / len(sporulation_scores)
            )
            sporulation_assessment = (
                sum(sporulation_assessments)
                / len(sporulation_assessments)
            )

            germination = stage_rows["ST009"]
            germination_status = germination["status"]
            germination_score = float(germination["score"])
            germination_assessment = float(
                germination["assessment_fraction"]
            )

            if (
                sporulation_assessment
                < args.minimum_assessment
            ):
                sporulation_status = "uncertain"

            if (
                germination_assessment
                < args.minimum_assessment
            ):
                germination_status = "uncertain"

            lifecycle_status = lifecycle_label(
                sporulation_status,
                germination_status,
            )

            overall_score = (
                sporulation_score + germination_score
            ) / 2.0

            writer.writerow(
                [
                    genome_id,
                    sporulation_status,
                    f"{sporulation_score:.3f}",
                    f"{sporulation_assessment:.3f}",
                    confidence_label(sporulation_assessment),
                    germination_status,
                    f"{germination_score:.3f}",
                    f"{germination_assessment:.3f}",
                    confidence_label(germination_assessment),
                    lifecycle_status,
                    f"{overall_score:.3f}",
                ]
            )

    print(f"Wrote: {args.output}")
    print(f"Genomes: {len(by_genome)}")


if __name__ == "__main__":
    main()
