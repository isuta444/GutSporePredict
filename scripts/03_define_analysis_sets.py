#!/usr/bin/env python3

"""
GutSporePredict v4.0-alpha1

Define nested comparative-genomics datasets from the selected
GTDB representative genomes.

Outputs:
    metadata/core_human_gut_set.tsv
    metadata/extended_gut_set.tsv
    metadata/evolutionary_set.tsv
    metadata/analysis_set_summary.tsv
    metadata/analysis_set_summary.json
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter
from pathlib import Path


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create nested GutSporePredict v4 analysis sets."
    )

    parser.add_argument(
        "--input",
        default="metadata/primary_genome_set.tsv",
        help="Primary genome-selection table.",
    )

    parser.add_argument(
        "--output-dir",
        default="metadata",
        help="Output directory.",
    )

    return parser.parse_args()


def truthy(value: str) -> bool:
    return str(value).strip().lower() in {
        "true",
        "t",
        "yes",
        "y",
        "1",
    }


def parse_float(value: str) -> float | None:
    value = str(value).strip()

    if not value:
        return None

    try:
        return float(value)
    except ValueError:
        return None


def write_tsv(
    path: Path,
    rows: list[dict[str, str]],
    fieldnames: list[str],
) -> None:
    with path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=fieldnames,
            extrasaction="ignore",
        )

        writer.writeheader()
        writer.writerows(rows)


def count_by_family(
    rows: list[dict[str, str]],
) -> Counter:
    return Counter(
        row.get("family", "")
        for row in rows
    )


def main() -> int:
    args = parse_arguments()

    input_path = Path(args.input)
    output_dir = Path(args.output_dir)

    if not input_path.exists():
        print(
            f"[ERROR] Input table not found: {input_path}",
            file=sys.stderr,
        )
        return 1

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    with input_path.open(
        encoding="utf-8",
        newline="",
    ) as handle:
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        if not reader.fieldnames:
            print(
                "[ERROR] Input header could not be read.",
                file=sys.stderr,
            )
            return 1

        fieldnames = list(reader.fieldnames)
        rows = list(reader)

    required_columns = {
        "gtdb_accession",
        "family",
        "quality_rank",
        "primary_selected",
        "host_origin",
        "body_site",
        "human_gut_score",
        "completeness",
        "contamination",
    }

    missing_columns = required_columns - set(fieldnames)

    if missing_columns:
        print(
            "[ERROR] Missing columns: "
            + ", ".join(sorted(missing_columns)),
            file=sys.stderr,
        )
        return 1

    evolutionary_set = []

    for row in rows:
        completeness = parse_float(
            row.get("completeness", "")
        )

        contamination = parse_float(
            row.get("contamination", "")
        )

        passes_quality = (
            completeness is not None
            and contamination is not None
            and completeness >= 90.0
            and contamination <= 5.0
        )

        if (
            truthy(row.get("primary_selected", ""))
            and passes_quality
        ):
            evolutionary_set.append(row)

    extended_gut_set = [
        row
        for row in evolutionary_set
        if row.get("body_site") == "gut"
    ]

    core_human_gut_set = [
        row
        for row in evolutionary_set
        if (
            row.get("host_origin") == "human"
            and row.get("body_site") == "gut"
        )
    ]

    strict_core_set = [
        row
        for row in core_human_gut_set
        if truthy(
            row.get("strict_quality_95_3", "")
        )
    ]

    set_definitions = {
        "core_human_gut": core_human_gut_set,
        "strict_core_human_gut": strict_core_set,
        "extended_gut": extended_gut_set,
        "evolutionary": evolutionary_set,
    }

    output_paths = {
        "core_human_gut": (
            output_dir / "core_human_gut_set.tsv"
        ),
        "strict_core_human_gut": (
            output_dir
            / "strict_core_human_gut_set.tsv"
        ),
        "extended_gut": (
            output_dir / "extended_gut_set.tsv"
        ),
        "evolutionary": (
            output_dir / "evolutionary_set.tsv"
        ),
    }

    for name, selected_rows in set_definitions.items():
        write_tsv(
            output_paths[name],
            selected_rows,
            fieldnames,
        )

    all_families = sorted({
        row.get("family", "")
        for row in evolutionary_set
        if row.get("family", "")
    })

    summary_rows = []

    for family in all_families:
        summary_row = {
            "family": family,
        }

        for set_name, selected_rows in set_definitions.items():
            summary_row[set_name] = sum(
                row.get("family") == family
                for row in selected_rows
            )

        summary_rows.append(summary_row)

    summary_fieldnames = [
        "family",
        "core_human_gut",
        "strict_core_human_gut",
        "extended_gut",
        "evolutionary",
    ]

    summary_tsv_path = (
        output_dir / "analysis_set_summary.tsv"
    )

    write_tsv(
        summary_tsv_path,
        summary_rows,
        summary_fieldnames,
    )

    summary_json = {
        "pipeline": "GutSporePredict",
        "version": "4.0-alpha1",
        "definitions": {
            "core_human_gut": (
                "GTDB representative; completeness >=90; "
                "contamination <=5; human and gut metadata"
            ),
            "strict_core_human_gut": (
                "Core human-gut set additionally satisfying "
                "completeness >=95 and contamination <=3"
            ),
            "extended_gut": (
                "GTDB representative; completeness >=90; "
                "contamination <=5; gut-associated metadata "
                "regardless of identified host"
            ),
            "evolutionary": (
                "All target-family GTDB representatives "
                "with completeness >=90 and contamination <=5"
            ),
        },
        "counts": {
            name: len(selected_rows)
            for name, selected_rows
            in set_definitions.items()
        },
        "family_counts": {
            name: dict(
                sorted(
                    count_by_family(
                        selected_rows
                    ).items()
                )
            )
            for name, selected_rows
            in set_definitions.items()
        },
    }

    summary_json_path = (
        output_dir / "analysis_set_summary.json"
    )

    summary_json_path.write_text(
        json.dumps(
            summary_json,
            indent=2,
            ensure_ascii=False,
        )
        + "\n",
        encoding="utf-8",
    )

    print("=" * 72)
    print("GutSporePredict v4.0-alpha1")
    print("Analysis-set definition")
    print("=" * 72)

    for name, selected_rows in set_definitions.items():
        print(
            f"{name:25s} "
            f"{len(selected_rows):>6,} genomes"
        )

    print()
    print("[Family counts]")

    for family in all_families:
        values = {
            name: count_by_family(
                selected_rows
            )[family]
            for name, selected_rows
            in set_definitions.items()
        }

        print(
            f"{family:25s} "
            f"core={values['core_human_gut']:>4,}  "
            f"strict={values['strict_core_human_gut']:>4,}  "
            f"gut={values['extended_gut']:>4,}  "
            f"evolutionary={values['evolutionary']:>4,}"
        )

    print()
    for path in output_paths.values():
        print(f"[OUTPUT] {path}")

    print(f"[OUTPUT] {summary_tsv_path}")
    print(f"[OUTPUT] {summary_json_path}")
    print("=" * 72)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
