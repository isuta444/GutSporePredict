#!/usr/bin/env python3
"""Aggregate module evaluations into sporulation stages."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import yaml

VALID_MODULE_STATUSES = {
    "complete",
    "partial",
    "absent",
    "uncertain",
}


@dataclass(frozen=True)
class StageDefinition:
    """Definition of one developmental stage."""

    stage_id: str
    name: str
    description: str
    module_ids: tuple[str, ...]
    aggregation: str = "all"


@dataclass(frozen=True)
class ModuleResult:
    """Module evaluation used for stage aggregation."""

    module_id: str
    status: str
    score: float
    genes_assessed: int
    genes_total: int


@dataclass(frozen=True)
class StageResult:
    """Aggregated result for one developmental stage."""

    stage: StageDefinition
    status: str
    score: float
    genes_assessed: int
    genes_total: int
    assessment_fraction: float
    module_results: tuple[ModuleResult, ...]


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(
        description=(
            "Aggregate GutSporePredict module evaluations into "
            "developmental-stage evaluations."
        )
    )
    parser.add_argument(
        "--module-evaluations",
        type=Path,
        required=True,
        help="Long-format module_evaluations.tsv file.",
    )
    parser.add_argument(
        "--stage-definitions",
        type=Path,
        default=Path("knowledge/stages.yaml"),
        help="YAML file defining stages and their modules.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory for stage evaluation outputs.",
    )
    return parser.parse_args()


def required_string(
    record: dict[str, object],
    field: str,
) -> str:
    """Return a required non-empty string."""

    value = record.get(field)

    if not isinstance(value, str) or not value.strip():
        raise ValueError(
            f"Stage field '{field}' must be a non-empty string."
        )

    return value.strip()


def load_stage_definitions(path: Path) -> list[StageDefinition]:
    """Load developmental-stage definitions from YAML."""

    with path.open(encoding="utf-8") as handle:
        data = yaml.safe_load(handle)

    if not isinstance(data, dict):
        raise ValueError(
            f"Stage YAML root must be a mapping: {path}"
        )

    raw_stages = data.get("stages")

    if not isinstance(raw_stages, list) or not raw_stages:
        raise ValueError(
            f"Stage YAML must contain a non-empty stages list: {path}"
        )

    stages: list[StageDefinition] = []
    seen_stage_ids: set[str] = set()

    for index, raw_stage in enumerate(raw_stages):
        if not isinstance(raw_stage, dict):
            raise ValueError(
                f"stages[{index}] must be a mapping."
            )

        stage_id = required_string(raw_stage, "stage_id")
        name = required_string(raw_stage, "name")
        description = required_string(
            raw_stage,
            "description",
        )
        raw_module_ids = raw_stage.get("module_ids")

        if (
            not isinstance(raw_module_ids, list)
            or not raw_module_ids
            or not all(
                isinstance(module_id, str)
                and module_id.strip()
                for module_id in raw_module_ids
            )
        ):
            raise ValueError(
                f"Stage '{stage_id}' must contain module_ids."
            )

        module_ids = tuple(
            module_id.strip()
            for module_id in raw_module_ids
        )

        if len(module_ids) != len(set(module_ids)):
            raise ValueError(
                f"Stage '{stage_id}' contains duplicate modules."
            )

        if stage_id in seen_stage_ids:
            raise ValueError(
                f"Duplicate stage identifier: {stage_id}"
            )

        seen_stage_ids.add(stage_id)

        aggregation_value = raw_stage.get(
            "aggregation",
            "all",
        )

        if (
            not isinstance(aggregation_value, str)
            or aggregation_value not in {"all", "any"}
        ):
            raise ValueError(
                f"Stage '{stage_id}' aggregation must be "
                "'all' or 'any'."
            )

        stages.append(
            StageDefinition(
                stage_id=stage_id,
                name=name,
                description=description,
                module_ids=module_ids,
                aggregation=aggregation_value,
            )
        )

    return stages


def load_module_evaluations(
    path: Path,
) -> dict[str, dict[str, ModuleResult]]:
    """Load module evaluations indexed by genome and module."""

    required_columns = {
        "genome_id",
        "module_id",
        "status",
        "score",
        "genes_assessed",
        "genes_total",
    }

    results: dict[str, dict[str, ModuleResult]] = {}

    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        if reader.fieldnames is None:
            raise ValueError(
                f"Empty module evaluation table: {path}"
            )

        missing_columns = required_columns.difference(
            reader.fieldnames
        )

        if missing_columns:
            raise ValueError(
                "Module evaluation table is missing columns: "
                + ", ".join(sorted(missing_columns))
            )

        for row in reader:
            genome_id = row["genome_id"].strip()
            module_id = row["module_id"].strip()
            status = row["status"].strip()

            if not genome_id or not module_id:
                raise ValueError(
                    "genome_id and module_id must not be empty."
                )

            if status not in VALID_MODULE_STATUSES:
                raise ValueError(
                    f"Invalid module status: {status}"
                )

            module_result = ModuleResult(
                module_id=module_id,
                status=status,
                score=float(row["score"]),
                genes_assessed=int(row["genes_assessed"]),
                genes_total=int(row["genes_total"]),
            )

            genome_results = results.setdefault(
                genome_id,
                {},
            )

            if module_id in genome_results:
                raise ValueError(
                    f"Duplicate module result: "
                    f"{genome_id}/{module_id}"
                )

            genome_results[module_id] = module_result

    if not results:
        raise ValueError(
            f"No module evaluations found: {path}"
        )

    return results


def determine_stage_status(
    module_results: tuple[ModuleResult, ...],
    aggregation: str,
) -> str:
    """Determine a stage status from its component modules."""

    statuses = [
        result.status
        for result in module_results
    ]

    if aggregation == "any":
        if "complete" in statuses:
            return "complete"

        if "partial" in statuses:
            return "partial"

        if "uncertain" in statuses:
            return "uncertain"

        return "absent"

    if all(status == "complete" for status in statuses):
        return "complete"

    if "absent" in statuses:
        return "absent"

    if "uncertain" in statuses:
        return "uncertain"

    return "partial"


def evaluate_stage(
    stage: StageDefinition,
    module_results: tuple[ModuleResult, ...],
) -> StageResult:
    """Aggregate component modules into one stage result."""

    genes_assessed = sum(
        result.genes_assessed
        for result in module_results
    )
    genes_total = sum(
        result.genes_total
        for result in module_results
    )

    assessment_fraction = (
        genes_assessed / genes_total
        if genes_total > 0
        else 0.0
    )

    if stage.aggregation == "any":
        score = max(
            (
                result.score
                for result in module_results
                if result.genes_assessed > 0
            ),
            default=0.0,
        )
    else:
        score_weight = sum(
            result.genes_assessed
            for result in module_results
        )

        score = (
            sum(
                result.score * result.genes_assessed
                for result in module_results
            )
            / score_weight
            if score_weight > 0
            else 0.0
        )

    return StageResult(
        stage=stage,
        status=determine_stage_status(
            module_results,
            stage.aggregation,
        ),
        score=score,
        genes_assessed=genes_assessed,
        genes_total=genes_total,
        assessment_fraction=assessment_fraction,
        module_results=module_results,
    )


def main() -> None:
    """Build long- and wide-format stage evaluation tables."""

    args = parse_args()
    stages = load_stage_definitions(
        args.stage_definitions
    )
    module_results_by_genome = load_module_evaluations(
        args.module_evaluations
    )

    stage_results_by_genome: dict[
        str,
        dict[str, StageResult],
    ] = {}

    for genome_id, genome_modules in (
        module_results_by_genome.items()
    ):
        genome_stage_results: dict[str, StageResult] = {}

        for stage in stages:
            missing_modules = [
                module_id
                for module_id in stage.module_ids
                if module_id not in genome_modules
            ]

            if missing_modules:
                raise ValueError(
                    f"Genome '{genome_id}' is missing module results: "
                    + ", ".join(missing_modules)
                )

            component_results = tuple(
                genome_modules[module_id]
                for module_id in stage.module_ids
            )

            genome_stage_results[stage.stage_id] = (
                evaluate_stage(
                    stage,
                    component_results,
                )
            )

        stage_results_by_genome[genome_id] = (
            genome_stage_results
        )

    args.output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    long_path = args.output_dir / "stage_evaluations.tsv"
    status_path = args.output_dir / "stage_status_matrix.tsv"
    score_path = args.output_dir / "stage_score_matrix.tsv"
    assessment_path = (
        args.output_dir / "stage_assessment_matrix.tsv"
    )

    with long_path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "genome_id",
                "stage_id",
                "stage_name",
                "module_ids",
                "status",
                "score",
                "genes_assessed",
                "genes_total",
                "assessment_fraction",
            ]
        )

        for genome_id in sorted(stage_results_by_genome):
            for stage in stages:
                result = stage_results_by_genome[
                    genome_id
                ][stage.stage_id]

                writer.writerow(
                    [
                        genome_id,
                        stage.stage_id,
                        stage.name,
                        ",".join(stage.module_ids),
                        result.status,
                        f"{result.score:.3f}",
                        result.genes_assessed,
                        result.genes_total,
                        f"{result.assessment_fraction:.3f}",
                    ]
                )

    def write_matrix(
        path: Path,
        value_name: str,
    ) -> None:
        with path.open(
            "w",
            encoding="utf-8",
            newline="",
        ) as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(
                [
                    "genome_id",
                    *[
                        stage.stage_id
                        for stage in stages
                    ],
                ]
            )

            for genome_id in sorted(
                stage_results_by_genome
            ):
                values = []

                for stage in stages:
                    result = stage_results_by_genome[
                        genome_id
                    ][stage.stage_id]

                    if value_name == "status":
                        value = result.status
                    elif value_name == "score":
                        value = f"{result.score:.3f}"
                    else:
                        value = (
                            f"{result.assessment_fraction:.3f}"
                        )

                    values.append(value)

                writer.writerow(
                    [
                        genome_id,
                        *values,
                    ]
                )

    write_matrix(status_path, "status")
    write_matrix(score_path, "score")
    write_matrix(assessment_path, "assessment")

    print(f"Wrote: {long_path}")
    print(f"Wrote: {status_path}")
    print(f"Wrote: {score_path}")
    print(f"Wrote: {assessment_path}")
    print(
        f"Genomes: {len(stage_results_by_genome)}, "
        f"stages: {len(stages)}"
    )


if __name__ == "__main__":
    main()
