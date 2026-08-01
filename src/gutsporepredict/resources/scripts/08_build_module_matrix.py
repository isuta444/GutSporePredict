#!/usr/bin/env python3
"""Evaluate biological modules from a gene presence/absence matrix."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from gutsporepredict.evidence.models import (
    EvidenceReason,
    EvidenceStatus,
    GeneEvidence,
)
from gutsporepredict.knowledge.evaluator import ModuleEvaluator
from gutsporepredict.knowledge.knowledge_registry import (
    KnowledgeRegistry,
)
from gutsporepredict.reference.models import ReferenceGene


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(
        description=(
            "Convert a gene presence/absence matrix into module "
            "evaluation tables."
        )
    )
    parser.add_argument(
        "--presence-matrix",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "--knowledge-dir",
        type=Path,
        default=Path("knowledge"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
    )
    return parser.parse_args()


def load_presence_matrix(
    path: Path,
) -> tuple[list[str], dict[str, dict[str, str]]]:
    """Load genome-by-gene three-state calls."""

    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        if reader.fieldnames is None:
            raise ValueError(f"Empty presence matrix: {path}")

        if "genome_id" not in reader.fieldnames:
            raise ValueError(
                f"Presence matrix must contain genome_id: {path}"
            )

        gene_ids = [
            column
            for column in reader.fieldnames
            if column != "genome_id"
        ]

        calls: dict[str, dict[str, str]] = {}
        allowed_calls = {"0", "1", "?"}

        for row in reader:
            genome_id = row["genome_id"].strip()

            if not genome_id:
                raise ValueError(
                    f"Empty genome_id in presence matrix: {path}"
                )

            genome_calls: dict[str, str] = {}

            for gene_id in gene_ids:
                call = row[gene_id].strip()

                if call not in allowed_calls:
                    raise ValueError(
                        f"Invalid gene call '{call}' for "
                        f"{genome_id}/{gene_id}. "
                        "Expected one of: 0, 1, ?"
                    )

                genome_calls[gene_id] = call

            calls[genome_id] = genome_calls

    return gene_ids, calls


def make_reference_gene(
    gene_id: str,
    module_id: str,
    pathway_id: str,
) -> ReferenceGene:
    """Create the minimal reference-gene object used by GeneEvidence."""

    return ReferenceGene(
        gene_id=gene_id,
        canonical_name=gene_id,
        pathway=pathway_id,
        module=module_id,
        stage="unknown",
        essentiality="unknown",
        phyletic_pattern="unknown",
        search_methods=("hmmer",),
        description="HMMER-derived presence/absence evidence.",
    )


def main() -> None:
    """Evaluate every loaded module for every genome."""

    args = parse_args()

    searched_gene_ids, presence_by_genome = load_presence_matrix(
        args.presence_matrix
    )
    searched_gene_id_set = set(searched_gene_ids)
    registry = KnowledgeRegistry.load(args.knowledge_dir)
    evaluator = ModuleEvaluator()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    long_path = args.output_dir / "module_evaluations.tsv"
    status_path = args.output_dir / "module_status_matrix.tsv"
    score_path = args.output_dir / "module_score_matrix.tsv"

    modules = [
        registry.modules[module_id]
        for module_id in sorted(registry.modules)
    ]

    results_by_genome = {}

    for genome_id, gene_calls in presence_by_genome.items():
        genome_results = {}

        for module in modules:
            evidence = []

            for module_gene in module.genes:
                reference_gene = make_reference_gene(
                    module_gene.gene_id,
                    module.module_id,
                    module.pathway_id,
                )

                if module_gene.gene_id not in searched_gene_id_set:
                    evidence.append(
                        GeneEvidence(
                            reference_gene=reference_gene,
                            status=EvidenceStatus.NOT_ASSESSED,
                            reason=EvidenceReason.NOT_SEARCHED,
                        )
                    )
                    continue

                call = gene_calls[module_gene.gene_id]

                if call == "1":
                    status = EvidenceStatus.PRESENT
                    reason = EvidenceReason.ACCEPTED_ASSIGNMENT
                elif call == "?":
                    status = EvidenceStatus.UNCERTAIN
                    reason = EvidenceReason.LOW_CONFIDENCE
                else:
                    status = EvidenceStatus.ABSENT
                    reason = EvidenceReason.NO_ACCEPTED_HIT

                evidence.append(
                    GeneEvidence(
                        reference_gene=reference_gene,
                        status=status,
                        reason=reason,
                    )
                )

            genome_results[module.module_id] = evaluator.evaluate(
                module,
                tuple(evidence),
            )

        results_by_genome[genome_id] = genome_results

    with long_path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "genome_id",
                "module_id",
                "module_name",
                "status",
                "score",
                "required_present",
                "required_total",
                "supporting_present",
                "supporting_total",
                "optional_present",
                "optional_total",
                "genes_assessed",
                "genes_total",
                "assessment_fraction",
            ]
        )

        for genome_id in sorted(results_by_genome):
            for module in modules:
                result = results_by_genome[genome_id][
                    module.module_id
                ]

                genes_assessed = sum(
                    evidence.assessed
                    for evidence in result.gene_evidence
                )
                genes_total = len(result.module.genes)
                assessment_fraction = genes_assessed / genes_total

                writer.writerow(
                    [
                        genome_id,
                        module.module_id,
                        module.name,
                        result.status.value,
                        f"{result.score:.3f}",
                        result.required_present,
                        result.required_total,
                        result.supporting_present,
                        result.supporting_total,
                        result.optional_present,
                        result.optional_total,
                        genes_assessed,
                        genes_total,
                        f"{assessment_fraction:.3f}",
                    ]
                )

    with status_path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "genome_id",
                *[module.module_id for module in modules],
            ]
        )

        for genome_id in sorted(results_by_genome):
            writer.writerow(
                [
                    genome_id,
                    *[
                        results_by_genome[genome_id][
                            module.module_id
                        ].status.value
                        for module in modules
                    ],
                ]
            )

    with score_path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "genome_id",
                *[module.module_id for module in modules],
            ]
        )

        for genome_id in sorted(results_by_genome):
            writer.writerow(
                [
                    genome_id,
                    *[
                        f"{results_by_genome[genome_id][module.module_id].score:.3f}"
                        for module in modules
                    ],
                ]
            )

    print(f"Wrote: {long_path}")
    print(f"Wrote: {status_path}")
    print(f"Wrote: {score_path}")
    print(
        f"Genomes: {len(results_by_genome)}, "
        f"modules: {len(modules)}"
    )


if __name__ == "__main__":
    main()
