"""Evaluate biological modules from gene-level evidence."""

from dataclasses import dataclass
from enum import Enum

from gutsporepredict.evidence.models import GeneEvidence
from gutsporepredict.knowledge.models import (
    GeneRequirement,
    ReferenceModule,
)


class ModuleStatus(str, Enum):
    """Evaluation result for one biological module."""

    COMPLETE = "complete"
    PARTIAL = "partial"
    ABSENT = "absent"
    UNCERTAIN = "uncertain"


@dataclass(frozen=True)
class ModuleEvaluation:
    """Evaluation result for one reference module."""

    module: ReferenceModule
    status: ModuleStatus

    required_present: int
    required_total: int

    supporting_present: int
    supporting_total: int

    optional_present: int
    optional_total: int

    score: float

    gene_evidence: tuple[GeneEvidence, ...]


class ModuleEvaluator:
    """Evaluate one biological module."""

    def evaluate(
        self,
        module: ReferenceModule,
        evidence: tuple[GeneEvidence, ...],
    ) -> ModuleEvaluation:
        """Evaluate a reference module."""

        evidence_by_gene = {
            item.gene_id: item
            for item in evidence
        }

        required_total = 0
        required_assessed = 0
        required_present = 0
        required_uncertain = 0

        supporting_total = 0
        supporting_present = 0

        optional_total = 0
        optional_present = 0

        collected: list[GeneEvidence] = []

        assessed_weight = 0.0
        observed_weight = 0.0

        for gene in module.genes:
            gene_evidence = evidence_by_gene.get(gene.gene_id)

            if gene_evidence is not None:
                collected.append(gene_evidence)

            assessed = (
                gene_evidence is not None
                and gene_evidence.assessed
            )
            present = (
                gene_evidence is not None
                and gene_evidence.present
            )

            if assessed:
                assessed_weight += gene.weight

            if present:
                observed_weight += gene.weight

            if gene.requirement is GeneRequirement.REQUIRED:
                required_total += 1

                if assessed:
                    required_assessed += 1

                if present:
                    required_present += 1

                if (
                    gene_evidence is not None
                    and gene_evidence.status.value == "uncertain"
                ):
                    required_uncertain += 1

            elif gene.requirement is GeneRequirement.SUPPORTING:
                supporting_total += 1

                if present:
                    supporting_present += 1

            else:
                optional_total += 1

                if present:
                    optional_present += 1

        if required_total == 0:
            status = ModuleStatus.UNCERTAIN
        elif required_assessed < required_total:
            status = ModuleStatus.UNCERTAIN
        elif required_present == required_total:
            status = ModuleStatus.COMPLETE
        elif required_uncertain > 0:
            status = ModuleStatus.UNCERTAIN
        elif required_present == 0:
            status = ModuleStatus.ABSENT
        else:
            status = ModuleStatus.PARTIAL

        score = (
            observed_weight / assessed_weight
            if assessed_weight > 0
            else 0.0
        )

        return ModuleEvaluation(
            module=module,
            status=status,
            required_present=required_present,
            required_total=required_total,
            supporting_present=supporting_present,
            supporting_total=supporting_total,
            optional_present=optional_present,
            optional_total=optional_total,
            score=score,
            gene_evidence=tuple(collected),
        )
