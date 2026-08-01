"""Data models for gene-level evidence evaluation."""

from dataclasses import dataclass
from enum import Enum

from gutsporepredict.assignment.models import (
    AssignmentConfidence,
    GeneAssignment,
)
from gutsporepredict.reference.models import ReferenceGene


class EvidenceStatus(str, Enum):
    """Assessment status for one reference gene."""

    PRESENT = "present"
    ABSENT = "absent"
    UNCERTAIN = "uncertain"
    NOT_ASSESSED = "not_assessed"
    MISSING_DUE_TO_GENOME_QUALITY = (
        "missing_due_to_genome_quality"
    )


class EvidenceReason(str, Enum):
    """Reason supporting a gene-evidence status."""

    ACCEPTED_ASSIGNMENT = "accepted_assignment"
    REJECTED_ASSIGNMENT = "rejected_assignment"
    NO_ACCEPTED_HIT = "no_accepted_hit"
    LOW_CONFIDENCE = "low_confidence"
    NOT_SEARCHED = "not_searched"
    GENOME_FRAGMENTED = "genome_fragmented"
    GENOME_INCOMPLETE = "genome_incomplete"
    UNKNOWN = "unknown"


@dataclass(frozen=True)
class GeneEvidence:
    """Evidence supporting the status of one reference gene."""

    reference_gene: ReferenceGene
    status: EvidenceStatus
    reason: EvidenceReason
    assignment: GeneAssignment | None = None
    notes: str | None = None

    @property
    def present(self) -> bool:
        """Return whether the reference gene is considered present."""

        return self.status is EvidenceStatus.PRESENT

    @property
    def assessed(self) -> bool:
        """Return whether the reference gene was assessed."""

        return self.status is not EvidenceStatus.NOT_ASSESSED

    @property
    def gene_id(self) -> str:
        """Return the stable GutSporePredict gene identifier."""

        return self.reference_gene.gene_id

    @property
    def canonical_name(self) -> str:
        """Return the canonical reference-gene name."""

        return self.reference_gene.canonical_name

    @property
    def confidence(self) -> AssignmentConfidence | None:
        """Return assignment confidence when an assignment exists."""

        if self.assignment is None:
            return None

        return self.assignment.confidence
