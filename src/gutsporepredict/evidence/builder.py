"""Construction of gene-level evidence from gene assignments."""

from collections.abc import Sequence

from gutsporepredict.assignment.models import (
    AssignmentConfidence,
    GeneAssignment,
)
from gutsporepredict.evidence.models import (
    EvidenceReason,
    EvidenceStatus,
    GeneEvidence,
)
from gutsporepredict.reference.models import ReferenceGene


class EvidenceBuilder:
    """Build gene-level evidence from assignment results."""

    _confidence_rank = {
        AssignmentConfidence.HIGH: 3,
        AssignmentConfidence.MEDIUM: 2,
        AssignmentConfidence.LOW: 1,
        AssignmentConfidence.REJECTED: 0,
    }

    def build(
        self,
        reference_gene: ReferenceGene,
        assignments: Sequence[GeneAssignment],
        *,
        assessed: bool = True,
    ) -> GeneEvidence:
        """Build evidence for one reference gene.

        Assignments belonging to other reference genes are ignored.
        The best matching assignment is selected using assignment
        confidence followed by search-hit bitscore.
        """

        if not assessed:
            return GeneEvidence(
                reference_gene=reference_gene,
                status=EvidenceStatus.NOT_ASSESSED,
                reason=EvidenceReason.NOT_SEARCHED,
            )

        matching_assignments = [
            assignment
            for assignment in assignments
            if assignment.gene_id == reference_gene.gene_id
        ]

        if not matching_assignments:
            return GeneEvidence(
                reference_gene=reference_gene,
                status=EvidenceStatus.ABSENT,
                reason=EvidenceReason.NO_ACCEPTED_HIT,
            )

        best_assignment = max(
            matching_assignments,
            key=self._assignment_sort_key,
        )

        if not best_assignment.accepted:
            return GeneEvidence(
                reference_gene=reference_gene,
                status=EvidenceStatus.ABSENT,
                reason=EvidenceReason.REJECTED_ASSIGNMENT,
                assignment=best_assignment,
            )

        reason = EvidenceReason.ACCEPTED_ASSIGNMENT

        if (
            best_assignment.confidence
            is AssignmentConfidence.LOW
        ):
            reason = EvidenceReason.LOW_CONFIDENCE

        return GeneEvidence(
            reference_gene=reference_gene,
            status=EvidenceStatus.PRESENT,
            reason=reason,
            assignment=best_assignment,
        )

    @classmethod
    def _assignment_sort_key(
        cls,
        assignment: GeneAssignment,
    ) -> tuple[int, float]:
        """Return the ordering key used to select an assignment."""

        return (
            cls._confidence_rank[assignment.confidence],
            assignment.search_hit.bitscore,
        )
