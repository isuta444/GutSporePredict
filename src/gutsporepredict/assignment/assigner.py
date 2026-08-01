"""Convert search hits into curated gene assignments."""

from collections.abc import Iterable

from gutsporepredict.assignment.exceptions import (
    ReferenceGeneNotFoundError,
)
from gutsporepredict.assignment.matcher import GeneMatcher
from gutsporepredict.assignment.models import (
    AssignmentConfidence,
    AssignmentMethod,
    GeneAssignment,
)
from gutsporepredict.assignment.rules import AssignmentRuleSet
from gutsporepredict.reference.models import ReferenceDatabase
from gutsporepredict.search.models import SearchHit

_CONFIDENCE_RANK = {
    AssignmentConfidence.REJECTED: 0,
    AssignmentConfidence.LOW: 1,
    AssignmentConfidence.MEDIUM: 2,
    AssignmentConfidence.HIGH: 3,
}


class GeneAssigner:
    """Assign search hits to curated reference genes."""

    def __init__(
        self,
        database: ReferenceDatabase,
        rules: AssignmentRuleSet | None = None,
    ) -> None:
        """Initialize the assigner."""

        self._matcher = GeneMatcher(database)
        self._rules = rules or AssignmentRuleSet()

    def assign(
        self,
        hit: SearchHit,
        *,
        strict_reference_matching: bool = False,
    ) -> GeneAssignment | None:
        """Convert a single search hit into a gene assignment.

        When strict_reference_matching is false, unresolved targets are
        ignored and None is returned. When true, an exception is raised.
        """

        reference_gene = self._matcher.match(hit)

        if reference_gene is None:
            if strict_reference_matching:
                raise ReferenceGeneNotFoundError(
                    "Search target could not be resolved: "
                    f"{hit.target_id}"
                )

            return None

        evaluation = self._rules.evaluate(
            hit,
            reference_gene,
        )

        return GeneAssignment(
            query_id=hit.query_id,
            reference_gene=reference_gene,
            search_hit=hit,
            confidence=evaluation.confidence,
            assignment_method=self._assignment_method(hit.method),
            passed_filters=evaluation.passed_filters,
            failed_filters=evaluation.failed_filters,
        )

    def assign_all(
        self,
        hits: Iterable[SearchHit],
        *,
        include_rejected: bool = False,
        strict_reference_matching: bool = False,
    ) -> list[GeneAssignment]:
        """Assign all resolvable search hits."""

        assignments: list[GeneAssignment] = []

        for hit in hits:
            assignment = self.assign(
                hit,
                strict_reference_matching=strict_reference_matching,
            )

            if assignment is None:
                continue

            if not include_rejected and not assignment.accepted:
                continue

            assignments.append(assignment)

        return assignments

    def best_assignments(
        self,
        hits: Iterable[SearchHit],
        *,
        include_rejected: bool = False,
    ) -> list[GeneAssignment]:
        """Return the best assignment for each query protein."""

        assignments = self.assign_all(
            hits,
            include_rejected=include_rejected,
        )

        best_by_query: dict[str, GeneAssignment] = {}

        for assignment in assignments:
            current = best_by_query.get(assignment.query_id)

            if current is None or self._is_better(
                assignment,
                current,
            ):
                best_by_query[assignment.query_id] = assignment

        return [
            best_by_query[query_id]
            for query_id in sorted(best_by_query)
        ]

    @staticmethod
    def _assignment_method(method: str) -> AssignmentMethod:
        normalized = method.strip().lower()

        try:
            return AssignmentMethod(normalized)
        except ValueError:
            return AssignmentMethod.UNKNOWN

    @staticmethod
    def _is_better(
        candidate: GeneAssignment,
        current: GeneAssignment,
    ) -> bool:
        candidate_key = (
            _CONFIDENCE_RANK[candidate.confidence],
            candidate.search_hit.bitscore,
            -candidate.search_hit.evalue,
            candidate.search_hit.identity,
            candidate.search_hit.query_coverage,
            candidate.search_hit.target_coverage,
        )
        current_key = (
            _CONFIDENCE_RANK[current.confidence],
            current.search_hit.bitscore,
            -current.search_hit.evalue,
            current.search_hit.identity,
            current.search_hit.query_coverage,
            current.search_hit.target_coverage,
        )

        return candidate_key > current_key
