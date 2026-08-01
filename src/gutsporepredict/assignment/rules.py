"""Rules for converting search hits into assignment confidence."""

from collections.abc import Sequence
from dataclasses import dataclass

from gutsporepredict.assignment.models import (
    AssignmentConfidence,
    RuleEvaluation,
)
from gutsporepredict.reference.models import ReferenceGene
from gutsporepredict.search.models import SearchHit


@dataclass(frozen=True)
class AssignmentThreshold:
    """Numerical thresholds for one assignment-confidence level."""

    minimum_identity: float
    minimum_query_coverage: float
    minimum_target_coverage: float
    maximum_evalue: float


HIGH_CONFIDENCE_THRESHOLD = AssignmentThreshold(
    minimum_identity=50.0,
    minimum_query_coverage=80.0,
    minimum_target_coverage=80.0,
    maximum_evalue=1e-20,
)

MEDIUM_CONFIDENCE_THRESHOLD = AssignmentThreshold(
    minimum_identity=35.0,
    minimum_query_coverage=60.0,
    minimum_target_coverage=60.0,
    maximum_evalue=1e-10,
)

LOW_CONFIDENCE_THRESHOLD = AssignmentThreshold(
    minimum_identity=25.0,
    minimum_query_coverage=40.0,
    minimum_target_coverage=40.0,
    maximum_evalue=1e-5,
)


class AssignmentRuleSet:
    """Evaluate search hits using configurable confidence thresholds."""

    def __init__(
        self,
        high: AssignmentThreshold = HIGH_CONFIDENCE_THRESHOLD,
        medium: AssignmentThreshold = MEDIUM_CONFIDENCE_THRESHOLD,
        low: AssignmentThreshold = LOW_CONFIDENCE_THRESHOLD,
    ) -> None:
        """Initialize the rule set."""

        self._thresholds: Sequence[
            tuple[AssignmentConfidence, AssignmentThreshold]
        ] = (
            (AssignmentConfidence.HIGH, high),
            (AssignmentConfidence.MEDIUM, medium),
            (AssignmentConfidence.LOW, low),
        )

        self._low_threshold = low

    def evaluate(
        self,
        hit: SearchHit,
        reference_gene: ReferenceGene,
    ) -> RuleEvaluation:
        """Evaluate a search hit."""

        # Reserved for future gene-specific rules.
        del reference_gene

        for confidence, threshold in self._thresholds:
            passed, failed = self._evaluate_threshold(
                hit,
                threshold,
            )

            if not failed:
                return RuleEvaluation(
                    confidence=confidence,
                    accepted=True,
                    passed_filters=tuple(passed),
                    failed_filters=(),
                )

        passed, failed = self._evaluate_threshold(
            hit,
            self._low_threshold,
        )

        return RuleEvaluation(
            confidence=AssignmentConfidence.REJECTED,
            accepted=False,
            passed_filters=tuple(passed),
            failed_filters=tuple(failed),
        )

    @staticmethod
    def _evaluate_threshold(
        hit: SearchHit,
        threshold: AssignmentThreshold,
    ) -> tuple[list[str], list[str]]:
        """Evaluate one threshold."""

        passed: list[str] = []
        failed: list[str] = []

        checks = (
            (
                "identity",
                hit.identity >= threshold.minimum_identity,
            ),
            (
                "query_coverage",
                hit.query_coverage
                >= threshold.minimum_query_coverage,
            ),
            (
                "target_coverage",
                hit.target_coverage
                >= threshold.minimum_target_coverage,
            ),
            (
                "evalue",
                hit.evalue <= threshold.maximum_evalue,
            ),
        )

        for name, condition in checks:
            if condition:
                passed.append(name)
            else:
                failed.append(name)

        return passed, failed
