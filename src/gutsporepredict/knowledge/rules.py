"""Rule models for biological inference."""

from dataclasses import dataclass

from gutsporepredict.knowledge.models import (
    KnowledgeEvidenceLevel,
)


@dataclass(frozen=True)
class ReferenceRule:
    """One biological inference rule."""

    rule_id: str
    module_id: str

    condition: str
    outcome: str

    evidence_level: KnowledgeEvidenceLevel

    literature: tuple[str, ...] = ()

    notes: str | None = None

    def __post_init__(self) -> None:
        """Validate the rule."""

        if not self.rule_id.strip():
            raise ValueError(
                "ReferenceRule.rule_id must not be empty."
            )

        if not self.module_id.strip():
            raise ValueError(
                "ReferenceRule.module_id must not be empty."
            )

        if not self.condition.strip():
            raise ValueError(
                "ReferenceRule.condition must not be empty."
            )

        if not self.outcome.strip():
            raise ValueError(
                "ReferenceRule.outcome must not be empty."
            )
