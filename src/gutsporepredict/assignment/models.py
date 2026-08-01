"""Models used by the gene-assignment subsystem."""

from dataclasses import dataclass
from enum import Enum

from gutsporepredict.reference.models import ReferenceGene
from gutsporepredict.search.models import SearchHit


class AssignmentConfidence(str, Enum):
    """Confidence assigned to a reference-gene match."""

    HIGH = "high"
    MEDIUM = "medium"
    LOW = "low"
    REJECTED = "rejected"


class AssignmentMethod(str, Enum):
    """Search method used to produce an assignment."""

    DIAMOND = "diamond"
    HMMER = "hmmer"
    MMSEQS = "mmseqs"
    CONSENSUS = "consensus"
    UNKNOWN = "unknown"


@dataclass(frozen=True)
class RuleEvaluation:
    """Result of evaluating one search hit against assignment rules."""

    confidence: AssignmentConfidence
    accepted: bool
    passed_filters: tuple[str, ...]
    failed_filters: tuple[str, ...]


@dataclass(frozen=True)
class GeneAssignment:
    """Connection between a query protein and a curated reference gene."""

    query_id: str
    reference_gene: ReferenceGene
    search_hit: SearchHit
    confidence: AssignmentConfidence
    assignment_method: AssignmentMethod
    passed_filters: tuple[str, ...]
    failed_filters: tuple[str, ...]

    @property
    def accepted(self) -> bool:
        """Return whether the assignment passed all mandatory filters."""

        return self.confidence is not AssignmentConfidence.REJECTED

    @property
    def gene_id(self) -> str:
        """Return the stable GutSporePredict gene identifier."""

        return self.reference_gene.gene_id

    @property
    def canonical_name(self) -> str:
        """Return the canonical name of the assigned reference gene."""

        return self.reference_gene.canonical_name
