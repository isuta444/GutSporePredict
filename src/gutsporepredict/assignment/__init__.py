"""Gene-assignment functionality for GutSporePredict."""

from gutsporepredict.assignment.assigner import GeneAssigner
from gutsporepredict.assignment.exceptions import (
    AssignmentError,
    ReferenceGeneNotFoundError,
    UnsupportedAssignmentMethodError,
)
from gutsporepredict.assignment.matcher import GeneMatcher
from gutsporepredict.assignment.models import (
    AssignmentConfidence,
    AssignmentMethod,
    GeneAssignment,
    RuleEvaluation,
)
from gutsporepredict.assignment.rules import (
    AssignmentRuleSet,
    AssignmentThreshold,
)

__all__ = [
    "AssignmentConfidence",
    "AssignmentError",
    "AssignmentMethod",
    "AssignmentRuleSet",
    "AssignmentThreshold",
    "GeneAssigner",
    "GeneAssignment",
    "GeneMatcher",
    "ReferenceGeneNotFoundError",
    "RuleEvaluation",
    "UnsupportedAssignmentMethodError",
]
