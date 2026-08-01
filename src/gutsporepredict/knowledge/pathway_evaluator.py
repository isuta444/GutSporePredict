"""Evaluate biological pathways from module evaluations."""

from dataclasses import dataclass
from enum import Enum

from gutsporepredict.knowledge.evaluator import ModuleEvaluation
from gutsporepredict.knowledge.models import ReferencePathway


class PathwayStatus(str, Enum):
    """Evaluation result for one biological pathway."""

    COMPLETE = "complete"
    PARTIAL = "partial"
    ABSENT = "absent"
    UNCERTAIN = "uncertain"


@dataclass(frozen=True)
class PathwayEvaluation:
    """Evaluation result for one biological pathway."""

    pathway: ReferencePathway
    status: PathwayStatus

    completed_modules: int
    partial_modules: int
    absent_modules: int

    score: float

    module_evaluations: tuple[ModuleEvaluation, ...]
