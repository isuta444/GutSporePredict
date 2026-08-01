"""Exceptions raised by gene-assignment components."""


class AssignmentError(RuntimeError):
    """Base exception for assignment failures."""


class ReferenceGeneNotFoundError(AssignmentError):
    """Raised when a search target cannot be resolved."""


class UnsupportedAssignmentMethodError(AssignmentError):
    """Raised when a search method cannot be interpreted."""
