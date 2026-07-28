"""Exceptions raised by reference-database components."""


class ReferenceDatabaseError(RuntimeError):
    """Base exception for reference-database failures."""


class ReferenceValidationError(ReferenceDatabaseError):
    """Raised when reference data are invalid."""


class ReferenceLoadError(ReferenceDatabaseError):
    """Raised when reference data cannot be loaded."""
