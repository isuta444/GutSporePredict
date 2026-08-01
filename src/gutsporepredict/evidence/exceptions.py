"""Exceptions raised by the evidence subsystem."""


class EvidenceError(Exception):
    """Base exception for evidence-related errors."""


class EvidenceBuildError(EvidenceError):
    """Raised when gene evidence cannot be constructed."""
