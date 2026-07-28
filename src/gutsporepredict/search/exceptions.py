"""Exceptions raised by sequence-search components."""


class SearchError(RuntimeError):
    """Base exception for sequence-search failures."""


class SearchExecutableNotFoundError(SearchError):
    """Raised when an external search executable is unavailable."""


class SearchExecutionError(SearchError):
    """Raised when an external search command fails."""


class SearchOutputError(SearchError):
    """Raised when search output is missing or malformed."""
