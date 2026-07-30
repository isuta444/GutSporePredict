"""Exceptions raised by the knowledge-base subsystem."""

from pathlib import Path


class KnowledgeBaseError(Exception):
    """Base exception for knowledge-base errors."""


class KnowledgeModelError(KnowledgeBaseError):
    """Raised when a knowledge-base model is invalid."""


class KnowledgeLoadError(KnowledgeBaseError):
    """Base exception for knowledge-file loading errors."""


class KnowledgeFileNotFoundError(KnowledgeLoadError):
    """Raised when a requested knowledge file does not exist."""

class KnowledgeDatabaseError(KnowledgeBaseError):
    """Raised when a knowledge database is invalid."""

    def __init__(self, path: str | Path) -> None:
        """Initialize the exception with the missing file path."""

        self.path = Path(path)
        super().__init__(
            f"Knowledge file does not exist: {self.path}"
        )


class KnowledgeFormatError(KnowledgeLoadError):
    """Raised when a knowledge file has an invalid format."""

    def __init__(
        self,
        path: str | Path,
        message: str,
    ) -> None:
        """Initialize the exception with file and format information."""

        self.path = Path(path)
        self.message = message
        super().__init__(
            f"Invalid knowledge file {self.path}: {message}"
        )
