"""Abstract interface for sequence-search engines."""

from abc import ABC, abstractmethod
from pathlib import Path

from gutsporepredict.search.models import SearchResult


class SearchEngine(ABC):
    """Abstract interface implemented by search engines."""

    @abstractmethod
    def search(
        self,
        query_fasta: str | Path,
        database: str | Path,
        output_file: str | Path,
    ) -> SearchResult:
        """Search query proteins against a reference database."""
