"""Sequence-search functionality."""

from gutsporepredict.search.base import SearchEngine
from gutsporepredict.search.diamond import (
    DIAMOND_OUTPUT_FIELDS,
    DiamondSearchEngine,
)
from gutsporepredict.search.exceptions import (
    SearchError,
    SearchExecutableNotFoundError,
    SearchExecutionError,
    SearchOutputError,
)
from gutsporepredict.search.models import (
    SearchHit,
    SearchResult,
)
from gutsporepredict.search.parser import (
    parse_diamond_output,
)

__all__ = [
    "DIAMOND_OUTPUT_FIELDS",
    "DiamondSearchEngine",
    "SearchEngine",
    "SearchError",
    "SearchExecutableNotFoundError",
    "SearchExecutionError",
    "SearchHit",
    "SearchOutputError",
    "SearchResult",
    "parse_diamond_output",
]
