"""Data models for sequence similarity searches."""

from dataclasses import dataclass, field
from pathlib import Path


@dataclass(frozen=True)
class SearchHit:
    """A single sequence-search hit."""

    query_id: str
    target_id: str
    identity: float
    alignment_length: int
    query_length: int
    target_length: int
    query_coverage: float
    target_coverage: float
    evalue: float
    bitscore: float
    method: str


@dataclass
class SearchResult:
    """Collection of hits produced by one search."""

    query_file: Path
    database: Path
    output_file: Path
    method: str
    hits: list[SearchHit] = field(default_factory=list)

    @property
    def hit_count(self) -> int:
        """Return the total number of hits."""

        return len(self.hits)

    @property
    def query_count(self) -> int:
        """Return the number of unique query sequences with hits."""

        return len({hit.query_id for hit in self.hits})

    def hits_for_query(self, query_id: str) -> list[SearchHit]:
        """Return hits belonging to one query sequence."""

        return [
            hit
            for hit in self.hits
            if hit.query_id == query_id
        ]

    def best_hit(self, query_id: str) -> SearchHit | None:
        """Return the highest-bitscore hit for one query."""

        query_hits = self.hits_for_query(query_id)

        if not query_hits:
            return None

        return max(query_hits, key=lambda hit: hit.bitscore)
