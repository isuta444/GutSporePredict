"""Resolve search-hit targets against the curated reference database."""

from gutsporepredict.reference.models import (
    ReferenceDatabase,
    ReferenceGene,
)
from gutsporepredict.search.models import SearchHit


class GeneMatcher:
    """Match search-hit target identifiers to reference genes."""

    def __init__(self, database: ReferenceDatabase) -> None:
        """Initialize the matcher with a loaded reference database."""

        self._database = database

    def match(self, hit: SearchHit) -> ReferenceGene | None:
        """Return the reference gene associated with a search hit."""

        target_id = self._normalize_target_id(hit.target_id)

        gene = self._database.gene_by_id(target_id)

        if gene is not None:
            return gene

        return self._database.resolve_name(target_id)

    @staticmethod
    def _normalize_target_id(target_id: str) -> str:
        """Normalize common FASTA-header-derived target identifiers."""

        normalized = target_id.strip()

        if not normalized:
            return normalized

        normalized = normalized.split()[0]

        if "|" in normalized:
            parts = [
                part.strip()
                for part in normalized.split("|")
                if part.strip()
            ]

            for part in parts:
                if part.upper().startswith("GSP"):
                    return part

            if parts:
                return parts[-1]

        return normalized
