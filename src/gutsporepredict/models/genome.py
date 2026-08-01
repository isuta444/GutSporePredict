"""Genome data model."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(slots=True)
class Genome:
    """Genome metadata used throughout GutSporePredict."""

    accession: str
    fasta: Path

    @property
    def name(self) -> str:
        """Return the genome name derived from the FASTA filename."""
        return self.fasta.stem
