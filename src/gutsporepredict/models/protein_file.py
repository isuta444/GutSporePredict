"""Protein FASTA input model."""

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class ProteinFile:
    """Protein FASTA file used as a GutSporePredict input."""

    accession: str
    fasta: Path
