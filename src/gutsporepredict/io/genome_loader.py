"""Genome loader."""

from pathlib import Path

from gutsporepredict.models.genome import Genome

class GenomeLoader:

    """Load genome FASTA files from a directory."""

    VALID_SUFFIXES = {".fa", ".fna", ".fasta"}

    def __init__(self, genome_dir: str | Path):

        self.genome_dir = Path(genome_dir)

    def load(self) -> list[Genome]:

        genomes = []

        for fasta in sorted(self.genome_dir.iterdir()):

            if fasta.suffix.lower() not in self.VALID_SUFFIXES:

                continue

            genomes.append(

                Genome(

                    accession=fasta.stem,

                    fasta=fasta,

                )

            )

        return genomes