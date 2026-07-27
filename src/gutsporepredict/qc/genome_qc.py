"""Genome quality control."""

from dataclasses import dataclass

from gutsporepredict.io.fasta_parser import FastaRecord


@dataclass
class GenomeQCResult:
    """Genome assembly statistics."""

    genome_size: int
    contig_count: int
    gc_content: float
    n50: int

class GenomeQC:
    """Calculate genome assembly statistics."""

    @staticmethod
    def calculate(records: list[FastaRecord]) -> GenomeQCResult:
        """Calculate genome QC metrics."""

        genome_size = GenomeQC.genome_size(records)
        contig_count = GenomeQC.contig_count(records)
        gc_content = GenomeQC.gc_content(records)
        n50 = GenomeQC.n50(records)

        return GenomeQCResult(
            genome_size=genome_size,
            contig_count=contig_count,
            gc_content=gc_content,
            n50=n50,
        )

    @staticmethod
    def genome_size(records: list[FastaRecord]) -> int:
        """Calculate genome size."""

        return sum(len(record.sequence) for record in records)

    @staticmethod
    def contig_count(records: list[FastaRecord]) -> int:
        """Count contigs."""

        return len(records)

    @staticmethod
    def gc_content(records: list[FastaRecord]) -> float:
        """Calculate GC content."""

        sequence = "".join(record.sequence.upper() for record in records)

        if not sequence:
            return 0.0

        gc = sequence.count("G") + sequence.count("C")

        return gc / len(sequence) * 100

    @staticmethod
    def n50(records: list[FastaRecord]) -> int:
        """Calculate N50."""

        lengths = sorted(
            (len(record.sequence) for record in records),
            reverse=True,
        )

        total = sum(lengths)
        threshold = total / 2

        cumulative = 0

        for length in lengths:

            cumulative += length

            if cumulative >= threshold:
                return length

        return 0
        