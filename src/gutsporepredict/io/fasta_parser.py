"""FASTA parser."""

from dataclasses import dataclass
from pathlib import Path


@dataclass
class FastaRecord:
    """Single FASTA record."""

    header: str
    sequence: str


class FastaParser:
    """Parser for FASTA files."""

    @staticmethod
    def parse(fasta_file: str | Path) -> list[FastaRecord]:
        """Parse a FASTA file."""

        fasta_file = Path(fasta_file)

        records: list[FastaRecord] = []

        header: str | None = None
        sequence: list[str] = []

        with fasta_file.open() as handle:

            for line in handle:

                line = line.strip()

                if not line:
                    continue

                if line.startswith(">"):

                    if header is not None:

                        records.append(
                            FastaRecord(
                                header=header,
                                sequence="".join(sequence),
                            )
                        )

                    header = line[1:]
                    sequence = []

                else:

                    sequence.append(line)

        if header is not None:

            records.append(
                FastaRecord(
                    header=header,
                    sequence="".join(sequence),
                )
            )

        return records