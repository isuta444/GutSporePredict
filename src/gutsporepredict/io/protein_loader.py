"""Protein FASTA loader."""

from pathlib import Path

from gutsporepredict.models.protein_file import ProteinFile


class ProteinLoader:
    """Load protein FASTA files from a directory."""

    VALID_SUFFIXES = {".faa"}

    def __init__(self, protein_dir: str | Path) -> None:
        self.protein_dir = Path(protein_dir)

    def load(self) -> list[ProteinFile]:
        """Load protein FASTA files."""

        protein_files: list[ProteinFile] = []

        for fasta in sorted(self.protein_dir.iterdir()):
            if not fasta.is_file():
                continue

            if fasta.suffix.lower() not in self.VALID_SUFFIXES:
                continue

            protein_files.append(
                ProteinFile(
                    accession=fasta.stem,
                    fasta=fasta,
                )
            )

        return protein_files
