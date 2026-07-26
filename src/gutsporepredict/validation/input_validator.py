"""Input validation utilities."""

from pathlib import Path


class InputValidator:
    """Validate user input before running the pipeline."""

    @staticmethod
    def validate_genome_directory(genome_dir: Path) -> None:
        """Validate that the genome directory exists."""

        if not genome_dir.exists():
            raise FileNotFoundError(
                f"Genome directory not found: {genome_dir}"
            )

        if not genome_dir.is_dir():
            raise NotADirectoryError(
                f"Not a directory: {genome_dir}"
            )

    @staticmethod
    def validate_fasta_files(genome_dir: Path) -> None:
        """Validate that the directory contains FASTA files."""

        fasta_extensions = {".fa", ".fasta", ".fna"}

        fasta_files = [
            path
            for path in genome_dir.iterdir()
            if path.is_file() and path.suffix.lower() in fasta_extensions
        ]

        if not fasta_files:
            raise FileNotFoundError(
                f"No FASTA files found in directory: {genome_dir}"
            )