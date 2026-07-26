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