"""Protein FASTA input validation."""

from pathlib import Path

from gutsporepredict.exceptions import InputValidationError


def validate_protein_directory(
    protein_dir: str | Path,
) -> Path:
    """Validate a directory containing protein FASTA files."""

    protein_dir = Path(protein_dir)

    if not protein_dir.exists():
        raise InputValidationError(
            f"Protein input directory does not exist: {protein_dir}"
        )

    if not protein_dir.is_dir():
        raise InputValidationError(
            f"Protein input path is not a directory: {protein_dir}"
        )

    return protein_dir


def validate_protein_fasta_files(
    protein_dir: str | Path,
) -> list[Path]:
    """Validate that the directory contains FAA files."""

    protein_dir = Path(protein_dir)

    fasta_files = sorted(
        path
        for path in protein_dir.iterdir()
        if path.is_file() and path.suffix.lower() == ".faa"
    )

    if not fasta_files:
        raise InputValidationError(
            f"No protein FASTA files (.faa) were found in: {protein_dir}"
        )

    empty_files = [
        path for path in fasta_files if path.stat().st_size == 0
    ]

    if empty_files:
        empty_names = ", ".join(path.name for path in empty_files)

        raise InputValidationError(
            f"Empty protein FASTA files were found: {empty_names}"
        )

    return fasta_files
