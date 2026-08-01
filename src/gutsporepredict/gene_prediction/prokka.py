"""Prokka gene predictor."""

import shutil
import subprocess
from pathlib import Path

from gutsporepredict.gene_prediction.base import (
    GenePredictionResult,
    GenePredictor,
)


class ProkkaError(RuntimeError):
    """Raised when Prokka gene prediction fails."""


class ProkkaPredictor(GenePredictor):
    """Run Prokka for prokaryotic gene prediction."""

    def __init__(
        self,
        executable: str = "prokka",
        extra_args: list[str] | None = None,
    ) -> None:
        self.executable = executable
        self.extra_args = list(extra_args or [])

    def predict(
        self,
        genome_fasta: str | Path,
        output_dir: str | Path,
        prefix: str | None = None,
    ) -> GenePredictionResult:
        """Run Prokka and return the generated output paths."""

        genome_fasta = Path(genome_fasta)
        output_dir = Path(output_dir)

        if not genome_fasta.exists():
            raise ProkkaError(
                f"Genome FASTA file does not exist: {genome_fasta}"
            )

        if not genome_fasta.is_file():
            raise ProkkaError(
                f"Genome FASTA path is not a file: {genome_fasta}"
            )

        executable_path = shutil.which(self.executable)

        if executable_path is None:
            raise ProkkaError(
                f"Prokka executable was not found: {self.executable}. "
                "Install Prokka or activate the appropriate environment."
            )

        genome_prefix = prefix or genome_fasta.stem

        output_dir.mkdir(parents=True, exist_ok=True)

        protein_fasta = output_dir / f"{genome_prefix}.faa"
        nucleotide_fasta = output_dir / f"{genome_prefix}.ffn"
        gff_file = output_dir / f"{genome_prefix}.gff"
        log_file = output_dir / "prokka.log"

        command = [
            executable_path,
            str(genome_fasta),
            "--outdir",
            str(output_dir),
            "--prefix",
            genome_prefix,
            "--force",
            *self.extra_args,
        ]

        try:
            completed_process = subprocess.run(
                command,
                capture_output=True,
                text=True,
                check=False,
            )
        except OSError as exc:
            raise ProkkaError(
                f"Failed to start Prokka: {exc}"
            ) from exc

        log_file.write_text(
            "COMMAND\n"
            f"{' '.join(command)}\n\n"
            "STDOUT\n"
            f"{completed_process.stdout}\n\n"
            "STDERR\n"
            f"{completed_process.stderr}\n",
            encoding="utf-8",
        )

        if completed_process.returncode != 0:
            raise ProkkaError(
                "Prokka failed with exit code "
                f"{completed_process.returncode}. "
                f"See log file: {log_file}"
            )

        required_outputs = [
            protein_fasta,
            nucleotide_fasta,
            gff_file,
        ]

        missing_outputs = [
            path for path in required_outputs if not path.exists()
        ]

        if missing_outputs:
            missing_text = ", ".join(
                str(path) for path in missing_outputs
            )
            raise ProkkaError(
                "Prokka finished, but required output files are missing: "
                f"{missing_text}. See log file: {log_file}"
            )

        return GenePredictionResult(
            input_fasta=genome_fasta,
            output_dir=output_dir,
            protein_fasta=protein_fasta,
            nucleotide_fasta=nucleotide_fasta,
            gff_file=gff_file,
            log_file=log_file,
        )
