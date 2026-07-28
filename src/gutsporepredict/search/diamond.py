"""DIAMOND protein-search engine."""

import shlex
import shutil
import subprocess
from pathlib import Path

from gutsporepredict.search.base import SearchEngine
from gutsporepredict.search.exceptions import (
    SearchExecutableNotFoundError,
    SearchExecutionError,
    SearchOutputError,
)
from gutsporepredict.search.models import SearchResult
from gutsporepredict.search.parser import parse_diamond_output

DIAMOND_OUTPUT_FIELDS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "qlen",
    "slen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
]


class DiamondSearchEngine(SearchEngine):
    """Run DIAMOND protein searches."""

    def __init__(
        self,
        executable: str = "diamond",
        threads: int = 1,
        evalue: float = 1e-5,
        max_target_seqs: int = 25,
        sensitivity: str | None = "sensitive",
        extra_args: list[str] | None = None,
    ) -> None:
        if threads < 1:
            raise ValueError("threads must be at least 1")

        if evalue <= 0:
            raise ValueError("evalue must be greater than 0")

        if max_target_seqs < 1:
            raise ValueError(
                "max_target_seqs must be at least 1"
            )

        self.executable = executable
        self.threads = threads
        self.evalue = evalue
        self.max_target_seqs = max_target_seqs
        self.sensitivity = sensitivity
        self.extra_args = list(extra_args or [])

    def _resolve_executable(self) -> str:
        executable_path = shutil.which(self.executable)

        if executable_path is None:
            raise SearchExecutableNotFoundError(
                f"DIAMOND executable was not found: "
                f"{self.executable}"
            )

        return executable_path

    @staticmethod
    def _validate_fasta(path: str | Path, label: str) -> Path:
        fasta = Path(path)

        if not fasta.exists():
            raise SearchOutputError(
                f"{label} FASTA file does not exist: {fasta}"
            )

        if not fasta.is_file():
            raise SearchOutputError(
                f"{label} FASTA path is not a file: {fasta}"
            )

        if fasta.stat().st_size == 0:
            raise SearchOutputError(
                f"{label} FASTA file is empty: {fasta}"
            )

        return fasta

    @staticmethod
    def _database_file(database: str | Path) -> Path:
        database = Path(database)

        if database.suffix == ".dmnd":
            return database

        return Path(f"{database}.dmnd")

    @staticmethod
    def _database_prefix(database: str | Path) -> Path:
        database = Path(database)

        if database.suffix == ".dmnd":
            return database.with_suffix("")

        return database

    def make_database(
        self,
        reference_fasta: str | Path,
        database: str | Path,
        log_file: str | Path | None = None,
    ) -> Path:
        """Create a DIAMOND database from reference proteins."""

        executable = self._resolve_executable()
        reference_fasta = self._validate_fasta(
            reference_fasta,
            "Reference",
        )
        database_prefix = self._database_prefix(database)
        database_file = self._database_file(database_prefix)

        database_prefix.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        if log_file is None:
            log_path = database_prefix.parent / (
                f"{database_prefix.name}.makedb.log"
            )
        else:
            log_path = Path(log_file)
            log_path.parent.mkdir(parents=True, exist_ok=True)

        command = [
            executable,
            "makedb",
            "--in",
            str(reference_fasta),
            "--db",
            str(database_prefix),
        ]

        completed = self._run_command(command, log_path)

        if completed.returncode != 0:
            raise SearchExecutionError(
                "DIAMOND makedb failed with exit code "
                f"{completed.returncode}. See: {log_path}"
            )

        if not database_file.exists():
            raise SearchOutputError(
                "DIAMOND makedb completed, but the database "
                f"was not created: {database_file}"
            )

        return database_file

    def search(
        self,
        query_fasta: str | Path,
        database: str | Path,
        output_file: str | Path,
    ) -> SearchResult:
        """Run DIAMOND blastp and parse its tabular output."""

        executable = self._resolve_executable()
        query_fasta = self._validate_fasta(
            query_fasta,
            "Query",
        )

        database_file = self._database_file(database)
        database_prefix = self._database_prefix(database)

        if not database_file.exists():
            raise SearchOutputError(
                f"DIAMOND database does not exist: {database_file}"
            )

        output_file = Path(output_file)
        output_file.parent.mkdir(parents=True, exist_ok=True)

        log_file = output_file.with_suffix(
            output_file.suffix + ".log"
        )

        command = [
            executable,
            "blastp",
            "--query",
            str(query_fasta),
            "--db",
            str(database_prefix),
            "--out",
            str(output_file),
            "--outfmt",
            "6",
            *DIAMOND_OUTPUT_FIELDS,
            "--evalue",
            str(self.evalue),
            "--max-target-seqs",
            str(self.max_target_seqs),
            "--threads",
            str(self.threads),
        ]

        if self.sensitivity:
            command.append(f"--{self.sensitivity}")

        command.extend(self.extra_args)

        completed = self._run_command(command, log_file)

        if completed.returncode != 0:
            raise SearchExecutionError(
                "DIAMOND blastp failed with exit code "
                f"{completed.returncode}. See: {log_file}"
            )

        if not output_file.exists():
            raise SearchOutputError(
                "DIAMOND completed, but the output file "
                f"was not created: {output_file}"
            )

        hits = parse_diamond_output(output_file)

        return SearchResult(
            query_file=query_fasta,
            database=database_file,
            output_file=output_file,
            method="diamond",
            hits=hits,
        )

    @staticmethod
    def _run_command(
        command: list[str],
        log_file: Path,
    ) -> subprocess.CompletedProcess[str]:
        log_file.parent.mkdir(parents=True, exist_ok=True)

        try:
            completed = subprocess.run(
                command,
                capture_output=True,
                text=True,
                check=False,
            )
        except OSError as exc:
            raise SearchExecutionError(
                f"Failed to start command: {exc}"
            ) from exc

        command_text = shlex.join(command)

        log_file.write_text(
            "COMMAND\n"
            f"{command_text}\n\n"
            "STDOUT\n"
            f"{completed.stdout}\n\n"
            "STDERR\n"
            f"{completed.stderr}\n",
            encoding="utf-8",
        )

        return completed
