"""End-to-end GutSporePredict pipeline runner."""

from __future__ import annotations

import shlex
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

from gutsporepredict.resources import resource_root


class PipelineError(RuntimeError):
    """Raised when the end-to-end pipeline cannot complete."""


@dataclass(frozen=True)
class PipelineConfig:
    """Configuration for one GutSporePredict analysis run."""

    genome_dir: Path
    output_dir: Path
    project_root: Path
    threads: int = 1
    minimum_assessment: float = 0.5


class PipelineRunner:
    """Run genome prediction through lifecycle classification."""

    GENOME_SUFFIXES = {
        ".fa",
        ".fna",
        ".fasta",
        ".fas",
    }

    def __init__(self, config: PipelineConfig) -> None:
        self.config = config

        configured_root = config.project_root.resolve()
        packaged_root = resource_root()

        if (
            configured_root / "scripts" / "06_run_hmmsearch.py"
        ).is_file():
            self.project_root = configured_root
        else:
            self.project_root = packaged_root

        self.protein_dir = config.output_dir / "01_proteins"
        self.prodigal_log_dir = config.output_dir / "logs" / "prodigal"

        self.sporulation_search_dir = (
            config.output_dir / "02_hmmsearch_sporulation"
        )
        self.germination_search_dir = (
            config.output_dir / "03_hmmsearch_germination"
        )
        self.spoiiq_clostridia_search_dir = (
            config.output_dir / "02b_hmmsearch_spoIIQ_clostridia"
        )

        self.sporulation_calls_dir = (
            config.output_dir / "04_gene_calls_sporulation"
        )
        self.germination_calls_dir = (
            config.output_dir / "05_gene_calls_germination"
        )

        self.combined_calls = (
            config.output_dir / "06_gene_calls" / "gene_calls.tsv"
        )
        self.module_dir = config.output_dir / "07_modules"
        self.stage_dir = config.output_dir / "08_stages"
        self.lifecycle_file = (
            config.output_dir / "lifecycle_summary.tsv"
        )
        self.pipeline_log = config.output_dir / "logs" / "pipeline.log"

    def run(self) -> Path:
        """Run the complete analysis and return the summary path."""

        self._validate_configuration()
        self._prepare_directories()

        print("[1/9] Predicting proteins with Prodigal")
        self._run_prodigal()

        print("[2/9] Searching sporulation HMM profiles")
        self._run_hmmsearch(
            targets=self._root(
                "config/gtdb_targets/02_sporulation_ready.tsv"
            ),
            output_dir=self.sporulation_search_dir,
        )

        print("[3/9] Searching Clostridia SpoIIQ profile")
        self._run_hmmsearch(
            targets=self._root(
                "config/gtdb_targets/04_spoIIQ_clostridia.tsv"
            ),
            output_dir=self.spoiiq_clostridia_search_dir,
        )

        print("[4/9] Searching germination HMM profiles")
        self._run_hmmsearch(
            targets=self._root(
                "config/gtdb_targets/03_germination_ready.tsv"
            ),
            output_dir=self.germination_search_dir,
        )

        print("[5/9] Building lineage-aware three-state gene calls")
        self._run_tristate_calls()

        print("[6/9] Combining sporulation and germination calls")
        self._run_python_script(
            "scripts/12_merge_gene_call_matrices.py",
            "--input",
            str(self.sporulation_calls_dir / "gene_calls.tsv"),
            "--input",
            str(self.germination_calls_dir / "gene_calls.tsv"),
            "--output",
            str(self.combined_calls),
        )

        print("[7/9] Evaluating biological modules")
        self._run_python_script(
            "scripts/08_build_module_matrix.py",
            "--presence-matrix",
            str(self.combined_calls),
            "--knowledge-dir",
            str(self._root("knowledge")),
            "--output-dir",
            str(self.module_dir),
        )

        print("[8/9] Evaluating developmental stages ST001-ST009")
        self._run_python_script(
            "scripts/11_build_stage_matrix.py",
            "--module-evaluations",
            str(self.module_dir / "module_evaluations.tsv"),
            "--stage-definitions",
            str(self._root("knowledge/stages.yaml")),
            "--output-dir",
            str(self.stage_dir),
        )

        print("[9/9] Building lifecycle summary")
        self._run_python_script(
            "scripts/13_build_lifecycle_summary.py",
            "--stage-evaluations",
            str(self.stage_dir / "stage_evaluations.tsv"),
            "--output",
            str(self.lifecycle_file),
            "--minimum-assessment",
            str(self.config.minimum_assessment),
        )

        if not self.lifecycle_file.is_file():
            raise PipelineError(
                "Pipeline finished without creating lifecycle summary: "
                f"{self.lifecycle_file}"
            )

        print()
        print("GutSporePredict analysis completed.")
        print(f"Main result: {self.lifecycle_file}")

        return self.lifecycle_file

    def _root(self, relative_path: str) -> Path:
        return self.project_root / relative_path

    def _validate_configuration(self) -> None:
        genome_dir = self.config.genome_dir

        if not genome_dir.exists():
            raise PipelineError(
                f"Genome directory does not exist: {genome_dir}"
            )

        if not genome_dir.is_dir():
            raise PipelineError(
                f"Genome path is not a directory: {genome_dir}"
            )

        if self.config.threads < 1:
            raise PipelineError("threads must be at least 1")

        if not 0.0 <= self.config.minimum_assessment <= 1.0:
            raise PipelineError(
                "minimum_assessment must be between 0 and 1"
            )

        genome_files = self._genome_files()

        if not genome_files:
            raise PipelineError(
                "No genome FASTA files were found in "
                f"{genome_dir}. Supported suffixes: "
                + ", ".join(sorted(self.GENOME_SUFFIXES))
            )

        required_paths = [
            self._root("scripts/06_run_hmmsearch.py"),
            self._root(
                "scripts/10_build_groupwise_tristate_matrix.py"
            ),
            self._root(
                "config/gtdb_targets/04_spoIIQ_clostridia.tsv"
            ),
            self._root(
                "database/gutspore/reference_v3/hmm/"
                "spoIIQ_Clostridia.hmm"
            ),
            self._root(
                "scripts/15_merge_spoIIQ_lineage_evidence.py"
            ),
            self._root("scripts/12_merge_gene_call_matrices.py"),
            self._root("scripts/08_build_module_matrix.py"),
            self._root("scripts/11_build_stage_matrix.py"),
            self._root("scripts/13_build_lifecycle_summary.py"),
            self._root(
                "config/gtdb_targets/02_sporulation_ready.tsv"
            ),
            self._root(
                "config/gtdb_targets/03_germination_ready.tsv"
            ),
            self._root("config/hmm/competition_groups.tsv"),
            self._root(
                "config/hmm/germination_competition_groups.tsv"
            ),
            self._root("knowledge/stages.yaml"),
            self._root("database/gutspore/reference_v3/hmm"),
        ]

        missing_paths = [
            path
            for path in required_paths
            if not path.exists()
        ]

        if missing_paths:
            formatted = "\n".join(
                f"  - {path}"
                for path in missing_paths
            )
            raise PipelineError(
                "Required GutSporePredict resources are missing:\n"
                f"{formatted}\n"
                "Use --project-root to specify the repository root."
            )

        for executable in ("prodigal", "hmmsearch"):
            if shutil.which(executable) is None:
                raise PipelineError(
                    f"Required executable was not found: {executable}"
                )

    def _prepare_directories(self) -> None:
        self.config.output_dir.mkdir(
            parents=True,
            exist_ok=True,
        )
        self.protein_dir.mkdir(
            parents=True,
            exist_ok=True,
        )
        self.prodigal_log_dir.mkdir(
            parents=True,
            exist_ok=True,
        )
        self.pipeline_log.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        self.pipeline_log.write_text(
            "GutSporePredict pipeline commands\n",
            encoding="utf-8",
        )

    def _genome_files(self) -> list[Path]:
        return sorted(
            path
            for path in self.config.genome_dir.iterdir()
            if (
                path.is_file()
                and path.suffix.lower() in self.GENOME_SUFFIXES
            )
        )

    def _run_prodigal(self) -> None:
        for genome in self._genome_files():
            accession = genome.stem

            protein_fasta = self.protein_dir / f"{accession}.faa"
            nucleotide_fasta = self.protein_dir / f"{accession}.ffn"
            gff_file = self.protein_dir / f"{accession}.gff"
            stdout_file = (
                self.prodigal_log_dir / f"{accession}.stdout.log"
            )
            stderr_file = (
                self.prodigal_log_dir / f"{accession}.stderr.log"
            )

            command = [
                "prodigal",
                "-i",
                str(genome),
                "-a",
                str(protein_fasta),
                "-d",
                str(nucleotide_fasta),
                "-o",
                str(gff_file),
                "-f",
                "gff",
                "-p",
                "single",
            ]

            self._record_command(command)
            print(f"  Prodigal: {accession}")

            with (
                stdout_file.open("w", encoding="utf-8") as stdout,
                stderr_file.open("w", encoding="utf-8") as stderr,
            ):
                completed = subprocess.run(
                    command,
                    stdout=stdout,
                    stderr=stderr,
                    text=True,
                    check=False,
                )

            if completed.returncode != 0:
                raise PipelineError(
                    f"Prodigal failed for {genome}. "
                    f"See: {stderr_file}"
                )

            if not protein_fasta.is_file():
                raise PipelineError(
                    "Prodigal did not create the expected protein "
                    f"FASTA: {protein_fasta}"
                )

    def _run_hmmsearch(
        self,
        targets: Path,
        output_dir: Path,
    ) -> None:
        self._run_python_script(
            "scripts/06_run_hmmsearch.py",
            "--targets",
            str(targets),
            "--hmm-dir",
            str(
                self._root(
                    "database/gutspore/reference_v3/hmm"
                )
            ),
            "--protein-dir",
            str(self.protein_dir),
            "--output-dir",
            str(output_dir),
            "--cpu",
            str(self.config.threads),
            "--evalue",
            "1e-5",
            "--domain-evalue",
            "1e-4",
        )

    def _run_tristate_calls(self) -> None:
        self._run_python_script(
            "scripts/10_build_groupwise_tristate_matrix.py",
            "--hmmsearch-dir",
            str(self.sporulation_search_dir),
            "--targets",
            str(
                self._root(
                    "config/gtdb_targets/"
                    "02_sporulation_ready.tsv"
                )
            ),
            "--competition-groups",
            str(
                self._root(
                    "config/hmm/competition_groups.tsv"
                )
            ),
            "--output-dir",
            str(self.sporulation_calls_dir),
        )

        self._run_python_script(
            "scripts/15_merge_spoIIQ_lineage_evidence.py",
            "--input-calls",
            str(
                self.sporulation_calls_dir
                / "gene_calls.tsv"
            ),
            "--input-details",
            str(
                self.sporulation_calls_dir
                / "gene_call_details.tsv"
            ),
            "--clostridia-search-dir",
            str(self.spoiiq_clostridia_search_dir),
            "--sporulation-search-dir",
            str(self.sporulation_search_dir),
            "--output-calls",
            str(
                self.sporulation_calls_dir
                / "gene_calls.tsv"
            ),
            "--output-details",
            str(
                self.sporulation_calls_dir
                / "gene_call_details.tsv"
            ),
            "--audit-output",
            str(
                self.sporulation_calls_dir
                / "spoIIQ_lineage_evidence.tsv"
            ),
        )

        self._run_python_script(
            "scripts/10_build_groupwise_tristate_matrix.py",
            "--hmmsearch-dir",
            str(self.germination_search_dir),
            "--targets",
            str(
                self._root(
                    "config/gtdb_targets/"
                    "03_germination_ready.tsv"
                )
            ),
            "--competition-groups",
            str(
                self._root(
                    "config/hmm/"
                    "germination_competition_groups.tsv"
                )
            ),
            "--output-dir",
            str(self.germination_calls_dir),
        )

    def _run_python_script(
        self,
        relative_script: str,
        *arguments: str,
    ) -> None:
        script = self._root(relative_script)
        command = [
            sys.executable,
            str(script),
            *arguments,
        ]

        self._record_command(command)

        completed = subprocess.run(
            command,
            text=True,
            capture_output=True,
            check=False,
            cwd=self.project_root,
        )

        if completed.stdout:
            print(completed.stdout, end="")

        if completed.returncode != 0:
            if completed.stderr:
                print(completed.stderr, file=sys.stderr, end="")

            raise PipelineError(
                f"Pipeline command failed: "
                f"{shlex.join(command)}"
            )

    def _record_command(self, command: list[str]) -> None:
        with self.pipeline_log.open(
            "a",
            encoding="utf-8",
        ) as handle:
            handle.write(shlex.join(command))
            handle.write("\n")
