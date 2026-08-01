"""Command-line interface for GutSporePredict."""

from __future__ import annotations

import argparse
from collections.abc import Sequence
from pathlib import Path

from gutsporepredict import __version__
from gutsporepredict.pipeline import (
    PipelineConfig,
    PipelineError,
    PipelineRunner,
)


def detect_project_root() -> Path:
    """Locate the GutSporePredict source repository."""

    candidate = Path(__file__).resolve().parents[2]

    required = [
        candidate / "scripts",
        candidate / "knowledge",
        candidate / "config",
        candidate / "database",
    ]

    if all(path.exists() for path in required):
        return candidate

    return Path.cwd()


def build_parser() -> argparse.ArgumentParser:
    """Create the top-level command-line parser."""

    parser = argparse.ArgumentParser(
        prog="gutsporepredict",
        description=(
            "Comparative genomics and evolutionary analysis of "
            "sporulation and germination potential."
        ),
    )

    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    subparsers = parser.add_subparsers(
        dest="command",
        title="commands",
    )

    subparsers.add_parser(
        "doctor",
        help="Check the GutSporePredict execution environment.",
    )

    run_parser = subparsers.add_parser(
        "run",
        help=(
            "Run protein prediction, HMM searches and lifecycle "
            "classification."
        ),
    )

    run_parser.add_argument(
        "--genomes",
        type=Path,
        required=True,
        help="Directory containing genome FASTA files.",
    )
    run_parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Output directory for the complete analysis.",
    )
    run_parser.add_argument(
        "--project-root",
        type=Path,
        default=detect_project_root(),
        help=(
            "GutSporePredict repository root containing scripts, "
            "knowledge, config and database directories. "
            "Default: automatically detected installation source."
        ),
    )
    run_parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="CPU threads used by each HMMER search. Default: 1.",
    )
    run_parser.add_argument(
        "--minimum-assessment",
        type=float,
        default=0.5,
        help=(
            "Minimum assessed fraction required for a definitive "
            "lifecycle call. Default: 0.5."
        ),
    )

    return parser


def run_doctor() -> int:
    """Run a minimal environment check."""

    print("GutSporePredict environment")
    print(f"Version: {__version__}")
    print("Status: package installation successful")
    return 0


def run_analysis(
    genome_dir: Path,
    output_dir: Path,
    project_root: Path,
    threads: int,
    minimum_assessment: float,
) -> int:
    """Run the complete GutSporePredict pipeline."""

    config = PipelineConfig(
        genome_dir=genome_dir.resolve(),
        output_dir=output_dir.resolve(),
        project_root=project_root.resolve(),
        threads=threads,
        minimum_assessment=minimum_assessment,
    )

    PipelineRunner(config).run()
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    """Run the GutSporePredict command-line interface."""

    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == "doctor":
        return run_doctor()

    if args.command == "run":
        try:
            return run_analysis(
                genome_dir=args.genomes,
                output_dir=args.output,
                project_root=args.project_root,
                threads=args.threads,
                minimum_assessment=args.minimum_assessment,
            )
        except (
            FileNotFoundError,
            NotADirectoryError,
            PipelineError,
            ValueError,
        ) as error:
            print(f"ERROR: {error}")
            return 1

    parser.print_help()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
