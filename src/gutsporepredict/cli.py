"""Command-line interface for GutSporePredict."""

from __future__ import annotations

import argparse
from collections.abc import Sequence
from pathlib import Path

from gutsporepredict import __version__
from gutsporepredict.io import GenomeLoader


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
        help="Load genome FASTA files and start the analysis pipeline.",
    )

    run_parser.add_argument(
        "--genomes",
        type=Path,
        required=True,
        help="Directory containing genome FASTA files.",
    )

    return parser


def run_doctor() -> int:
    """Run a minimal environment check."""
    print("GutSporePredict environment")
    print(f"Version: {__version__}")
    print("Status: package installation successful")
    return 0


def run_pipeline(genome_dir: Path) -> int:
    """Load genomes and run the initial pipeline step."""
    print(f"Loading genomes from: {genome_dir}")

    loader = GenomeLoader(genome_dir)
    genomes = loader.load()

    print(f"Loaded {len(genomes)} genome(s)")

    for genome in genomes:
        print(f"- {genome.accession}")

    print("Done.")
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    """Run the GutSporePredict command-line interface."""
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == "doctor":
        return run_doctor()

    if args.command == "run":
        return run_pipeline(args.genomes)

    parser.print_help()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())