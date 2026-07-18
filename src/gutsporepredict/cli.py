"""Command-line interface for GutSporePredict."""

from __future__ import annotations

import argparse
from collections.abc import Sequence

from gutsporepredict import __version__


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

    return parser


def run_doctor() -> int:
    """Run a minimal environment check."""
    print("GutSporePredict environment")
    print(f"Version: {__version__}")
    print("Status: package installation successful")
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    """Run the GutSporePredict command-line interface."""
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == "doctor":
        return run_doctor()

    parser.print_help()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
