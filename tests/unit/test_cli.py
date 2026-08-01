"""Tests for the GutSporePredict command-line interface."""

import pytest

from gutsporepredict.cli import main


def test_doctor_command(
    capsys: pytest.CaptureFixture[str],
) -> None:
    """The doctor command should report an installed package."""
    exit_code = main(["doctor"])
    captured = capsys.readouterr()

    assert exit_code == 0
    assert "GutSporePredict environment" in captured.out
    assert "4.0.0b1" in captured.out


def test_no_command_displays_help(
    capsys: pytest.CaptureFixture[str],
) -> None:
    """Calling the CLI without a subcommand should display help."""
    exit_code = main([])
    captured = capsys.readouterr()

    assert exit_code == 0
    assert "usage:" in captured.out


def test_run_command_accepts_public_options() -> None:
    """The run parser should expose the public pipeline options."""

    from gutsporepredict.cli import build_parser

    parser = build_parser()
    args = parser.parse_args(
        [
            "run",
            "--genomes",
            "input/genomes",
            "--output",
            "results/run1",
            "--threads",
            "2",
            "--minimum-assessment",
            "0.6",
        ]
    )

    assert args.command == "run"
    assert str(args.genomes) == "input/genomes"
    assert str(args.output) == "results/run1"
    assert args.threads == 2
    assert args.minimum_assessment == 0.6
