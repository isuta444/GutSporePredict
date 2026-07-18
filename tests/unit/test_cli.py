"""Tests for the GutSporePredict command-line interface."""

from gutsporepredict.cli import main


def test_doctor_command(
    capsys,
) -> None:
    """The doctor command should report an installed package."""
    exit_code = main(["doctor"])
    captured = capsys.readouterr()

    assert exit_code == 0
    assert "GutSporePredict environment" in captured.out
    assert "4.0.0a1" in captured.out


def test_no_command_displays_help(
    capsys,
) -> None:
    """Calling the CLI without a subcommand should display help."""
    exit_code = main([])
    captured = capsys.readouterr()

    assert exit_code == 0
    assert "usage:" in captured.out
