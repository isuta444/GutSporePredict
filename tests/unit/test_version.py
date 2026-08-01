"""Tests for package version information."""

from gutsporepredict import __version__


def test_version() -> None:
    """The Version 4 beta identifier should be exposed."""
    assert __version__ == "4.0.0b1"
