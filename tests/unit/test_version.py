"""Tests for package version information."""

from gutsporepredict import __version__


def test_version() -> None:
    """The initial Version 4 alpha identifier should be exposed."""
    assert __version__ == "4.0.0a1"
