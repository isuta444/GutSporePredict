"""Bundled runtime resources for GutSporePredict."""

from pathlib import Path


def resource_root() -> Path:
    """Return the root directory containing bundled runtime resources."""
    return Path(__file__).resolve().parent


__all__ = ["resource_root"]
