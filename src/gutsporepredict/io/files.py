"""Common file-system utilities."""

from __future__ import annotations

from pathlib import Path

from gutsporepredict.exceptions import InputValidationError


def require_file(
    path: str | Path,
    *,
    allow_empty: bool = False,
) -> Path:
    """Validate that a required file exists.

    Parameters
    ----------
    path
        File path to validate.
    allow_empty
        Permit a zero-byte file when True.

    Returns
    -------
    pathlib.Path
        Resolved validated file path.

    Raises
    ------
    InputValidationError
        If the path does not exist, is not a file, or is empty.
    """
    resolved = Path(path).expanduser().resolve()

    if not resolved.exists():
        raise InputValidationError(
            f"Required file does not exist: {resolved}"
        )

    if not resolved.is_file():
        raise InputValidationError(
            f"Expected a file but found another path type: {resolved}"
        )

    if not allow_empty and resolved.stat().st_size == 0:
        raise InputValidationError(
            f"Required file is empty: {resolved}"
        )

    return resolved


def ensure_directory(path: str | Path) -> Path:
    """Create a directory if necessary and return its resolved path."""
    resolved = Path(path).expanduser().resolve()
    resolved.mkdir(parents=True, exist_ok=True)
    return resolved
