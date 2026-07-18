"""Tests for common file-system utilities."""

from pathlib import Path

import pytest

from gutsporepredict.exceptions import InputValidationError
from gutsporepredict.io.files import ensure_directory, require_file


def test_require_file_accepts_nonempty_file(
    tmp_path: Path,
) -> None:
    """A nonempty regular file should pass validation."""
    test_file = tmp_path / "genome.fna"
    test_file.write_text(">contig1\nATGC\n", encoding="utf-8")

    result = require_file(test_file)

    assert result == test_file.resolve()


def test_require_file_rejects_missing_file(
    tmp_path: Path,
) -> None:
    """A missing file should raise an expected validation error."""
    missing = tmp_path / "missing.fna"

    with pytest.raises(
        InputValidationError,
        match="does not exist",
    ):
        require_file(missing)


def test_require_file_rejects_empty_file(
    tmp_path: Path,
) -> None:
    """An empty file should be rejected by default."""
    empty = tmp_path / "empty.fna"
    empty.touch()

    with pytest.raises(
        InputValidationError,
        match="empty",
    ):
        require_file(empty)


def test_ensure_directory_creates_nested_directory(
    tmp_path: Path,
) -> None:
    """Nested output directories should be created automatically."""
    output = tmp_path / "results" / "gene_prediction"

    result = ensure_directory(output)

    assert result.is_dir()
    assert result == output.resolve()
