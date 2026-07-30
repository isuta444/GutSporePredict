"""Tests for the YAML knowledge loader."""

from pathlib import Path

import pytest

from gutsporepredict.knowledge.exceptions import (
    KnowledgeFileNotFoundError,
    KnowledgeFormatError,
)
from gutsporepredict.knowledge.loader import load_yaml


def write_yaml(
    directory: Path,
    filename: str,
    content: str,
) -> Path:
    """Write a temporary YAML file and return its path."""

    path = directory / filename
    path.write_text(content, encoding="utf-8")
    return path


def test_load_yaml_returns_top_level_mapping(
    tmp_path: Path,
) -> None:
    """A valid YAML mapping should be returned as a dictionary."""

    path = write_yaml(
        tmp_path,
        "modules.yaml",
        """
knowledge_version: "1.0"
modules:
  - module_id: sporulation_initiation
    name: Sporulation initiation
""",
    )

    result = load_yaml(path)

    assert result["knowledge_version"] == "1.0"

    modules = result["modules"]

    assert isinstance(modules, list)
    assert modules[0]["module_id"] == "sporulation_initiation"


def test_load_yaml_returns_empty_mapping_for_empty_file(
    tmp_path: Path,
) -> None:
    """An empty YAML file should produce an empty dictionary."""

    path = write_yaml(
        tmp_path,
        "empty.yaml",
        "",
    )

    assert load_yaml(path) == {}


def test_load_yaml_raises_for_missing_file(
    tmp_path: Path,
) -> None:
    """A missing knowledge file should raise a dedicated error."""

    path = tmp_path / "missing.yaml"

    with pytest.raises(
        KnowledgeFileNotFoundError,
        match="does not exist",
    ):
        load_yaml(path)


def test_load_yaml_rejects_invalid_yaml(
    tmp_path: Path,
) -> None:
    """Malformed YAML should raise a knowledge-format error."""

    path = write_yaml(
        tmp_path,
        "invalid.yaml",
        """
modules:
  - module_id: sporulation
    genes: [spo0A, spo0F
""",
    )

    with pytest.raises(
        KnowledgeFormatError,
        match="Could not parse YAML",
    ):
        load_yaml(path)


def test_load_yaml_rejects_top_level_sequence(
    tmp_path: Path,
) -> None:
    """The top-level YAML value must be a mapping."""

    path = write_yaml(
        tmp_path,
        "sequence.yaml",
        """
- module_id: sporulation_initiation
- module_id: engulfment
""",
    )

    with pytest.raises(
        KnowledgeFormatError,
        match="top-level YAML value must be a mapping",
    ):
        load_yaml(path)


def test_load_yaml_rejects_non_string_top_level_keys(
    tmp_path: Path,
) -> None:
    """Top-level knowledge keys must be strings."""

    path = write_yaml(
        tmp_path,
        "numeric_key.yaml",
        """
1:
  name: invalid
""",
    )

    with pytest.raises(
        KnowledgeFormatError,
        match="top-level YAML keys must be strings",
    ):
        load_yaml(path)
