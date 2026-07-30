"""YAML loading utilities for the GutSporePredict knowledge base."""

from pathlib import Path
from typing import TypeAlias

import yaml

from gutsporepredict.knowledge.exceptions import (
    KnowledgeFileNotFoundError,
    KnowledgeFormatError,
    KnowledgeLoadError,
)

KnowledgeMapping: TypeAlias = dict[str, object]


def load_yaml(path: str | Path) -> KnowledgeMapping:
    """Load one YAML knowledge file as a string-keyed mapping.

    The loader is intentionally responsible only for reading and
    basic structural validation. Interpretation of individual fields
    belongs to the knowledge models, while biological consistency
    checks belong to the knowledge database.

    Empty YAML files are represented as empty dictionaries.

    Args:
        path:
            Path to the YAML file.

    Returns:
        A dictionary containing the top-level YAML mapping.

    Raises:
        KnowledgeFileNotFoundError:
            If the requested file does not exist.

        KnowledgeFormatError:
            If the YAML cannot be parsed, its top-level value is not a
            mapping, or one of its top-level keys is not a string.

        KnowledgeLoadError:
            If the file cannot be read for another operating-system
            reason.
    """

    file_path = Path(path)

    if not file_path.exists():
        raise KnowledgeFileNotFoundError(
            f"Knowledge file does not exist: {file_path}"
        )

    if not file_path.is_file():
        raise KnowledgeFormatError(
            file_path,
            "Expected a regular file.",
        )

    try:
        text = file_path.read_text(encoding="utf-8")
    except OSError as exc:
        raise KnowledgeLoadError(
            f"Could not read knowledge file {file_path}: {exc}"
        ) from exc

    try:
        data = yaml.safe_load(text)
    except yaml.YAMLError as exc:
        raise KnowledgeFormatError(
            file_path,
            f"Could not parse YAML: {exc}",
        ) from exc

    if data is None:
        return {}

    if not isinstance(data, dict):
        raise KnowledgeFormatError(
            file_path,
            "The top-level YAML value must be a mapping.",
        )

    if not all(isinstance(key, str) for key in data):
        raise KnowledgeFormatError(
            file_path,
            "All top-level YAML keys must be strings.",
        )

    return dict(data)
