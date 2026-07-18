"""Logging configuration for GutSporePredict."""

from __future__ import annotations

import logging
from pathlib import Path

DEFAULT_LOG_FORMAT = (
    "%(asctime)s | %(levelname)s | %(name)s | %(message)s"
)


def configure_logging(
    *,
    level: int = logging.INFO,
    log_file: str | Path | None = None,
) -> None:
    """Configure console and optional file logging.

    This function clears existing root handlers so that repeated calls
    from tests or command-line entry points remain predictable.
    """
    root_logger = logging.getLogger()
    root_logger.setLevel(level)
    root_logger.handlers.clear()

    formatter = logging.Formatter(DEFAULT_LOG_FORMAT)

    console_handler = logging.StreamHandler()
    console_handler.setLevel(level)
    console_handler.setFormatter(formatter)
    root_logger.addHandler(console_handler)

    if log_file is not None:
        path = Path(log_file).expanduser().resolve()
        path.parent.mkdir(parents=True, exist_ok=True)

        file_handler = logging.FileHandler(
            path,
            encoding="utf-8",
        )
        file_handler.setLevel(level)
        file_handler.setFormatter(formatter)
        root_logger.addHandler(file_handler)
