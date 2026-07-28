"""Tests for DIAMOND tabular-output parsing."""

from pathlib import Path

import pytest

from gutsporepredict.search.exceptions import SearchOutputError
from gutsporepredict.search.parser import parse_diamond_output


def test_parse_diamond_output(tmp_path: Path) -> None:
    output_file = tmp_path / "hits.tsv"
    output_file.write_text(
        "query_1\ttarget_A\t90.0\t90\t100\t120\t"
        "1\t90\t5\t94\t1e-30\t200.0\n"
        "query_1\ttarget_B\t75.0\t80\t100\t100\t"
        "1\t80\t1\t80\t1e-15\t150.0\n",
        encoding="utf-8",
    )

    hits = parse_diamond_output(output_file)

    assert len(hits) == 2

    first = hits[0]

    assert first.query_id == "query_1"
    assert first.target_id == "target_A"
    assert first.identity == pytest.approx(90.0)
    assert first.query_coverage == pytest.approx(90.0)
    assert first.target_coverage == pytest.approx(75.0)
    assert first.evalue == pytest.approx(1e-30)
    assert first.bitscore == pytest.approx(200.0)
    assert first.method == "diamond"


def test_parse_empty_output(tmp_path: Path) -> None:
    output_file = tmp_path / "empty.tsv"
    output_file.write_text("", encoding="utf-8")

    hits = parse_diamond_output(output_file)

    assert hits == []


def test_reject_invalid_column_count(tmp_path: Path) -> None:
    output_file = tmp_path / "invalid.tsv"
    output_file.write_text(
        "query_1\ttarget_A\t90.0\n",
        encoding="utf-8",
    )

    with pytest.raises(
        SearchOutputError,
        match="Expected 12 columns",
    ):
        parse_diamond_output(output_file)
