"""Tests for sequence-search data models."""

from pathlib import Path

from gutsporepredict.search.models import (
    SearchHit,
    SearchResult,
)


def make_hit(
    target_id: str,
    bitscore: float,
) -> SearchHit:
    return SearchHit(
        query_id="query_1",
        target_id=target_id,
        identity=90.0,
        alignment_length=90,
        query_length=100,
        target_length=100,
        query_coverage=90.0,
        target_coverage=90.0,
        evalue=1e-20,
        bitscore=bitscore,
        method="diamond",
    )


def test_search_result_summary() -> None:
    result = SearchResult(
        query_file=Path("query.faa"),
        database=Path("reference.dmnd"),
        output_file=Path("hits.tsv"),
        method="diamond",
        hits=[
            make_hit("target_A", 100.0),
            make_hit("target_B", 200.0),
        ],
    )

    assert result.hit_count == 2
    assert result.query_count == 1
    assert result.best_hit("query_1") is not None
    assert result.best_hit("query_1").target_id == "target_B"
    assert result.best_hit("missing") is None
