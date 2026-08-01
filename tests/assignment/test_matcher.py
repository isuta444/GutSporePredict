"""Tests for reference-gene matching."""

from gutsporepredict.assignment.matcher import GeneMatcher
from gutsporepredict.reference.loader import ReferenceLoader
from gutsporepredict.search.models import SearchHit


def make_hit(target_id: str) -> SearchHit:
    return SearchHit(
        query_id="query_001",
        target_id=target_id,
        identity=80.0,
        alignment_length=200,
        query_length=220,
        target_length=210,
        query_coverage=90.9,
        target_coverage=95.2,
        evalue=1e-50,
        bitscore=300.0,
        method="diamond",
    )


def load_matcher() -> GeneMatcher:
    database = ReferenceLoader().load(
        "database/reference/genes.tsv",
        "database/reference/aliases.tsv",
    )

    return GeneMatcher(database)


def test_match_gene_id() -> None:
    gene = load_matcher().match(make_hit("GSP0001"))

    assert gene is not None
    assert gene.canonical_name == "spo0A"


def test_match_canonical_name() -> None:
    gene = load_matcher().match(make_hit("spoIIAC"))

    assert gene is not None
    assert gene.gene_id == "GSP0004"


def test_match_alias() -> None:
    gene = load_matcher().match(make_hit("sigF"))

    assert gene is not None
    assert gene.gene_id == "GSP0004"


def test_match_pipe_delimited_target() -> None:
    gene = load_matcher().match(
        make_hit("reference|GSP0001|spo0A")
    )

    assert gene is not None
    assert gene.gene_id == "GSP0001"


def test_unknown_target_returns_none() -> None:
    gene = load_matcher().match(
        make_hit("unknown_target")
    )

    assert gene is None
