"""Tests for gene-assignment confidence rules."""

from gutsporepredict.assignment.models import AssignmentConfidence
from gutsporepredict.assignment.rules import AssignmentRuleSet
from gutsporepredict.reference.loader import ReferenceLoader
from gutsporepredict.search.models import SearchHit


def make_hit(
    *,
    identity: float,
    query_coverage: float,
    target_coverage: float,
    evalue: float,
) -> SearchHit:
    return SearchHit(
        query_id="query_001",
        target_id="GSP0001",
        identity=identity,
        alignment_length=180,
        query_length=200,
        target_length=200,
        query_coverage=query_coverage,
        target_coverage=target_coverage,
        evalue=evalue,
        bitscore=250.0,
        method="diamond",
    )


def reference_gene():
    database = ReferenceLoader().load(
        "database/reference/genes.tsv",
        "database/reference/aliases.tsv",
    )

    gene = database.gene_by_id("GSP0001")
    assert gene is not None
    return gene


def test_high_confidence_assignment() -> None:
    result = AssignmentRuleSet().evaluate(
        make_hit(
            identity=70.0,
            query_coverage=90.0,
            target_coverage=90.0,
            evalue=1e-60,
        ),
        reference_gene(),
    )

    assert result.accepted
    assert result.confidence is AssignmentConfidence.HIGH
    assert result.failed_filters == ()


def test_medium_confidence_assignment() -> None:
    result = AssignmentRuleSet().evaluate(
        make_hit(
            identity=42.0,
            query_coverage=70.0,
            target_coverage=70.0,
            evalue=1e-15,
        ),
        reference_gene(),
    )

    assert result.accepted
    assert result.confidence is AssignmentConfidence.MEDIUM


def test_low_confidence_assignment() -> None:
    result = AssignmentRuleSet().evaluate(
        make_hit(
            identity=28.0,
            query_coverage=50.0,
            target_coverage=50.0,
            evalue=1e-6,
        ),
        reference_gene(),
    )

    assert result.accepted
    assert result.confidence is AssignmentConfidence.LOW


def test_rejected_assignment() -> None:
    result = AssignmentRuleSet().evaluate(
        make_hit(
            identity=20.0,
            query_coverage=30.0,
            target_coverage=30.0,
            evalue=0.01,
        ),
        reference_gene(),
    )

    assert not result.accepted
    assert result.confidence is AssignmentConfidence.REJECTED
    assert "identity" in result.failed_filters
