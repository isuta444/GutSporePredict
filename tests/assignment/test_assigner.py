"""Tests for the gene-assignment engine."""

import pytest

from gutsporepredict.assignment.assigner import GeneAssigner
from gutsporepredict.assignment.exceptions import (
    ReferenceGeneNotFoundError,
)
from gutsporepredict.assignment.models import (
    AssignmentConfidence,
    AssignmentMethod,
)
from gutsporepredict.reference.loader import ReferenceLoader
from gutsporepredict.search.models import SearchHit


def make_hit(
    *,
    query_id: str = "query_001",
    target_id: str = "GSP0001",
    identity: float = 70.0,
    query_coverage: float = 90.0,
    target_coverage: float = 90.0,
    evalue: float = 1e-50,
    bitscore: float = 300.0,
    method: str = "diamond",
) -> SearchHit:
    return SearchHit(
        query_id=query_id,
        target_id=target_id,
        identity=identity,
        alignment_length=180,
        query_length=200,
        target_length=200,
        query_coverage=query_coverage,
        target_coverage=target_coverage,
        evalue=evalue,
        bitscore=bitscore,
        method=method,
    )


def load_assigner() -> GeneAssigner:
    database = ReferenceLoader().load(
        "database/reference/genes.tsv",
        "database/reference/aliases.tsv",
    )

    return GeneAssigner(database)


def test_assign_search_hit() -> None:
    assignment = load_assigner().assign(make_hit())

    assert assignment is not None
    assert assignment.gene_id == "GSP0001"
    assert assignment.canonical_name == "spo0A"
    assert assignment.confidence is AssignmentConfidence.HIGH
    assert assignment.assignment_method is AssignmentMethod.DIAMOND
    assert assignment.accepted


def test_unresolved_target_returns_none() -> None:
    assignment = load_assigner().assign(
        make_hit(target_id="unknown")
    )

    assert assignment is None


def test_strict_matching_rejects_unknown_target() -> None:
    with pytest.raises(
        ReferenceGeneNotFoundError,
        match="could not be resolved",
    ):
        load_assigner().assign(
            make_hit(target_id="unknown"),
            strict_reference_matching=True,
        )


def test_assign_all_excludes_rejected_by_default() -> None:
    assignments = load_assigner().assign_all(
        [
            make_hit(),
            make_hit(
                query_id="query_002",
                target_id="GSP0004",
                identity=10.0,
                query_coverage=20.0,
                target_coverage=20.0,
                evalue=1.0,
            ),
        ]
    )

    assert len(assignments) == 1
    assert assignments[0].query_id == "query_001"


def test_assign_all_can_include_rejected() -> None:
    assignments = load_assigner().assign_all(
        [
            make_hit(
                identity=10.0,
                query_coverage=20.0,
                target_coverage=20.0,
                evalue=1.0,
            )
        ],
        include_rejected=True,
    )

    assert len(assignments) == 1
    assert (
        assignments[0].confidence
        is AssignmentConfidence.REJECTED
    )


def test_best_assignment_is_selected_per_query() -> None:
    assignments = load_assigner().best_assignments(
        [
            make_hit(
                target_id="GSP0001",
                identity=40.0,
                query_coverage=70.0,
                target_coverage=70.0,
                evalue=1e-15,
                bitscore=180.0,
            ),
            make_hit(
                target_id="GSP0004",
                identity=75.0,
                query_coverage=95.0,
                target_coverage=95.0,
                evalue=1e-60,
                bitscore=350.0,
            ),
        ]
    )

    assert len(assignments) == 1
    assert assignments[0].gene_id == "GSP0004"
    assert (
        assignments[0].confidence
        is AssignmentConfidence.HIGH
    )
