"""Tests for gene-evidence construction."""

from gutsporepredict.assignment.models import (
    AssignmentConfidence,
    AssignmentMethod,
    GeneAssignment,
)
from gutsporepredict.evidence.builder import EvidenceBuilder
from gutsporepredict.evidence.models import (
    EvidenceReason,
    EvidenceStatus,
)
from gutsporepredict.reference.models import ReferenceGene
from gutsporepredict.search.models import SearchHit


def make_reference_gene(
    gene_id: str = "GSP0001",
    canonical_name: str = "spo0A",
) -> ReferenceGene:
    """Create a reference gene for testing."""

    return ReferenceGene(
        gene_id=gene_id,
        canonical_name=canonical_name,
        pathway="sporulation",
        module="initiation",
        stage="stage_0",
        essentiality="core",
        phyletic_pattern="broad",
        search_methods=("diamond",),
        description="Test reference gene.",
    )


def make_assignment(
    reference_gene: ReferenceGene,
    confidence: AssignmentConfidence,
    *,
    bitscore: float = 200.0,
    query_id: str = "query_1",
) -> GeneAssignment:
    """Create a gene assignment for testing."""

    search_hit = SearchHit(
        query_id=query_id,
        target_id=reference_gene.gene_id,
        identity=60.0,
        alignment_length=200,
        query_length=220,
        target_length=210,
        query_coverage=90.0,
        target_coverage=95.0,
        evalue=1e-40,
        bitscore=bitscore,
        method="diamond",
    )

    failed_filters: tuple[str, ...] = ()

    if confidence is AssignmentConfidence.REJECTED:
        failed_filters = ("identity",)

    return GeneAssignment(
        query_id=query_id,
        reference_gene=reference_gene,
        search_hit=search_hit,
        confidence=confidence,
        assignment_method=AssignmentMethod.DIAMOND,
        passed_filters=("query_coverage",),
        failed_filters=failed_filters,
    )


def test_build_present_evidence() -> None:
    """An accepted assignment should produce present evidence."""

    reference_gene = make_reference_gene()
    assignment = make_assignment(
        reference_gene,
        AssignmentConfidence.HIGH,
    )

    evidence = EvidenceBuilder().build(
        reference_gene,
        [assignment],
    )

    assert evidence.status is EvidenceStatus.PRESENT
    assert evidence.reason is EvidenceReason.ACCEPTED_ASSIGNMENT
    assert evidence.assignment is assignment
    assert evidence.confidence is AssignmentConfidence.HIGH


def test_build_low_confidence_evidence() -> None:
    """A low-confidence accepted assignment should remain present."""

    reference_gene = make_reference_gene()
    assignment = make_assignment(
        reference_gene,
        AssignmentConfidence.LOW,
    )

    evidence = EvidenceBuilder().build(
        reference_gene,
        [assignment],
    )

    assert evidence.status is EvidenceStatus.PRESENT
    assert evidence.reason is EvidenceReason.LOW_CONFIDENCE


def test_build_rejected_evidence() -> None:
    """A rejected assignment should produce absent evidence."""

    reference_gene = make_reference_gene()
    assignment = make_assignment(
        reference_gene,
        AssignmentConfidence.REJECTED,
    )

    evidence = EvidenceBuilder().build(
        reference_gene,
        [assignment],
    )

    assert evidence.status is EvidenceStatus.ABSENT
    assert evidence.reason is EvidenceReason.REJECTED_ASSIGNMENT
    assert evidence.assignment is assignment


def test_build_evidence_without_matching_hit() -> None:
    """No matching assignment should produce absent evidence."""

    reference_gene = make_reference_gene()
    other_gene = make_reference_gene(
        gene_id="GSP0002",
        canonical_name="spoIIE",
    )
    other_assignment = make_assignment(
        other_gene,
        AssignmentConfidence.HIGH,
    )

    evidence = EvidenceBuilder().build(
        reference_gene,
        [other_assignment],
    )

    assert evidence.status is EvidenceStatus.ABSENT
    assert evidence.reason is EvidenceReason.NO_ACCEPTED_HIT
    assert evidence.assignment is None


def test_build_not_assessed_evidence() -> None:
    """An unsearched gene should not be classified as absent."""

    reference_gene = make_reference_gene()

    evidence = EvidenceBuilder().build(
        reference_gene,
        [],
        assessed=False,
    )

    assert evidence.status is EvidenceStatus.NOT_ASSESSED
    assert evidence.reason is EvidenceReason.NOT_SEARCHED
    assert evidence.assignment is None


def test_builder_selects_highest_confidence_assignment() -> None:
    """Confidence should take priority when selecting evidence."""

    reference_gene = make_reference_gene()
    low_assignment = make_assignment(
        reference_gene,
        AssignmentConfidence.LOW,
        bitscore=500.0,
        query_id="query_low",
    )
    high_assignment = make_assignment(
        reference_gene,
        AssignmentConfidence.HIGH,
        bitscore=200.0,
        query_id="query_high",
    )

    evidence = EvidenceBuilder().build(
        reference_gene,
        [low_assignment, high_assignment],
    )

    assert evidence.assignment is high_assignment


def test_builder_uses_bitscore_as_tiebreaker() -> None:
    """Bitscore should resolve assignments with equal confidence."""

    reference_gene = make_reference_gene()
    lower_bitscore = make_assignment(
        reference_gene,
        AssignmentConfidence.MEDIUM,
        bitscore=150.0,
        query_id="query_lower",
    )
    higher_bitscore = make_assignment(
        reference_gene,
        AssignmentConfidence.MEDIUM,
        bitscore=250.0,
        query_id="query_higher",
    )

    evidence = EvidenceBuilder().build(
        reference_gene,
        [lower_bitscore, higher_bitscore],
    )

    assert evidence.assignment is higher_bitscore
