"""Tests for evidence data models."""

from gutsporepredict.evidence.models import (
    EvidenceReason,
    EvidenceStatus,
    GeneEvidence,
)
from gutsporepredict.reference.models import ReferenceGene


def make_reference_gene() -> ReferenceGene:
    """Create a reference gene for testing."""

    return ReferenceGene(
        gene_id="GSP0001",
        canonical_name="spo0A",
        pathway="sporulation",
        module="initiation",
        stage="stage_0",
        essentiality="core",
        phyletic_pattern="broad",
        search_methods=("diamond", "hmmer"),
        description="Master regulator of sporulation.",
    )


def test_present_evidence_properties() -> None:
    """Present evidence should expose reference-gene metadata."""

    evidence = GeneEvidence(
        reference_gene=make_reference_gene(),
        status=EvidenceStatus.PRESENT,
        reason=EvidenceReason.ACCEPTED_ASSIGNMENT,
    )

    assert evidence.present
    assert evidence.assessed
    assert evidence.gene_id == "GSP0001"
    assert evidence.canonical_name == "spo0A"
    assert evidence.confidence is None


def test_not_assessed_evidence() -> None:
    """Not-assessed evidence should be distinguishable from absence."""

    evidence = GeneEvidence(
        reference_gene=make_reference_gene(),
        status=EvidenceStatus.NOT_ASSESSED,
        reason=EvidenceReason.NOT_SEARCHED,
    )

    assert not evidence.present
    assert not evidence.assessed
