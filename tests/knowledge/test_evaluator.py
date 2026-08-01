"""Tests for module evaluation."""

from gutsporepredict.evidence.models import (
    EvidenceReason,
    EvidenceStatus,
    GeneEvidence,
)
from gutsporepredict.knowledge.evaluator import (
    ModuleEvaluator,
    ModuleStatus,
)
from gutsporepredict.knowledge.models import (
    GeneRequirement,
    KnowledgeEvidenceLevel,
    ModuleGene,
    ReferenceModule,
)
from gutsporepredict.reference.models import ReferenceGene


def make_reference_gene(gene_id: str) -> ReferenceGene:
    """Create a minimal reference gene for testing."""

    return ReferenceGene(
        gene_id=gene_id,
        canonical_name=gene_id,
        pathway="sporulation",
        module="test_module",
        stage="test",
        essentiality="required",
        phyletic_pattern="core",
        search_methods=("diamond",),
        description="Test reference gene.",
    )


def make_gene_evidence(
    gene_id: str,
    *,
    present: bool,
) -> GeneEvidence:
    """Create gene-level evidence for testing."""

    if present:
        status = EvidenceStatus.PRESENT
        reason = EvidenceReason.ACCEPTED_ASSIGNMENT
    else:
        status = EvidenceStatus.ABSENT
        reason = EvidenceReason.NO_ACCEPTED_HIT

    return GeneEvidence(
        reference_gene=make_reference_gene(gene_id),
        status=status,
        reason=reason,
    )


def make_module() -> ReferenceModule:
    """Create a module containing two required genes."""

    return ReferenceModule(
        module_id="test_module",
        name="Test module",
        pathway_id="sporulation",
        description="Module used for evaluator tests.",
        genes=(
            ModuleGene(
                gene_id="spo0A",
                requirement=GeneRequirement.REQUIRED,
            ),
            ModuleGene(
                gene_id="spoIIE",
                requirement=GeneRequirement.REQUIRED,
            ),
        ),
        evidence_level=KnowledgeEvidenceLevel.CURATED,
    )


def test_complete_module() -> None:
    """A module is complete when all required genes are present."""

    result = ModuleEvaluator().evaluate(
        make_module(),
        (
            make_gene_evidence("spo0A", present=True),
            make_gene_evidence("spoIIE", present=True),
        ),
    )

    assert result.status is ModuleStatus.COMPLETE
    assert result.required_present == 2
    assert result.required_total == 2
    assert result.score == 1.0


def test_partial_module() -> None:
    """A module is partial when some required genes are present."""

    result = ModuleEvaluator().evaluate(
        make_module(),
        (
            make_gene_evidence("spo0A", present=True),
            make_gene_evidence("spoIIE", present=False),
        ),
    )

    assert result.status is ModuleStatus.PARTIAL
    assert result.required_present == 1
    assert result.required_total == 2
    assert result.score == 0.5


def test_absent_module() -> None:
    """A module is absent when no required genes are present."""

    result = ModuleEvaluator().evaluate(
        make_module(),
        (
            make_gene_evidence("spo0A", present=False),
            make_gene_evidence("spoIIE", present=False),
        ),
    )

    assert result.status is ModuleStatus.ABSENT
    assert result.required_present == 0
    assert result.required_total == 2
    assert result.score == 0.0


def test_uncertain_required_gene_makes_module_uncertain() -> None:
    """An uncertain required gene prevents a definitive module call."""

    from gutsporepredict.evidence.models import (
        EvidenceReason,
        EvidenceStatus,
        GeneEvidence,
    )
    from gutsporepredict.reference.models import ReferenceGene

    module = make_module()

    uncertain_gene = ReferenceGene(
        gene_id="spoIIE",
        canonical_name="spoIIE",
        pathway="sporulation",
        module=module.module_id,
        stage="unknown",
        essentiality="unknown",
        phyletic_pattern="unknown",
        search_methods=("hmmer",),
        description="Test reference gene.",
    )

    result = ModuleEvaluator().evaluate(
        module,
        (
            make_gene_evidence("spo0A", present=True),
            GeneEvidence(
                reference_gene=uncertain_gene,
                status=EvidenceStatus.UNCERTAIN,
                reason=EvidenceReason.LOW_CONFIDENCE,
            ),
        ),
    )

    assert result.status is ModuleStatus.UNCERTAIN
    assert result.required_present == 1
    assert result.required_total == 2
