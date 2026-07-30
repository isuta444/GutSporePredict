"""Tests for knowledge-base data models."""

import pytest

from gutsporepredict.knowledge.exceptions import KnowledgeModelError
from gutsporepredict.knowledge.models import (
    GeneRequirement,
    KnowledgeEvidenceLevel,
    ModuleGene,
    ReferenceModule,
    ReferencePathway,
    ReferencePhenotype,
)


def make_module_gene(
    gene_id: str,
    requirement: GeneRequirement,
    *,
    weight: float = 1.0,
) -> ModuleGene:
    """Create a module-gene definition for testing."""

    return ModuleGene(
        gene_id=gene_id,
        requirement=requirement,
        weight=weight,
    )


def test_module_gene_accepts_positive_weight() -> None:
    """A module gene should retain its biological weight."""

    gene = make_module_gene(
        "GSP0001",
        GeneRequirement.REQUIRED,
        weight=3.0,
    )

    assert gene.gene_id == "GSP0001"
    assert gene.weight == 3.0


def test_module_gene_rejects_non_positive_weight() -> None:
    """Gene weights must be greater than zero."""

    with pytest.raises(
        KnowledgeModelError,
        match="greater than zero",
    ):
        make_module_gene(
            "GSP0001",
            GeneRequirement.REQUIRED,
            weight=0.0,
        )


def test_module_gene_from_dict() -> None:
    """ModuleGene should be created from raw knowledge data."""

    gene = ModuleGene.from_dict(
        {
            "gene_id": "spo0A",
            "requirement": "required",
            "weight": 5,
            "taxonomic_scope": ["Bacillota"],
            "notes": "Master regulator.",
        }
    )

    assert gene.gene_id == "spo0A"
    assert gene.requirement is GeneRequirement.REQUIRED
    assert gene.weight == 5.0
    assert gene.taxonomic_scope == ("Bacillota",)
    assert gene.notes == "Master regulator."


def test_module_gene_from_dict_uses_default_weight() -> None:
    """Module genes should have a default weight of one."""

    gene = ModuleGene.from_dict(
        {
            "gene_id": "spo0A",
            "requirement": "required",
        }
    )

    assert gene.weight == 1.0


def test_module_gene_from_dict_rejects_invalid_requirement() -> None:
    """Unknown gene requirements should be rejected."""

    with pytest.raises(
        KnowledgeModelError,
        match="must be one of",
    ):
        ModuleGene.from_dict(
            {
                "gene_id": "spo0A",
                "requirement": "essential",
            }
        )


def test_reference_module_groups_genes_by_requirement() -> None:
    """A module should expose genes grouped by requirement."""

    module = ReferenceModule(
        module_id="sporulation_initiation",
        name="Sporulation initiation",
        pathway_id="sporulation",
        description="Initiation of the sporulation programme.",
        genes=(
            make_module_gene(
                "GSP0001",
                GeneRequirement.REQUIRED,
                weight=5.0,
            ),
            make_module_gene(
                "GSP0002",
                GeneRequirement.SUPPORTING,
                weight=2.0,
            ),
            make_module_gene(
                "GSP0003",
                GeneRequirement.OPTIONAL,
            ),
        ),
        evidence_level=KnowledgeEvidenceLevel.CURATED,
    )

    assert module.required_gene_ids == ("GSP0001",)
    assert module.supporting_gene_ids == ("GSP0002",)
    assert module.optional_gene_ids == ("GSP0003",)
    assert module.total_weight == 8.0


def test_reference_module_rejects_duplicate_genes() -> None:
    """A gene must not be listed twice in one module."""

    with pytest.raises(
        KnowledgeModelError,
        match="duplicate gene",
    ):
        ReferenceModule(
            module_id="sporulation_initiation",
            name="Sporulation initiation",
            pathway_id="sporulation",
            description="Test module.",
            genes=(
                make_module_gene(
                    "GSP0001",
                    GeneRequirement.REQUIRED,
                ),
                make_module_gene(
                    "GSP0001",
                    GeneRequirement.OPTIONAL,
                ),
            ),
            evidence_level=KnowledgeEvidenceLevel.CURATED,
        )


def test_reference_module_from_dict() -> None:
    """ReferenceModule should parse nested gene records."""

    module = ReferenceModule.from_dict(
        {
            "module_id": "sporulation_initiation",
            "name": "Sporulation initiation",
            "pathway_id": "sporulation",
            "description": "Initiation of sporulation.",
            "genes": [
                {
                    "gene_id": "spo0A",
                    "requirement": "required",
                    "weight": 5.0,
                },
                {
                    "gene_id": "spo0F",
                    "requirement": "supporting",
                },
            ],
            "evidence_level": "curated",
            "literature": ["PMID:12345678"],
            "taxonomic_scope": ["Bacillota"],
        }
    )

    assert module.module_id == "sporulation_initiation"
    assert module.required_gene_ids == ("spo0A",)
    assert module.supporting_gene_ids == ("spo0F",)
    assert module.evidence_level is KnowledgeEvidenceLevel.CURATED
    assert module.literature == ("PMID:12345678",)


def test_reference_module_from_dict_rejects_non_mapping_gene() -> None:
    """Nested module genes must be mappings."""

    with pytest.raises(
        KnowledgeModelError,
        match="must contain only mappings",
    ):
        ReferenceModule.from_dict(
            {
                "module_id": "sporulation_initiation",
                "name": "Sporulation initiation",
                "pathway_id": "sporulation",
                "description": "Test module.",
                "genes": ["spo0A"],
                "evidence_level": "curated",
            }
        )


def test_reference_pathway_retains_module_order() -> None:
    """Pathway module order should follow biological progression."""

    pathway = ReferencePathway(
        pathway_id="sporulation",
        name="Sporulation",
        description="Formation of a dormant endospore.",
        module_ids=(
            "sporulation_initiation",
            "asymmetric_division",
            "engulfment",
            "cortex_formation",
            "coat_assembly",
        ),
        evidence_level=KnowledgeEvidenceLevel.CURATED,
    )

    assert pathway.module_ids[0] == "sporulation_initiation"
    assert pathway.module_ids[-1] == "coat_assembly"


def test_reference_pathway_rejects_duplicate_modules() -> None:
    """A pathway must not contain duplicate module identifiers."""

    with pytest.raises(
        KnowledgeModelError,
        match="duplicate module",
    ):
        ReferencePathway(
            pathway_id="sporulation",
            name="Sporulation",
            description="Test pathway.",
            module_ids=(
                "sporulation_initiation",
                "sporulation_initiation",
            ),
            evidence_level=KnowledgeEvidenceLevel.CURATED,
        )


def test_reference_pathway_from_dict() -> None:
    """ReferencePathway should be created from raw data."""

    pathway = ReferencePathway.from_dict(
        {
            "pathway_id": "sporulation",
            "name": "Sporulation",
            "description": "Formation of a mature endospore.",
            "module_ids": [
                "sporulation_initiation",
                "asymmetric_division",
            ],
            "evidence_level": "experimental",
            "literature": ["PMID:87654321"],
        }
    )

    assert pathway.pathway_id == "sporulation"
    assert pathway.module_ids == (
        "sporulation_initiation",
        "asymmetric_division",
    )
    assert (
        pathway.evidence_level
        is KnowledgeEvidenceLevel.EXPERIMENTAL
    )


def test_reference_phenotype_links_to_pathways() -> None:
    """A phenotype should retain its supporting pathways."""

    phenotype = ReferencePhenotype(
        phenotype_id="endospore_formation",
        name="Endospore formation",
        description="Capacity to produce a mature endospore.",
        pathway_ids=("sporulation",),
        evidence_level=KnowledgeEvidenceLevel.CURATED,
    )

    assert phenotype.pathway_ids == ("sporulation",)


def test_reference_phenotype_requires_a_pathway() -> None:
    """A phenotype without a pathway definition is invalid."""

    with pytest.raises(
        KnowledgeModelError,
        match="at least one pathway",
    ):
        ReferencePhenotype(
            phenotype_id="endospore_formation",
            name="Endospore formation",
            description="Test phenotype.",
            pathway_ids=(),
            evidence_level=KnowledgeEvidenceLevel.PROVISIONAL,
        )


def test_reference_phenotype_from_dict() -> None:
    """ReferencePhenotype should be created from raw data."""

    phenotype = ReferencePhenotype.from_dict(
        {
            "phenotype_id": "endospore_formation",
            "name": "Endospore formation",
            "description": "Capacity to form an endospore.",
            "pathway_ids": ["sporulation"],
            "evidence_level": "curated",
        }
    )

    assert phenotype.phenotype_id == "endospore_formation"
    assert phenotype.pathway_ids == ("sporulation",)


def test_from_dict_rejects_missing_required_string() -> None:
    """Missing required text fields should raise a model error."""

    with pytest.raises(
        KnowledgeModelError,
        match="ReferencePathway.name",
    ):
        ReferencePathway.from_dict(
            {
                "pathway_id": "sporulation",
                "description": "Test pathway.",
                "module_ids": ["sporulation_initiation"],
                "evidence_level": "curated",
            }
        )


def test_from_dict_rejects_string_instead_of_sequence() -> None:
    """A string must not be interpreted as a sequence of IDs."""

    with pytest.raises(
        KnowledgeModelError,
        match="must be a sequence of strings",
    ):
        ReferencePhenotype.from_dict(
            {
                "phenotype_id": "endospore_formation",
                "name": "Endospore formation",
                "description": "Test phenotype.",
                "pathway_ids": "sporulation",
                "evidence_level": "curated",
            }
        )
