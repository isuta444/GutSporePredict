"""Tests for the module loader."""

from gutsporepredict.knowledge.module_loader import load_module


def test_load_sp001() -> None:
    """SP001 should load correctly."""

    module = load_module("knowledge/modules/SP001.yaml")

    assert module.module_id == "SP001"
    assert module.name == "Sporulation initiation"
    assert module.pathway_id == "sporulation"

    assert "spo0A" in module.required_gene_ids
    assert "sigH" in module.supporting_gene_ids
    assert module.optional_gene_ids == ()
