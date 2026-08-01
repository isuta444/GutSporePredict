"""Tests for the knowledge registry."""

from gutsporepredict.knowledge.knowledge_registry import (
    KnowledgeRegistry,
)


def test_registry_load() -> None:
    """Registry should load all knowledge."""

    registry = KnowledgeRegistry.load()

    assert "spo0A" in registry.genes
    assert "SP001" in registry.modules

    assert (
        registry.modules["SP001"].required_gene_ids
        == ("spo0A",)
    )
