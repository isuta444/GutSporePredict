"""Tests for the gene loader."""

from gutsporepredict.knowledge.gene_loader import (
    load_gene,
    load_genes,
)


def test_load_spo0a() -> None:
    """spo0A should load correctly."""

    gene = load_gene("knowledge/genes/spo0A.yaml")

    assert gene.gene_id == "spo0A"
    assert gene.symbol == "spo0A"
    assert gene.full_name == "Stage 0 sporulation protein A"

def test_load_gene_directory() -> None:
    """All genes should load."""

    genes = load_genes("knowledge/genes")

    assert "spo0A" in genes
    assert genes["spo0A"].symbol == "spo0A"
