"""Integrity tests for the biological knowledge base."""

from gutsporepredict.knowledge.gene_loader import load_genes
from gutsporepredict.knowledge.module_loader import load_module


def test_sp001_gene_references_exist() -> None:
    """Every gene referenced by SP001 should exist in the gene database."""

    genes = load_genes("knowledge/genes")
    module = load_module("knowledge/modules/SP001.yaml")

    missing_gene_ids = sorted(
        gene.gene_id
        for gene in module.genes
        if gene.gene_id not in genes
    )

    assert missing_gene_ids == []
