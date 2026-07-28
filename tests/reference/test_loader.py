"""Tests for reference-database loading."""

from pathlib import Path

import pytest

from gutsporepredict.reference.exceptions import ReferenceLoadError
from gutsporepredict.reference.loader import ReferenceLoader


def test_load_reference_database() -> None:
    loader = ReferenceLoader()

    database = loader.load(
        "database/reference/genes.tsv",
        "database/reference/aliases.tsv",
    )

    assert len(database.genes) == 13
    assert len(database.aliases) == 12

    spo0a = database.gene_by_id("GSP0001")

    assert spo0a is not None
    assert spo0a.canonical_name == "spo0A"
    assert spo0a.essentiality == "essential"
    assert spo0a.search_methods == (
        "diamond",
        "hmmer",
    )


def test_resolve_canonical_name_and_alias() -> None:
    loader = ReferenceLoader()

    database = loader.load(
        "database/reference/genes.tsv",
        "database/reference/aliases.tsv",
    )

    assert database.resolve_name("spoIIAC") is not None
    assert database.resolve_name("sigF") is not None
    assert (
        database.resolve_name("sigF").gene_id
        == "GSP0004"
    )
    assert database.resolve_name("unknown") is None


def test_reject_missing_gene_columns(
    tmp_path: Path,
) -> None:
    genes_path = tmp_path / "genes.tsv"
    genes_path.write_text(
        "gene_id\tcanonical_name\n"
        "GSP0001\tspo0A\n",
        encoding="utf-8",
    )

    aliases_path = tmp_path / "aliases.tsv"
    aliases_path.write_text(
        "gene_id\talias\n",
        encoding="utf-8",
    )

    with pytest.raises(
        ReferenceLoadError,
        match="Missing columns",
    ):
        ReferenceLoader().load(
            genes_path,
            aliases_path,
        )
