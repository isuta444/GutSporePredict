"""Tests for reference-database validation."""

import pytest

from gutsporepredict.reference.exceptions import (
    ReferenceValidationError,
)
from gutsporepredict.reference.loader import ReferenceLoader
from gutsporepredict.reference.models import (
    GeneAlias,
    ReferenceDatabase,
    ReferenceGene,
)
from gutsporepredict.reference.validator import (
    ReferenceValidator,
)


def make_gene(
    gene_id: str = "GSP0001",
    canonical_name: str = "spo0A",
) -> ReferenceGene:
    return ReferenceGene(
        gene_id=gene_id,
        canonical_name=canonical_name,
        pathway="sporulation",
        module="initiation",
        stage="stage0",
        essentiality="essential",
        phyletic_pattern="core",
        search_methods=("diamond", "hmmer"),
        description="Test gene",
    )


def test_validate_reference_database() -> None:
    database = ReferenceLoader().load(
        "database/reference/genes.tsv",
        "database/reference/aliases.tsv",
    )

    ReferenceValidator().validate(database)


def test_reject_duplicate_gene_id() -> None:
    database = ReferenceDatabase(
        genes=(
            make_gene(),
            make_gene(canonical_name="spo0B"),
        ),
        aliases=(),
    )

    with pytest.raises(
        ReferenceValidationError,
        match="Duplicate gene_id",
    ):
        ReferenceValidator().validate(database)


def test_reject_unknown_alias_gene_id() -> None:
    database = ReferenceDatabase(
        genes=(make_gene(),),
        aliases=(
            GeneAlias(
                gene_id="GSP9999",
                alias="missing_gene",
            ),
        ),
    )

    with pytest.raises(
        ReferenceValidationError,
        match="unknown gene_id",
    ):
        ReferenceValidator().validate(database)


def test_reject_invalid_search_method() -> None:
    invalid_gene = ReferenceGene(
        gene_id="GSP0001",
        canonical_name="spo0A",
        pathway="sporulation",
        module="initiation",
        stage="stage0",
        essentiality="essential",
        phyletic_pattern="core",
        search_methods=("blast",),
        description="Test gene",
    )

    database = ReferenceDatabase(
        genes=(invalid_gene,),
        aliases=(),
    )

    with pytest.raises(
        ReferenceValidationError,
        match="Invalid search methods",
    ):
        ReferenceValidator().validate(database)
