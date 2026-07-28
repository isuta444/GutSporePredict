"""Validation of the GutSporePredict reference database."""

import re

from gutsporepredict.reference.exceptions import (
    ReferenceValidationError,
)
from gutsporepredict.reference.models import ReferenceDatabase

GENE_ID_PATTERN = re.compile(r"^GSP\d{4,}$")

ALLOWED_PATHWAYS = {
    "sporulation",
    "germination",
}

ALLOWED_ESSENTIALITY = {
    "essential",
    "accessory",
    "lineage_specific",
    "unknown",
}

ALLOWED_PHYLETIC_PATTERNS = {
    "core",
    "strict",
    "extended",
    "evolutionary",
    "lineage_specific",
    "unknown",
}

ALLOWED_SEARCH_METHODS = {
    "diamond",
    "hmmer",
}


class ReferenceValidator:
    """Validate loaded genes and aliases."""

    def validate(
        self,
        database: ReferenceDatabase,
    ) -> None:
        """Raise an exception if reference data are invalid."""

        errors: list[str] = []

        self._validate_genes(database, errors)
        self._validate_aliases(database, errors)

        if errors:
            formatted = "\n".join(
                f"- {error}" for error in errors
            )
            raise ReferenceValidationError(
                "Reference database validation failed:\n"
                f"{formatted}"
            )

    @staticmethod
    def _validate_genes(
        database: ReferenceDatabase,
        errors: list[str],
    ) -> None:
        gene_ids: set[str] = set()
        canonical_names: set[str] = set()

        for gene in database.genes:
            if gene.gene_id in gene_ids:
                errors.append(
                    f"Duplicate gene_id: {gene.gene_id}"
                )

            gene_ids.add(gene.gene_id)

            normalized_name = gene.canonical_name.lower()

            if normalized_name in canonical_names:
                errors.append(
                    "Duplicate canonical_name: "
                    f"{gene.canonical_name}"
                )

            canonical_names.add(normalized_name)

            if not GENE_ID_PATTERN.fullmatch(gene.gene_id):
                errors.append(
                    f"Invalid gene_id format: {gene.gene_id}"
                )

            if not gene.canonical_name:
                errors.append(
                    f"Empty canonical_name: {gene.gene_id}"
                )

            if gene.pathway not in ALLOWED_PATHWAYS:
                errors.append(
                    f"Invalid pathway for {gene.gene_id}: "
                    f"{gene.pathway}"
                )

            if gene.essentiality not in ALLOWED_ESSENTIALITY:
                errors.append(
                    f"Invalid essentiality for {gene.gene_id}: "
                    f"{gene.essentiality}"
                )

            if (
                gene.phyletic_pattern
                not in ALLOWED_PHYLETIC_PATTERNS
            ):
                errors.append(
                    "Invalid phyletic_pattern for "
                    f"{gene.gene_id}: "
                    f"{gene.phyletic_pattern}"
                )

            if not gene.module:
                errors.append(
                    f"Empty module: {gene.gene_id}"
                )

            if not gene.stage:
                errors.append(
                    f"Empty stage: {gene.gene_id}"
                )

            if not gene.search_methods:
                errors.append(
                    f"No search methods: {gene.gene_id}"
                )

            invalid_methods = (
                set(gene.search_methods)
                - ALLOWED_SEARCH_METHODS
            )

            if invalid_methods:
                errors.append(
                    f"Invalid search methods for {gene.gene_id}: "
                    f"{', '.join(sorted(invalid_methods))}"
                )

    @staticmethod
    def _validate_aliases(
        database: ReferenceDatabase,
        errors: list[str],
    ) -> None:
        valid_gene_ids = {
            gene.gene_id for gene in database.genes
        }
        canonical_names = {
            gene.canonical_name.lower()
            for gene in database.genes
        }
        alias_names: set[str] = set()

        for alias in database.aliases:
            normalized_alias = alias.alias.lower()

            if alias.gene_id not in valid_gene_ids:
                errors.append(
                    "Alias references unknown gene_id: "
                    f"{alias.gene_id}"
                )

            if normalized_alias in alias_names:
                errors.append(
                    f"Duplicate alias: {alias.alias}"
                )

            alias_names.add(normalized_alias)

            matching_gene = database.gene_by_id(
                alias.gene_id
            )

            if (
                normalized_alias in canonical_names
                and (
                    matching_gene is None
                    or matching_gene.canonical_name.lower()
                    != normalized_alias
                )
            ):
                errors.append(
                    "Alias conflicts with canonical name: "
                    f"{alias.alias}"
                )
