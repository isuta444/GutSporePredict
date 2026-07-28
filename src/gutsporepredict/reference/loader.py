"""Load the GutSporePredict reference database."""

import csv
from collections.abc import Sequence
from pathlib import Path

from gutsporepredict.reference.exceptions import ReferenceLoadError
from gutsporepredict.reference.models import (
    GeneAlias,
    ReferenceDatabase,
    ReferenceGene,
)

GENE_COLUMNS = (
    "gene_id",
    "canonical_name",
    "pathway",
    "module",
    "stage",
    "essentiality",
    "phyletic_pattern",
    "search_methods",
    "description",
)

ALIAS_COLUMNS = (
    "gene_id",
    "alias",
)


class ReferenceLoader:
    """Load genes and aliases from tab-separated files."""

    def load(
        self,
        genes_path: str | Path,
        aliases_path: str | Path,
    ) -> ReferenceDatabase:
        """Load a complete reference database."""

        genes = self.load_genes(genes_path)
        aliases = self.load_aliases(aliases_path)

        return ReferenceDatabase(
            genes=tuple(genes),
            aliases=tuple(aliases),
        )

    def load_genes(
        self,
        path: str | Path,
    ) -> list[ReferenceGene]:
        """Load reference genes from a TSV file."""

        path = self._validate_file(path, "Genes")

        with path.open(
            encoding="utf-8",
            newline="",
        ) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            self._validate_columns(
                reader.fieldnames,
                GENE_COLUMNS,
                path,
            )

            genes: list[ReferenceGene] = []

            for line_number, row in enumerate(reader, start=2):
                search_methods = tuple(
                    method.strip().lower()
                    for method in row["search_methods"].split(",")
                    if method.strip()
                )

                genes.append(
                    ReferenceGene(
                        gene_id=row["gene_id"].strip(),
                        canonical_name=row[
                            "canonical_name"
                        ].strip(),
                        pathway=row["pathway"].strip().lower(),
                        module=row["module"].strip(),
                        stage=row["stage"].strip(),
                        essentiality=row[
                            "essentiality"
                        ].strip().lower(),
                        phyletic_pattern=row[
                            "phyletic_pattern"
                        ].strip().lower(),
                        search_methods=search_methods,
                        description=row["description"].strip(),
                    )
                )

                if not genes[-1].gene_id:
                    raise ReferenceLoadError(
                        f"Empty gene_id at {path}:{line_number}"
                    )

        return genes

    def load_aliases(
        self,
        path: str | Path,
    ) -> list[GeneAlias]:
        """Load gene aliases from a TSV file."""

        path = self._validate_file(path, "Aliases")

        with path.open(
            encoding="utf-8",
            newline="",
        ) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            self._validate_columns(
                reader.fieldnames,
                ALIAS_COLUMNS,
                path,
            )

            aliases: list[GeneAlias] = []

            for line_number, row in enumerate(reader, start=2):
                alias = GeneAlias(
                    gene_id=row["gene_id"].strip(),
                    alias=row["alias"].strip(),
                )

                if not alias.gene_id or not alias.alias:
                    raise ReferenceLoadError(
                        "Empty alias field at "
                        f"{path}:{line_number}"
                    )

                aliases.append(alias)

        return aliases

    @staticmethod
    def _validate_file(
        path: str | Path,
        label: str,
    ) -> Path:
        path = Path(path)

        if not path.exists():
            raise ReferenceLoadError(
                f"{label} file does not exist: {path}"
            )

        if not path.is_file():
            raise ReferenceLoadError(
                f"{label} path is not a file: {path}"
            )

        return path

    @staticmethod
    def _validate_columns(
        fieldnames: Sequence[str] | None,
        required_columns: tuple[str, ...],
        path: Path,
    ) -> None:
        available = set(fieldnames or [])
        missing = [
            column
            for column in required_columns
            if column not in available
        ]

        if missing:
            raise ReferenceLoadError(
                f"Missing columns in {path}: "
                f"{', '.join(missing)}"
            )
