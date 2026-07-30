"""Knowledge database for curated biological reference models."""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path
from typing import TypeVar

from gutsporepredict.knowledge.exceptions import (
    KnowledgeDatabaseError,
    KnowledgeFormatError,
)
from gutsporepredict.knowledge.loader import load_yaml
from gutsporepredict.knowledge.models import (
    KnowledgeRecord,
    ReferenceGene,
    ReferenceModule,
    ReferencePathway,
    ReferencePhenotype,
)

T = TypeVar("T")


class KnowledgeDatabase:
    """Store and validate curated biological knowledge."""

    def __init__(
        self,
        *,
        modules: Iterable[ReferenceModule],
        pathways: Iterable[ReferencePathway],
        phenotypes: Iterable[ReferencePhenotype],
        genes: Iterable[ReferenceGene] = (),
    ) -> None:
        """Initialize the knowledge database."""

        self._genes = self._index_by_id(
            genes,
            "gene_id",
            "gene",
        )
        self._modules = self._index_by_id(
            modules,
            "module_id",
            "module",
        )
        self._pathways = self._index_by_id(
            pathways,
            "pathway_id",
            "pathway",
        )
        self._phenotypes = self._index_by_id(
            phenotypes,
            "phenotype_id",
            "phenotype",
        )

        self._validate_references()

    @classmethod
    def from_directory(
        cls,
        directory: str | Path,
    ) -> KnowledgeDatabase:
        """Load a knowledge database from four YAML files."""

        root = Path(directory)

        genes = tuple(
            ReferenceGene.from_dict(record)
            for record in cls._load_records(
                root / "genes.yaml",
                "genes",
            )
        )
        modules = tuple(
            ReferenceModule.from_dict(record)
            for record in cls._load_records(
                root / "modules.yaml",
                "modules",
            )
        )
        pathways = tuple(
            ReferencePathway.from_dict(record)
            for record in cls._load_records(
                root / "pathways.yaml",
                "pathways",
            )
        )
        phenotypes = tuple(
            ReferencePhenotype.from_dict(record)
            for record in cls._load_records(
                root / "phenotypes.yaml",
                "phenotypes",
            )
        )

        return cls(
            genes=genes,
            modules=modules,
            pathways=pathways,
            phenotypes=phenotypes,
        )

    @property
    def genes(self) -> dict[str, ReferenceGene]:
        """Return genes indexed by gene identifier."""

        return self._genes

    @property
    def modules(self) -> dict[str, ReferenceModule]:
        """Return modules indexed by module identifier."""

        return self._modules

    @property
    def pathways(self) -> dict[str, ReferencePathway]:
        """Return pathways indexed by pathway identifier."""

        return self._pathways

    @property
    def phenotypes(self) -> dict[str, ReferencePhenotype]:
        """Return phenotypes indexed by phenotype identifier."""

        return self._phenotypes

    def statistics(self) -> dict[str, int]:
        """Return counts of genes and reference models."""

        if self._genes:
            gene_count = len(self._genes)
        else:
            gene_count = len(
                {
                    gene.gene_id
                    for module in self._modules.values()
                    for gene in module.genes
                }
            )

        return {
            "genes": gene_count,
            "modules": len(self._modules),
            "pathways": len(self._pathways),
            "phenotypes": len(self._phenotypes),
        }

    @staticmethod
    def _load_records(
        path: Path,
        collection_name: str,
    ) -> tuple[KnowledgeRecord, ...]:
        """Load a sequence of records from one YAML collection."""

        data = load_yaml(path)
        value = data.get(collection_name)

        if isinstance(value, str) or not isinstance(value, Sequence):
            raise KnowledgeFormatError(
                path,
                f"'{collection_name}' must be a sequence of mappings.",
            )

        records: list[KnowledgeRecord] = []

        for index, item in enumerate(value):
            if not isinstance(item, Mapping):
                raise KnowledgeFormatError(
                    path,
                    f"'{collection_name}[{index}]' must be a mapping.",
                )

            if not all(isinstance(key, str) for key in item):
                raise KnowledgeFormatError(
                    path,
                    f"'{collection_name}[{index}]' must use string keys.",
                )

            records.append(item)

        return tuple(records)

    @staticmethod
    def _index_by_id(
        items: Iterable[T],
        attribute: str,
        item_type: str,
    ) -> dict[str, T]:
        """Index objects by an identifier attribute."""

        indexed: dict[str, T] = {}

        for item in items:
            identifier = getattr(item, attribute)

            if identifier in indexed:
                raise KnowledgeDatabaseError(
                    f"Duplicate {item_type} identifier: {identifier}"
                )

            indexed[identifier] = item

        return indexed

    def has_gene(self, gene_id: str) -> bool:
        """Return True if the gene exists."""

        return gene_id in self._genes

    def get_gene(self, gene_id: str) -> ReferenceGene:
        """Return a gene by its identifier."""

        try:
            return self._genes[gene_id]
        except KeyError as exc:
            raise KnowledgeDatabaseError(
                f"Unknown gene: {gene_id}"
            ) from exc

    def has_module(self, module_id: str) -> bool:
        """Return True if the module exists."""

        return module_id in self._modules

    def get_module(self, module_id: str) -> ReferenceModule:
        """Return a module by its identifier."""

        try:
            return self._modules[module_id]
        except KeyError as exc:
            raise KnowledgeDatabaseError(
                f"Unknown module: {module_id}"
            ) from exc

    def get_pathway(self, pathway_id: str) -> ReferencePathway:
        """Return a pathway by its identifier."""

        try:
            return self._pathways[pathway_id]
        except KeyError as exc:
            raise KnowledgeDatabaseError(
                f"Unknown pathway: {pathway_id}"
            ) from exc

    def get_phenotype(
        self,
        phenotype_id: str,
    ) -> ReferencePhenotype:
        """Return a phenotype by its identifier."""

        try:
            return self._phenotypes[phenotype_id]
        except KeyError as exc:
            raise KnowledgeDatabaseError(
                f"Unknown phenotype: {phenotype_id}"
            ) from exc

    def _validate_references(self) -> None:
        """Validate cross-references between knowledge objects."""

        if self._genes:
            for module in self._modules.values():
                for module_gene in module.genes:
                    if module_gene.gene_id not in self._genes:
                        raise KnowledgeDatabaseError(
                            f"Module '{module.module_id}' references "
                            f"unknown gene '{module_gene.gene_id}'."
                        )

        for module in self._modules.values():
            if module.pathway_id not in self._pathways:
                raise KnowledgeDatabaseError(
                    f"Module '{module.module_id}' references "
                    f"unknown pathway '{module.pathway_id}'."
                )

        for pathway in self._pathways.values():
            for module_id in pathway.module_ids:
                if module_id not in self._modules:
                    raise KnowledgeDatabaseError(
                        f"Pathway '{pathway.pathway_id}' references "
                        f"unknown module '{module_id}'."
                    )

        for phenotype in self._phenotypes.values():
            for pathway_id in phenotype.pathway_ids:
                if pathway_id not in self._pathways:
                    raise KnowledgeDatabaseError(
                        f"Phenotype '{phenotype.phenotype_id}' "
                        f"references unknown pathway '{pathway_id}'."
                    )
