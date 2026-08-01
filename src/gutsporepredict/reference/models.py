"""Models for the GutSporePredict reference database."""

from dataclasses import dataclass, field


@dataclass(frozen=True)
class ReferenceGene:
    """Curated reference-gene definition."""

    gene_id: str
    canonical_name: str
    pathway: str
    module: str
    stage: str
    essentiality: str
    phyletic_pattern: str
    search_methods: tuple[str, ...]
    description: str


@dataclass(frozen=True)
class GeneAlias:
    """Alternative name associated with a reference gene."""

    gene_id: str
    alias: str


@dataclass(frozen=True)
class ReferenceDatabase:
    """Loaded reference genes and aliases."""

    genes: tuple[ReferenceGene, ...]
    aliases: tuple[GeneAlias, ...]

    _gene_index: dict[str, ReferenceGene] = field(
        init=False,
        repr=False,
    )
    _name_index: dict[str, ReferenceGene] = field(
        init=False,
        repr=False,
    )

    def __post_init__(self) -> None:
        """Build lookup indexes."""

        gene_index = {
            gene.gene_id: gene
            for gene in self.genes
        }

        name_index = {
            gene.canonical_name.lower(): gene
            for gene in self.genes
        }

        for alias in self.aliases:
            gene = gene_index.get(alias.gene_id)
            if gene is not None:
                name_index.setdefault(
                    alias.alias.lower(),
                    gene,
                )

        object.__setattr__(
            self,
            "_gene_index",
            gene_index,
        )
        object.__setattr__(
            self,
            "_name_index",
            name_index,
        )

    def gene_by_id(
        self,
        gene_id: str,
    ) -> ReferenceGene | None:
        """Return a reference gene by its stable identifier."""

        return self._gene_index.get(gene_id)

    def resolve_name(
        self,
        name: str,
    ) -> ReferenceGene | None:
        """Resolve a canonical name or alias."""

        return self._name_index.get(
            name.strip().lower()
        )
