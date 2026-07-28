"""Models for the GutSporePredict reference database."""

from dataclasses import dataclass


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

    def gene_by_id(self, gene_id: str) -> ReferenceGene | None:
        """Return a reference gene by its stable identifier."""

        return next(
            (
                gene
                for gene in self.genes
                if gene.gene_id == gene_id
            ),
            None,
        )

    def resolve_name(self, name: str) -> ReferenceGene | None:
        """Resolve a canonical name or alias to a reference gene."""

        normalized_name = name.strip().lower()

        for gene in self.genes:
            if gene.canonical_name.lower() == normalized_name:
                return gene

        alias_gene_ids = {
            alias.gene_id
            for alias in self.aliases
            if alias.alias.lower() == normalized_name
        }

        if len(alias_gene_ids) != 1:
            return None

        return self.gene_by_id(next(iter(alias_gene_ids)))
