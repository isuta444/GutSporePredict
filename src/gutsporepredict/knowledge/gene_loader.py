"""Load reference genes from YAML files."""

from pathlib import Path

import yaml

from gutsporepredict.knowledge.models import ReferenceGene


def load_gene(path: str | Path) -> ReferenceGene:
    """Load one reference gene from YAML."""

    path = Path(path)

    with path.open(encoding="utf-8") as fh:
        data = yaml.safe_load(fh)

    return ReferenceGene.from_dict(data)


def load_genes(directory: str | Path) -> dict[str, ReferenceGene]:
    """Load all reference genes in a directory."""

    directory = Path(directory)

    genes: dict[str, ReferenceGene] = {}

    for path in sorted(directory.glob("*.yaml")):
        gene = load_gene(path)

        if gene.gene_id in genes:
            raise ValueError(
                f"Duplicate gene_id: {gene.gene_id}"
            )

        genes[gene.gene_id] = gene

    return genes
