"""Central registry for the GutSporePredict knowledge base."""

from dataclasses import dataclass
from pathlib import Path

from gutsporepredict.knowledge.gene_loader import load_genes
from gutsporepredict.knowledge.models import (
    ReferenceGene,
    ReferenceModule,
)
from gutsporepredict.knowledge.module_loader import load_module


@dataclass(frozen=True)
class KnowledgeRegistry:
    """Central access point for biological knowledge."""

    genes: dict[str, ReferenceGene]
    modules: dict[str, ReferenceModule]

    @classmethod
    def load(
        cls,
        knowledge_dir: str | Path = "knowledge",
    ) -> "KnowledgeRegistry":
        """Load all knowledge."""

        knowledge_dir = Path(knowledge_dir)

        genes = load_genes(
            knowledge_dir / "genes"
        )

        modules: dict[str, ReferenceModule] = {}

        for path in sorted(
            (knowledge_dir / "modules").glob("*.yaml")
        ):
            module = load_module(path)
            modules[module.module_id] = module

        return cls(
            genes=genes,
            modules=modules,
        )
