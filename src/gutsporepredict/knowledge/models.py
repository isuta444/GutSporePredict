"""Data models for the GutSporePredict knowledge base."""

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from enum import Enum
from typing import TypeAlias, TypeVar

from gutsporepredict.knowledge.exceptions import KnowledgeModelError

KnowledgeRecord: TypeAlias = Mapping[str, object]
EnumType = TypeVar("EnumType", bound=Enum)


class KnowledgeEvidenceLevel(str, Enum):
    """Strength of biological support for a knowledge-base entry."""

    EXPERIMENTAL = "experimental"
    CURATED = "curated"
    COMPUTATIONAL = "computational"
    PROVISIONAL = "provisional"


class GeneRequirement(str, Enum):
    """Role of a gene within a biological module."""

    REQUIRED = "required"
    SUPPORTING = "supporting"
    OPTIONAL = "optional"


def _required_string(
    data: KnowledgeRecord,
    field: str,
    model_name: str,
) -> str:
    """Return a required non-empty string field."""

    value = data.get(field)

    if not isinstance(value, str) or not value.strip():
        raise KnowledgeModelError(
            f"{model_name}.{field} must be a non-empty string."
        )

    return value.strip()


def _optional_string(
    data: KnowledgeRecord,
    field: str,
    model_name: str,
) -> str | None:
    """Return an optional string field."""

    value = data.get(field)

    if value is None:
        return None

    if not isinstance(value, str):
        raise KnowledgeModelError(
            f"{model_name}.{field} must be a string or null."
        )

    stripped = value.strip()
    return stripped or None


def _string_tuple(
    data: KnowledgeRecord,
    field: str,
    model_name: str,
) -> tuple[str, ...]:
    """Convert a sequence of strings into a tuple."""

    value = data.get(field, ())

    if value is None:
        return ()

    if isinstance(value, str) or not isinstance(value, Sequence):
        raise KnowledgeModelError(
            f"{model_name}.{field} must be a sequence of strings."
        )

    result: list[str] = []

    for item in value:
        if not isinstance(item, str) or not item.strip():
            raise KnowledgeModelError(
                f"{model_name}.{field} must contain only "
                "non-empty strings."
            )

        result.append(item.strip())

    return tuple(result)


def _enum_value(
    data: KnowledgeRecord,
    field: str,
    enum_type: type[EnumType],
    model_name: str,
) -> EnumType:
    """Convert a string field into an enum value."""

    value = data.get(field)

    if not isinstance(value, str):
        raise KnowledgeModelError(
            f"{model_name}.{field} must be a string."
        )

    try:
        return enum_type(value.strip().lower())
    except ValueError as exc:
        allowed = ", ".join(member.value for member in enum_type)

        raise KnowledgeModelError(
            f"{model_name}.{field} must be one of: {allowed}."
        ) from exc


def _positive_float(
    data: KnowledgeRecord,
    field: str,
    model_name: str,
    *,
    default: float,
) -> float:
    """Return a positive floating-point value."""

    value = data.get(field, default)

    if isinstance(value, bool) or not isinstance(value, int | float):
        raise KnowledgeModelError(
            f"{model_name}.{field} must be a number."
        )

    result = float(value)

    if result <= 0:
        raise KnowledgeModelError(
            f"{model_name}.{field} must be greater than zero."
        )

    return result


def _record_sequence(
    data: KnowledgeRecord,
    field: str,
    model_name: str,
) -> tuple[KnowledgeRecord, ...]:
    """Return a required sequence of mapping records."""

    value = data.get(field)

    if isinstance(value, str) or not isinstance(value, Sequence):
        raise KnowledgeModelError(
            f"{model_name}.{field} must be a sequence of mappings."
        )

    records: list[KnowledgeRecord] = []

    for item in value:
        if not isinstance(item, Mapping):
            raise KnowledgeModelError(
                f"{model_name}.{field} must contain only mappings."
            )

        records.append(item)

    return tuple(records)


@dataclass(frozen=True)
class ModuleGene:
    """Membership and importance of one gene in a module."""

    gene_id: str
    requirement: GeneRequirement
    weight: float = 1.0
    taxonomic_scope: tuple[str, ...] = ()
    notes: str | None = None

    def __post_init__(self) -> None:
        """Validate the module-gene definition."""

        if not self.gene_id.strip():
            raise KnowledgeModelError(
                "ModuleGene.gene_id must not be empty."
            )

        if self.weight <= 0:
            raise KnowledgeModelError(
                "ModuleGene.weight must be greater than zero."
            )

    @classmethod
    def from_dict(cls, data: KnowledgeRecord) -> "ModuleGene":
        """Create a module-gene definition from a mapping."""

        model_name = cls.__name__

        return cls(
            gene_id=_required_string(data, "gene_id", model_name),
            requirement=_enum_value(
                data,
                "requirement",
                GeneRequirement,
                model_name,
            ),
            weight=_positive_float(
                data,
                "weight",
                model_name,
                default=1.0,
            ),
            taxonomic_scope=_string_tuple(
                data,
                "taxonomic_scope",
                model_name,
            ),
            notes=_optional_string(data, "notes", model_name),
        )

@dataclass(frozen=True)
class ReferenceGene:
    """Curated definition of a reference gene."""

    gene_id: str
    symbol: str
    full_name: str
    aliases: tuple[str, ...] = ()
    description: str | None = None

    def __post_init__(self) -> None:
        """Validate the reference-gene definition."""

        if not self.gene_id.strip():
            raise KnowledgeModelError(
                "ReferenceGene.gene_id must not be empty."
            )

        if not self.symbol.strip():
            raise KnowledgeModelError(
                "ReferenceGene.symbol must not be empty."
            )

        if not self.full_name.strip():
            raise KnowledgeModelError(
                "ReferenceGene.full_name must not be empty."
            )

    @classmethod
    def from_dict(
        cls,
        data: KnowledgeRecord,
    ) -> "ReferenceGene":
        """Create a reference gene from a mapping."""

        model_name = cls.__name__

        return cls(
            gene_id=_required_string(
                data,
                "gene_id",
                model_name,
            ),
            symbol=_required_string(
                data,
                "symbol",
                model_name,
            ),
            full_name=_required_string(
                data,
                "full_name",
                model_name,
            ),
            aliases=_string_tuple(
                data,
                "aliases",
                model_name,
            ),
            description=_optional_string(
                data,
                "description",
                model_name,
            ),
        )

@dataclass(frozen=True)
class ReferenceModule:
    """Curated definition of a biological module."""

    module_id: str
    name: str
    pathway_id: str
    description: str
    genes: tuple[ModuleGene, ...]
    evidence_level: KnowledgeEvidenceLevel
    literature: tuple[str, ...] = ()
    taxonomic_scope: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        """Validate the reference-module definition."""

        if not self.module_id.strip():
            raise KnowledgeModelError(
                "ReferenceModule.module_id must not be empty."
            )

        if not self.name.strip():
            raise KnowledgeModelError(
                "ReferenceModule.name must not be empty."
            )

        if not self.pathway_id.strip():
            raise KnowledgeModelError(
                "ReferenceModule.pathway_id must not be empty."
            )

        if not self.genes:
            raise KnowledgeModelError(
                "ReferenceModule must contain at least one gene."
            )

        gene_ids = [gene.gene_id for gene in self.genes]

        if len(gene_ids) != len(set(gene_ids)):
            raise KnowledgeModelError(
                "ReferenceModule contains duplicate gene identifiers."
            )

    @classmethod
    def from_dict(cls, data: KnowledgeRecord) -> "ReferenceModule":
        """Create a reference module from a mapping."""

        model_name = cls.__name__
        gene_records = _record_sequence(data, "genes", model_name)

        return cls(
            module_id=_required_string(data, "module_id", model_name),
            name=_required_string(data, "name", model_name),
            pathway_id=_required_string(
                data,
                "pathway_id",
                model_name,
            ),
            description=_required_string(
                data,
                "description",
                model_name,
            ),
            genes=tuple(
                ModuleGene.from_dict(record)
                for record in gene_records
            ),
            evidence_level=_enum_value(
                data,
                "evidence_level",
                KnowledgeEvidenceLevel,
                model_name,
            ),
            literature=_string_tuple(
                data,
                "literature",
                model_name,
            ),
            taxonomic_scope=_string_tuple(
                data,
                "taxonomic_scope",
                model_name,
            ),
        )

    @property
    def required_gene_ids(self) -> tuple[str, ...]:
        """Return identifiers of genes required by the module."""

        return tuple(
            gene.gene_id
            for gene in self.genes
            if gene.requirement is GeneRequirement.REQUIRED
        )

    @property
    def supporting_gene_ids(self) -> tuple[str, ...]:
        """Return identifiers of genes supporting the module."""

        return tuple(
            gene.gene_id
            for gene in self.genes
            if gene.requirement is GeneRequirement.SUPPORTING
        )

    @property
    def optional_gene_ids(self) -> tuple[str, ...]:
        """Return identifiers of optional module genes."""

        return tuple(
            gene.gene_id
            for gene in self.genes
            if gene.requirement is GeneRequirement.OPTIONAL
        )

    @property
    def total_weight(self) -> float:
        """Return the total weight of all genes in the module."""

        return sum(gene.weight for gene in self.genes)


@dataclass(frozen=True)
class ReferencePathway:
    """Curated definition of a biological pathway."""

    pathway_id: str
    name: str
    description: str
    module_ids: tuple[str, ...]
    evidence_level: KnowledgeEvidenceLevel
    literature: tuple[str, ...] = ()
    taxonomic_scope: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        """Validate the reference-pathway definition."""

        if not self.pathway_id.strip():
            raise KnowledgeModelError(
                "ReferencePathway.pathway_id must not be empty."
            )

        if not self.name.strip():
            raise KnowledgeModelError(
                "ReferencePathway.name must not be empty."
            )

        if not self.module_ids:
            raise KnowledgeModelError(
                "ReferencePathway must contain at least one module."
            )

        if len(self.module_ids) != len(set(self.module_ids)):
            raise KnowledgeModelError(
                "ReferencePathway contains duplicate module identifiers."
            )

    @classmethod
    def from_dict(cls, data: KnowledgeRecord) -> "ReferencePathway":
        """Create a reference pathway from a mapping."""

        model_name = cls.__name__

        return cls(
            pathway_id=_required_string(
                data,
                "pathway_id",
                model_name,
            ),
            name=_required_string(data, "name", model_name),
            description=_required_string(
                data,
                "description",
                model_name,
            ),
            module_ids=_string_tuple(
                data,
                "module_ids",
                model_name,
            ),
            evidence_level=_enum_value(
                data,
                "evidence_level",
                KnowledgeEvidenceLevel,
                model_name,
            ),
            literature=_string_tuple(
                data,
                "literature",
                model_name,
            ),
            taxonomic_scope=_string_tuple(
                data,
                "taxonomic_scope",
                model_name,
            ),
        )


@dataclass(frozen=True)
class ReferencePhenotype:
    """Curated phenotype inferred from one or more pathways."""

    phenotype_id: str
    name: str
    description: str
    pathway_ids: tuple[str, ...]
    evidence_level: KnowledgeEvidenceLevel
    literature: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        """Validate the reference-phenotype definition."""

        if not self.phenotype_id.strip():
            raise KnowledgeModelError(
                "ReferencePhenotype.phenotype_id must not be empty."
            )

        if not self.name.strip():
            raise KnowledgeModelError(
                "ReferencePhenotype.name must not be empty."
            )

        if not self.pathway_ids:
            raise KnowledgeModelError(
                "ReferencePhenotype must contain at least one pathway."
            )

        if len(self.pathway_ids) != len(set(self.pathway_ids)):
            raise KnowledgeModelError(
                "ReferencePhenotype contains duplicate pathway identifiers."
            )

    @classmethod
    def from_dict(cls, data: KnowledgeRecord) -> "ReferencePhenotype":
        """Create a reference phenotype from a mapping."""

        model_name = cls.__name__

        return cls(
            phenotype_id=_required_string(
                data,
                "phenotype_id",
                model_name,
            ),
            name=_required_string(data, "name", model_name),
            description=_required_string(
                data,
                "description",
                model_name,
            ),
            pathway_ids=_string_tuple(
                data,
                "pathway_ids",
                model_name,
            ),
            evidence_level=_enum_value(
                data,
                "evidence_level",
                KnowledgeEvidenceLevel,
                model_name,
            ),
            literature=_string_tuple(
                data,
                "literature",
                model_name,
            ),
        )
