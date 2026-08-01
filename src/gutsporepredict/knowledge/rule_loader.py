"""Load biological rules from YAML files."""

from pathlib import Path

import yaml

from gutsporepredict.knowledge.models import KnowledgeEvidenceLevel
from gutsporepredict.knowledge.rules import ReferenceRule


def load_rule(path: str | Path) -> ReferenceRule:
    """Load one rule from a YAML file."""

    with open(path, encoding="utf-8") as fh:
        data = yaml.safe_load(fh)

    return ReferenceRule(
        rule_id=data["rule_id"],
        module_id=data["module_id"],
        condition=data["condition"],
        outcome=data["outcome"],
        evidence_level=KnowledgeEvidenceLevel(
            data["evidence_level"]
        ),
        literature=tuple(data.get("literature", ())),
        notes=data.get("notes"),
    )
