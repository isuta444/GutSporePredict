"""Load reference modules from YAML files."""

from pathlib import Path

import yaml

from gutsporepredict.knowledge.models import ReferenceModule


def load_module(path: str | Path) -> ReferenceModule:
    """Load one reference module from YAML."""

    with open(path, encoding="utf-8") as fh:
        data = yaml.safe_load(fh)

    return ReferenceModule.from_dict(data)
