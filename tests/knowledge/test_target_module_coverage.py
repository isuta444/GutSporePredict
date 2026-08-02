"""Validate assignment of searched genes to knowledge modules."""

import csv
from collections import Counter
from pathlib import Path

import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[2]


def read_target_genes(path: Path) -> set[str]:
    """Read gene identifiers from a target TSV file."""
    with path.open(
        encoding="utf-8",
        newline="",
    ) as handle:
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )
        return {
            row["gene_id"].strip()
            for row in reader
            if row["gene_id"].strip()
        }


def test_all_target_genes_have_one_module() -> None:
    """Every searched gene should belong to exactly one module."""
    target_genes = (
        read_target_genes(
            PROJECT_ROOT
            / "config"
            / "gtdb_targets"
            / "02_sporulation_ready.tsv"
        )
        | read_target_genes(
            PROJECT_ROOT
            / "config"
            / "gtdb_targets"
            / "03_germination_ready.tsv"
        )
    )

    assignments: list[str] = []

    for path in sorted(
        (
            PROJECT_ROOT
            / "knowledge"
            / "modules"
        ).glob("*.yaml")
    ):
        raw = yaml.safe_load(
            path.read_text(encoding="utf-8")
        )

        assignments.extend(
            gene["gene_id"]
            for gene in raw["genes"]
        )

    counts = Counter(assignments)

    assert len(target_genes) == 45
    assert set(counts) == target_genes
    assert all(count == 1 for count in counts.values())
