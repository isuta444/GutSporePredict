"""Tests for phylogenetic genome-manifest construction."""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from gutsporepredict.phylogeny.manifest import (
    GenomeManifestBuilder,
)


def write_metadata(path: Path) -> None:
    """Write minimal GTDB-style metadata for testing."""

    columns = [
        "gtdb_accession",
        "ncbi_accession",
        "domain",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species",
    ]

    rows = [
        {
            "gtdb_accession": "GB_GCA_000001.1",
            "ncbi_accession": "GCA_000001.1",
            "domain": "Bacteria",
            "phylum": "Bacillota",
            "class": "Clostridia",
            "order": "Oscillospirales",
            "family": "Acutalibacteraceae",
            "genus": "Examplegenus",
            "species": "Examplegenus species1",
        },
        {
            "gtdb_accession": "GB_GCA_000002.1",
            "ncbi_accession": "GCA_000002.1",
            "domain": "Bacteria",
            "phylum": "Bacillota",
            "class": "Clostridia",
            "order": "Lachnospirales",
            "family": "Lachnospiraceae",
            "genus": "Othergenus",
            "species": "Othergenus species2",
        },
    ]

    with path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)


def test_build_manifest_links_metadata_to_genomes(
    tmp_path: Path,
) -> None:
    """Only metadata rows with matching FASTA files are written."""

    metadata_path = tmp_path / "metadata.tsv"
    genome_dir = tmp_path / "genomes"
    output_path = tmp_path / "output" / "genome_manifest.tsv"

    genome_dir.mkdir()
    write_metadata(metadata_path)

    genome_path = genome_dir / "GCA_000001.1.fna"
    genome_path.write_text(
        ">contig1\nATGCGT\n",
        encoding="utf-8",
    )

    rows = GenomeManifestBuilder().build(
        metadata_path=metadata_path,
        genome_dir=genome_dir,
        output_path=output_path,
    )

    assert len(rows) == 1
    assert rows[0].genome_id == "GCA_000001.1"
    assert rows[0].family == "Acutalibacteraceae"
    assert rows[0].genome_path == str(genome_path.resolve())

    with output_path.open(
        encoding="utf-8",
        newline="",
    ) as handle:
        output_rows = list(
            csv.DictReader(
                handle,
                delimiter="\t",
            )
        )

    assert len(output_rows) == 1
    assert output_rows[0]["class"] == "Clostridia"
    assert output_rows[0]["species"] == "Examplegenus species1"


def test_build_manifest_rejects_missing_metadata(
    tmp_path: Path,
) -> None:
    """A missing metadata file raises an informative error."""

    genome_dir = tmp_path / "genomes"
    genome_dir.mkdir()

    with pytest.raises(
        FileNotFoundError,
        match="Metadata TSV was not found",
    ):
        GenomeManifestBuilder().build(
            metadata_path=tmp_path / "missing.tsv",
            genome_dir=genome_dir,
            output_path=tmp_path / "manifest.tsv",
        )


def test_build_manifest_rejects_no_matching_genomes(
    tmp_path: Path,
) -> None:
    """Metadata and FASTA sets must share at least one accession."""

    metadata_path = tmp_path / "metadata.tsv"
    genome_dir = tmp_path / "genomes"

    genome_dir.mkdir()
    write_metadata(metadata_path)

    (genome_dir / "GCA_999999.1.fna").write_text(
        ">contig1\nATGCGT\n",
        encoding="utf-8",
    )

    with pytest.raises(
        ValueError,
        match="No metadata accessions matched",
    ):
        GenomeManifestBuilder().build(
            metadata_path=metadata_path,
            genome_dir=genome_dir,
            output_path=tmp_path / "manifest.tsv",
        )
