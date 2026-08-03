"""Build a genome manifest for phylogenetic analyses."""

from __future__ import annotations

import csv
from dataclasses import asdict, dataclass
from pathlib import Path

SUPPORTED_FASTA_SUFFIXES = (
    ".fna",
    ".fa",
    ".fasta",
    ".fas",
)


@dataclass(frozen=True)
class GenomeManifestRow:
    """One genome linked to its taxonomy and FASTA file."""

    genome_id: str
    ncbi_accession: str
    gtdb_accession: str
    domain: str
    phylum: str
    class_name: str
    order: str
    family: str
    genus: str
    species: str
    genome_path: str


class GenomeManifestBuilder:
    """Link analysis metadata with local genome FASTA files."""

    OUTPUT_COLUMNS = [
        "genome_id",
        "ncbi_accession",
        "gtdb_accession",
        "domain",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species",
        "genome_path",
    ]

    REQUIRED_METADATA_COLUMNS = {
        "gtdb_accession",
        "ncbi_accession",
        "domain",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species",
    }

    def build(
        self,
        metadata_path: Path,
        genome_dir: Path,
        output_path: Path,
    ) -> list[GenomeManifestRow]:
        """Build and write a genome manifest."""

        metadata_path = metadata_path.resolve()
        genome_dir = genome_dir.resolve()
        output_path = output_path.resolve()

        self._validate_inputs(
            metadata_path=metadata_path,
            genome_dir=genome_dir,
        )

        genome_files = self._index_genome_files(genome_dir)
        rows = self._load_rows(
            metadata_path=metadata_path,
            genome_files=genome_files,
        )

        output_path.parent.mkdir(
            parents=True,
            exist_ok=True,
        )
        self._write_rows(
            rows=rows,
            output_path=output_path,
        )

        return rows

    def _validate_inputs(
        self,
        metadata_path: Path,
        genome_dir: Path,
    ) -> None:
        """Validate metadata and genome input paths."""

        if not metadata_path.is_file():
            raise FileNotFoundError(f"Metadata TSV was not found: {metadata_path}")

        if not genome_dir.exists():
            raise FileNotFoundError(f"Genome directory was not found: {genome_dir}")

        if not genome_dir.is_dir():
            raise NotADirectoryError(f"Genome path is not a directory: {genome_dir}")

    def _index_genome_files(
        self,
        genome_dir: Path,
    ) -> dict[str, Path]:
        """Index genome FASTA files by accession-like filename stem."""

        indexed: dict[str, Path] = {}

        for path in sorted(genome_dir.iterdir()):
            if not path.is_file():
                continue

            if path.suffix.lower() not in SUPPORTED_FASTA_SUFFIXES:
                continue

            genome_id = path.stem

            if genome_id in indexed:
                raise ValueError(
                    f"Multiple genome FASTA files have the same stem: {genome_id}"
                )

            indexed[genome_id] = path.resolve()

        if not indexed:
            raise ValueError(f"No supported genome FASTA files found in: {genome_dir}")

        return indexed

    def _load_rows(
        self,
        metadata_path: Path,
        genome_files: dict[str, Path],
    ) -> list[GenomeManifestRow]:
        """Read metadata and retain genomes available locally."""

        rows: list[GenomeManifestRow] = []
        seen_accessions: set[str] = set()

        with metadata_path.open(
            encoding="utf-8",
            newline="",
        ) as handle:
            reader = csv.DictReader(
                handle,
                delimiter="\t",
            )

            fieldnames = set(reader.fieldnames or [])
            missing_columns = self.REQUIRED_METADATA_COLUMNS - fieldnames

            if missing_columns:
                raise ValueError(
                    "Metadata TSV is missing required columns: "
                    + ", ".join(sorted(missing_columns))
                )

            for metadata in reader:
                accession = metadata["ncbi_accession"].strip()

                if not accession:
                    continue

                genome_path = genome_files.get(accession)

                if genome_path is None:
                    continue

                if accession in seen_accessions:
                    raise ValueError(
                        f"Duplicate NCBI accession in metadata: {accession}"
                    )

                seen_accessions.add(accession)

                rows.append(
                    GenomeManifestRow(
                        genome_id=accession,
                        ncbi_accession=accession,
                        gtdb_accession=metadata["gtdb_accession"].strip(),
                        domain=metadata["domain"].strip(),
                        phylum=metadata["phylum"].strip(),
                        class_name=metadata["class"].strip(),
                        order=metadata["order"].strip(),
                        family=metadata["family"].strip(),
                        genus=metadata["genus"].strip(),
                        species=metadata["species"].strip(),
                        genome_path=str(genome_path),
                    )
                )

        if not rows:
            raise ValueError("No metadata accessions matched local genome FASTA files.")

        return rows

    def _write_rows(
        self,
        rows: list[GenomeManifestRow],
        output_path: Path,
    ) -> None:
        """Write manifest rows as TSV."""

        with output_path.open(
            "w",
            encoding="utf-8",
            newline="",
        ) as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=self.OUTPUT_COLUMNS,
                delimiter="\t",
            )
            writer.writeheader()

            for row in rows:
                values = asdict(row)
                values["class"] = values.pop("class_name")
                writer.writerow(values)
