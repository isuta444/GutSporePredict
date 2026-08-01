#!/usr/bin/env python3
"""Run HMMER searches for multiple genes across multiple protein FASTA files."""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
from pathlib import Path


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(
        description=(
            "Run one HMM profile against every protein FASTA file "
            "for all genes listed in a target TSV."
        )
    )
    parser.add_argument(
        "--targets",
        type=Path,
        required=True,
        help="TSV containing a gene_id column.",
    )
    parser.add_argument(
        "--hmm-dir",
        type=Path,
        required=True,
        help="Directory containing <gene_id>.hmm files.",
    )
    parser.add_argument(
        "--protein-dir",
        type=Path,
        required=True,
        help="Directory containing protein FASTA files.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Output directory with one subdirectory per genome.",
    )
    parser.add_argument(
        "--cpu",
        type=int,
        default=1,
        help="CPU count passed to each hmmsearch process.",
    )
    parser.add_argument(
        "--evalue",
        type=float,
        default=1e-5,
        help="Full-sequence E-value cutoff.",
    )
    parser.add_argument(
        "--domain-evalue",
        type=float,
        default=1e-4,
        help="Domain E-value cutoff.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing HMMER output files.",
    )
    return parser.parse_args()


def load_gene_ids(path: Path) -> list[str]:
    """Load unique gene IDs while retaining input order."""

    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        if reader.fieldnames is None or "gene_id" not in reader.fieldnames:
            raise ValueError(
                f"Target TSV must contain a gene_id column: {path}"
            )

        gene_ids: list[str] = []
        seen: set[str] = set()

        for row in reader:
            gene_id = row["gene_id"].strip()

            if not gene_id or gene_id in seen:
                continue

            seen.add(gene_id)
            gene_ids.append(gene_id)

    if not gene_ids:
        raise ValueError(f"No gene IDs found in: {path}")

    return gene_ids


def find_protein_fastas(directory: Path) -> list[Path]:
    """Return protein FASTA files from the input directory."""

    if not directory.is_dir():
        raise NotADirectoryError(
            f"Protein directory does not exist: {directory}"
        )

    fastas = sorted(
        path
        for path in directory.iterdir()
        if path.is_file() and path.suffix.lower() == ".faa"
    )

    if not fastas:
        raise ValueError(
            f"No .faa files found in protein directory: {directory}"
        )

    return fastas


def run_hmmsearch(
    executable: str,
    hmm_path: Path,
    protein_fasta: Path,
    output_dir: Path,
    *,
    cpu: int,
    evalue: float,
    domain_evalue: float,
    force: bool,
) -> None:
    """Run one HMM profile against one protein FASTA."""

    gene_id = hmm_path.stem
    tblout = output_dir / f"{gene_id}.tblout"
    domtblout = output_dir / f"{gene_id}.domtblout"
    stdout_path = output_dir / f"{gene_id}.stdout.txt"
    stderr_path = output_dir / f"{gene_id}.stderr.txt"

    outputs = (tblout, domtblout, stdout_path, stderr_path)

    if not force and tblout.exists() and domtblout.exists():
        print(
            f"Skipping existing: {gene_id} vs {protein_fasta.stem}"
        )
        return

    output_dir.mkdir(parents=True, exist_ok=True)

    command = [
        executable,
        "--cpu",
        str(cpu),
        "--noali",
        "-E",
        str(evalue),
        "--domE",
        str(domain_evalue),
        "--tblout",
        str(tblout),
        "--domtblout",
        str(domtblout),
        str(hmm_path),
        str(protein_fasta),
    ]

    print(f"Searching {gene_id} in {protein_fasta.stem}")

    with stdout_path.open("w", encoding="utf-8") as stdout_handle, (
        stderr_path.open("w", encoding="utf-8")
    ) as stderr_handle:
        completed = subprocess.run(
            command,
            stdout=stdout_handle,
            stderr=stderr_handle,
            text=True,
            check=False,
        )

    if completed.returncode != 0:
        for path in outputs:
            if path.exists() and path.stat().st_size == 0:
                path.unlink()

        raise RuntimeError(
            "hmmsearch failed for "
            f"{gene_id} against {protein_fasta}. "
            f"See: {stderr_path}"
        )


def main() -> None:
    """Run all gene-by-genome HMM searches."""

    args = parse_args()

    if args.cpu < 1:
        raise ValueError("--cpu must be at least 1")

    executable = shutil.which("hmmsearch")

    if executable is None:
        raise RuntimeError(
            "hmmsearch was not found in the active environment."
        )

    gene_ids = load_gene_ids(args.targets)
    protein_fastas = find_protein_fastas(args.protein_dir)

    missing_hmms = [
        args.hmm_dir / f"{gene_id}.hmm"
        for gene_id in gene_ids
        if not (args.hmm_dir / f"{gene_id}.hmm").is_file()
    ]

    if missing_hmms:
        missing_text = "\n".join(str(path) for path in missing_hmms)
        raise FileNotFoundError(
            f"Missing HMM profiles:\n{missing_text}"
        )

    comparisons = 0

    for protein_fasta in protein_fastas:
        genome_output_dir = args.output_dir / protein_fasta.stem

        for gene_id in gene_ids:
            run_hmmsearch(
                executable,
                args.hmm_dir / f"{gene_id}.hmm",
                protein_fasta,
                genome_output_dir,
                cpu=args.cpu,
                evalue=args.evalue,
                domain_evalue=args.domain_evalue,
                force=args.force,
            )
            comparisons += 1

    print()
    print(f"Genes: {len(gene_ids)}")
    print(f"Genomes: {len(protein_fastas)}")
    print(f"Comparisons: {comparisons}")
    print(f"Output: {args.output_dir}")


if __name__ == "__main__":
    main()
