#!/usr/bin/env python3

"""
GutSporePredict v4.0-alpha1
Genome selection and quality ranking from GTDB bacterial metadata.

Outputs:
    metadata/all_target_genomes.tsv
    metadata/primary_genome_set.tsv
    metadata/genome_selection_summary.tsv
    metadata/genome_selection_summary.json
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Optional


RANK_ORDER = {
    "A": 1,
    "B": 2,
    "C": 3,
    "D": 4,
    "E": 5,
}


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Select and rank target GTDB genomes for "
            "GutSporePredict comparative genomics."
        )
    )

    parser.add_argument(
        "--metadata",
        default="database/gtdb/bac120_metadata.tsv.gz",
        help="GTDB bacterial metadata file.",
    )

    parser.add_argument(
        "--config",
        default="config/genome_selection.json",
        help="Genome-selection configuration JSON.",
    )

    parser.add_argument(
        "--output-dir",
        default="metadata",
        help="Output directory.",
    )

    return parser.parse_args()


def clean_text(value: Optional[str]) -> str:
    if value is None:
        return ""

    value = str(value).strip()

    if value.lower() in {
        "",
        "na",
        "nan",
        "none",
        "null",
        "not available",
    }:
        return ""

    return value


def parse_float(value: Optional[str]) -> Optional[float]:
    value = clean_text(value)

    if not value:
        return None

    try:
        result = float(value)
    except ValueError:
        return None

    if not math.isfinite(result):
        return None

    return result


def first_existing_value(
    row: dict[str, str],
    candidate_columns: list[str],
) -> tuple[str, str]:
    for column in candidate_columns:
        if column in row:
            value = clean_text(row.get(column))

            if value:
                return value, column

    return "", ""


def truthy(value: Optional[str]) -> bool:
    return clean_text(value).lower() in {
        "t",
        "true",
        "yes",
        "y",
        "1",
    }


def parse_taxonomy(taxonomy: str) -> dict[str, str]:
    rank_prefixes = {
        "d__": "domain",
        "p__": "phylum",
        "c__": "class",
        "o__": "order",
        "f__": "family",
        "g__": "genus",
        "s__": "species",
    }

    result = {
        "domain": "",
        "phylum": "",
        "class": "",
        "order": "",
        "family": "",
        "genus": "",
        "species": "",
    }

    for item in clean_text(taxonomy).split(";"):
        item = item.strip()

        for prefix, rank in rank_prefixes.items():
            if item.startswith(prefix):
                result[rank] = item[len(prefix):].strip()
                break

    return result


def remove_gtdb_prefix(accession: str) -> str:
    accession = clean_text(accession)

    for prefix in ("RS_", "GB_"):
        if accession.startswith(prefix):
            return accession[len(prefix):]

    return accession


def detect_representative(row: dict[str, str]) -> bool:
    value, _ = first_existing_value(
        row,
        [
            "gtdb_representative",
            "species_representative",
        ],
    )

    return truthy(value)


def detect_type_material(row: dict[str, str]) -> tuple[bool, str]:
    evidence = []

    gtdb_value = clean_text(
        row.get("gtdb_type_designation_ncbi_taxa")
    )

    gtdb_source = clean_text(
        row.get("gtdb_type_designation_ncbi_taxa_sources")
    )

    ncbi_value = clean_text(
        row.get("ncbi_type_material_designation")
    )

    positive_terms = [
        "type strain",
        "type material",
        "type species",
        "type subspecies",
        "heterotypic synonym",
        "neotype",
    ]

    is_type = False

    if any(term in gtdb_value.lower() for term in positive_terms):
        is_type = True
        evidence.append(
            f"gtdb_type_designation_ncbi_taxa:{gtdb_value}"
        )

        if gtdb_source:
            evidence.append(
                f"gtdb_source:{gtdb_source}"
            )

    if any(term in ncbi_value.lower() for term in positive_terms):
        is_type = True
        evidence.append(
            f"ncbi_type_material_designation:{ncbi_value}"
        )

    return is_type, " | ".join(evidence)


def detect_assembly_level(row: dict[str, str]) -> str:
    value, _ = first_existing_value(
        row,
        [
            "ncbi_assembly_level",
            "assembly_level",
        ],
    )

    return value


def assembly_is_high_quality(assembly_level: str) -> bool:
    value = clean_text(assembly_level).lower()

    return value in {
        "complete genome",
        "chromosome",
        "complete",
    }


def detect_genome_category(row: dict[str, str]) -> str:
    category = clean_text(
        row.get("ncbi_genome_category")
    ).lower()

    isolate_name = clean_text(
        row.get("ncbi_isolate")
    )

    if "derived from metagenome" in category:
        return "MAG"

    if "derived from single cell" in category:
        return "SAG"

    if isolate_name:
        return "isolate"

    if category in {"", "none", "na"}:
        return "isolate_or_unspecified"

    return category


def classify_source_origin(
    row: dict[str, str],
) -> tuple[str, str, int, str]:
    source = clean_text(
        row.get("ncbi_isolation_source")
    )

    isolate = clean_text(
        row.get("ncbi_isolate")
    )

    combined = " | ".join(
        value for value in [source, isolate] if value
    )

    lower = combined.lower()

    gut_terms = [
        "feces",
        "faeces",
        "fecal",
        "faecal",
        "stool",
        "gut",
        "intestinal",
        "intestine",
        "colon",
        "colonic",
        "gastrointestinal",
        "rectal",
    ]

    human_terms = [
        "human",
        "homo sapiens",
        "patient",
    ]

    nonhuman_terms = [
        "mouse",
        "mice",
        "murine",
        "rat",
        "bovine",
        "cow",
        "cattle",
        "pig",
        "swine",
        "chicken",
        "avian",
        "dog",
        "canine",
        "cat",
        "feline",
        "horse",
        "equine",
        "sheep",
        "goat",
        "termite",
    ]

    is_gut = any(term in lower for term in gut_terms)
    is_human = any(term in lower for term in human_terms)
    is_nonhuman = any(term in lower for term in nonhuman_terms)

    if is_human:
        host_origin = "human"
    elif is_nonhuman:
        host_origin = "nonhuman"
    else:
        host_origin = "unknown"

    if is_gut:
        body_site = "gut"
    elif combined:
        body_site = "non_gut_or_unspecified"
    else:
        body_site = "unknown"

    if host_origin == "human" and body_site == "gut":
        human_gut_score = 3
    elif body_site == "gut":
        human_gut_score = 2
    elif host_origin == "human":
        human_gut_score = 1
    else:
        human_gut_score = 0

    return (
        host_origin,
        body_site,
        human_gut_score,
        combined,
    )


def determine_quality_rank(
    representative: bool,
    type_material: bool,
    high_quality_assembly: bool,
    completeness: Optional[float],
    contamination: Optional[float],
    minimum_completeness: float,
    maximum_contamination: float,
) -> str:
    passes_quality = (
        completeness is not None
        and contamination is not None
        and completeness >= minimum_completeness
        and contamination <= maximum_contamination
    )

    if not passes_quality:
        return "E"

    if representative and type_material and high_quality_assembly:
        return "A"

    if representative and high_quality_assembly:
        return "B"

    if representative:
        return "C"

    return "D"


def calculate_priority_score(
    rank: str,
    representative: bool,
    type_material: bool,
    high_quality_assembly: bool,
    genome_category: str,
    human_gut_score: int,
    completeness: Optional[float],
    contamination: Optional[float],
) -> float:
    score = {
        "A": 100.0,
        "B": 85.0,
        "C": 70.0,
        "D": 50.0,
        "E": 0.0,
    }[rank]

    if representative:
        score += 10.0

    if type_material:
        score += 8.0

    if high_quality_assembly:
        score += 6.0

    if genome_category == "isolate":
        score += 4.0
    elif genome_category == "MAG":
        score -= 2.0
    elif genome_category == "SAG":
        score -= 3.0

    score += human_gut_score * 2.0

    if completeness is not None:
        score += min(completeness, 100.0) / 20.0

    if contamination is not None:
        score -= contamination

    return round(score, 3)


def choose_completeness_and_contamination(
    row: dict[str, str],
    prefer_checkm2: bool,
) -> tuple[
    Optional[float],
    Optional[float],
    str,
    str,
]:
    if prefer_checkm2:
        completeness_columns = [
            "checkm2_completeness",
            "checkm_completeness",
        ]

        contamination_columns = [
            "checkm2_contamination",
            "checkm_contamination",
        ]
    else:
        completeness_columns = [
            "checkm_completeness",
            "checkm2_completeness",
        ]

        contamination_columns = [
            "checkm_contamination",
            "checkm2_contamination",
        ]

    completeness_raw, completeness_source = first_existing_value(
        row,
        completeness_columns,
    )

    contamination_raw, contamination_source = first_existing_value(
        row,
        contamination_columns,
    )

    return (
        parse_float(completeness_raw),
        parse_float(contamination_raw),
        completeness_source,
        contamination_source,
    )


def write_tsv(
    path: Path,
    records: list[dict[str, object]],
    fieldnames: list[str],
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=fieldnames,
            extrasaction="ignore",
        )

        writer.writeheader()

        for record in records:
            writer.writerow(record)


def main() -> int:
    args = parse_arguments()

    metadata_path = Path(args.metadata)
    config_path = Path(args.config)
    output_dir = Path(args.output_dir)

    if not metadata_path.exists():
        print(
            f"[ERROR] GTDB metadata not found: {metadata_path}",
            file=sys.stderr,
        )
        return 1

    if not config_path.exists():
        print(
            f"[ERROR] Configuration not found: {config_path}",
            file=sys.stderr,
        )
        return 1

    config = json.loads(
        config_path.read_text(encoding="utf-8")
    )

    target_families = set(config["target_families"])

    thresholds = config["quality_thresholds"]

    minimum_completeness = float(
        thresholds["minimum_completeness"]
    )

    maximum_contamination = float(
        thresholds["maximum_contamination"]
    )

    strict_minimum_completeness = float(
        thresholds["strict_minimum_completeness"]
    )

    strict_maximum_contamination = float(
        thresholds["strict_maximum_contamination"]
    )

    selection_config = config["selection"]

    prefer_checkm2 = bool(
        selection_config.get("prefer_checkm2", True)
    )

    representative_primary = bool(
        selection_config.get(
            "representatives_only_for_primary_set",
            True,
        )
    )

    human_keywords = [
        str(keyword).lower()
        for keyword in config.get(
            "human_gut_keywords",
            [],
        )
    ]

    output_dir.mkdir(parents=True, exist_ok=True)

    all_records: list[dict[str, object]] = []
    family_counts = Counter()
    rank_counts = Counter()
    category_counts = Counter()
    primary_family_counts = Counter()
    source_counts = Counter()

    print("=" * 72)
    print("GutSporePredict v4.0-alpha1")
    print("Genome Selection Policy")
    print("=" * 72)
    print(f"Metadata: {metadata_path}")
    print(f"Target families: {len(target_families)}")
    print(
        f"Quality threshold: completeness >= "
        f"{minimum_completeness}, contamination <= "
        f"{maximum_contamination}"
    )
    print()

    with gzip.open(
        metadata_path,
        "rt",
        encoding="utf-8",
        errors="replace",
        newline="",
    ) as handle:
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        if not reader.fieldnames:
            print(
                "[ERROR] Metadata header could not be read.",
                file=sys.stderr,
            )
            return 1

        required_columns = {
            "accession",
            "gtdb_taxonomy",
        }

        missing = required_columns - set(reader.fieldnames)

        if missing:
            print(
                "[ERROR] Required columns missing: "
                + ", ".join(sorted(missing)),
                file=sys.stderr,
            )
            return 1

        for row in reader:
            taxonomy = parse_taxonomy(
                row.get("gtdb_taxonomy", "")
            )

            family = taxonomy["family"]

            if family not in target_families:
                continue

            accession = clean_text(
                row.get("accession")
            )

            ncbi_accession = remove_gtdb_prefix(
                accession
            )

            (
                completeness,
                contamination,
                completeness_source,
                contamination_source,
            ) = choose_completeness_and_contamination(
                row,
                prefer_checkm2,
            )

            representative = detect_representative(row)

            type_material, type_evidence = (
                detect_type_material(row)
            )

            assembly_level = detect_assembly_level(row)

            high_quality_assembly = (
                assembly_is_high_quality(
                    assembly_level
                )
            )

            genome_category = detect_genome_category(row)

            (
                host_origin,
                body_site,
                human_gut_score,
                source_text,
            ) = classify_source_origin(row)

            quality_rank = determine_quality_rank(
                representative=representative,
                type_material=type_material,
                high_quality_assembly=high_quality_assembly,
                completeness=completeness,
                contamination=contamination,
                minimum_completeness=minimum_completeness,
                maximum_contamination=maximum_contamination,
            )

            strict_quality = (
                completeness is not None
                and contamination is not None
                and completeness
                >= strict_minimum_completeness
                and contamination
                <= strict_maximum_contamination
            )

            primary_selected = (
                quality_rank != "E"
                and (
                    representative
                    if representative_primary
                    else True
                )
            )

            priority_score = calculate_priority_score(
                rank=quality_rank,
                representative=representative,
                type_material=type_material,
                high_quality_assembly=high_quality_assembly,
                genome_category=genome_category,
                human_gut_score=human_gut_score,
                completeness=completeness,
                contamination=contamination,
            )

            record = {
                "gtdb_accession": accession,
                "ncbi_accession": ncbi_accession,
                "domain": taxonomy["domain"],
                "phylum": taxonomy["phylum"],
                "class": taxonomy["class"],
                "order": taxonomy["order"],
                "family": taxonomy["family"],
                "genus": taxonomy["genus"],
                "species": taxonomy["species"],
                "quality_rank": quality_rank,
                "priority_score": priority_score,
                "primary_selected": str(
                    primary_selected
                ).lower(),
                "gtdb_representative": str(
                    representative
                ).lower(),
                "type_material": str(
                    type_material
                ).lower(),
                "type_evidence": type_evidence,
                "assembly_level": assembly_level,
                "high_quality_assembly": str(
                    high_quality_assembly
                ).lower(),
                "genome_category": genome_category,
                "completeness": (
                    ""
                    if completeness is None
                    else round(completeness, 3)
                ),
                "contamination": (
                    ""
                    if contamination is None
                    else round(contamination, 3)
                ),
                "strict_quality_95_3": str(
                    strict_quality
                ).lower(),
                "completeness_source": (
                    completeness_source
                ),
                "contamination_source": (
                    contamination_source
                ),
                "host_origin": host_origin,
                "body_site": body_site,
                "human_gut_score": human_gut_score,
                "source_metadata": source_text,
                "genome_representation": clean_text(
                    row.get("ncbi_genome_representation")
                ),
                "ncbi_bioproject": clean_text(
                    row.get("ncbi_bioproject")
                ),
                "ncbi_biosample": clean_text(
                    row.get("ncbi_biosample")
                ),
            }

            all_records.append(record)

            family_counts[family] += 1
            rank_counts[quality_rank] += 1
            category_counts[genome_category] += 1

            if primary_selected:
                primary_family_counts[family] += 1

            source_counts[
                completeness_source
                or "missing_completeness"
            ] += 1

    all_records.sort(
        key=lambda record: (
            record["family"],
            RANK_ORDER.get(
                str(record["quality_rank"]),
                99,
            ),
            -float(record["priority_score"]),
            str(record["species"]),
            str(record["gtdb_accession"]),
        )
    )

    primary_records = [
        record
        for record in all_records
        if record["primary_selected"] == "true"
    ]

    fieldnames = [
        "gtdb_accession",
        "ncbi_accession",
        "domain",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species",
        "quality_rank",
        "priority_score",
        "primary_selected",
        "gtdb_representative",
        "type_material",
        "type_evidence",
        "assembly_level",
        "high_quality_assembly",
        "genome_category",
        "completeness",
        "contamination",
        "strict_quality_95_3",
        "completeness_source",
        "contamination_source",
        "host_origin",
        "body_site",
        "human_gut_score",
        "source_metadata",
        "genome_representation",
        "ncbi_bioproject",
        "ncbi_biosample",
    ]

    all_path = output_dir / "all_target_genomes.tsv"

    primary_path = (
        output_dir / "primary_genome_set.tsv"
    )

    write_tsv(
        all_path,
        all_records,
        fieldnames,
    )

    write_tsv(
        primary_path,
        primary_records,
        fieldnames,
    )

    summary_rows = []

    for family in sorted(target_families):
        family_records = [
            record
            for record in all_records
            if record["family"] == family
        ]

        row = {
            "family": family,
            "all_target_genomes": len(
                family_records
            ),
            "primary_genomes": sum(
                record["primary_selected"]
                == "true"
                for record in family_records
            ),
            "rank_A": sum(
                record["quality_rank"] == "A"
                for record in family_records
            ),
            "rank_B": sum(
                record["quality_rank"] == "B"
                for record in family_records
            ),
            "rank_C": sum(
                record["quality_rank"] == "C"
                for record in family_records
            ),
            "rank_D": sum(
                record["quality_rank"] == "D"
                for record in family_records
            ),
            "rank_E": sum(
                record["quality_rank"] == "E"
                for record in family_records
            ),
            "human_gut_score_2_or_3": sum(
                int(record["human_gut_score"]) >= 2
                for record in family_records
            ),
        }

        summary_rows.append(row)

    summary_tsv_path = (
        output_dir
        / "genome_selection_summary.tsv"
    )

    summary_fields = [
        "family",
        "all_target_genomes",
        "primary_genomes",
        "rank_A",
        "rank_B",
        "rank_C",
        "rank_D",
        "rank_E",
        "human_gut_score_2_or_3",
    ]

    write_tsv(
        summary_tsv_path,
        summary_rows,
        summary_fields,
    )

    summary_json = {
        "pipeline": "GutSporePredict",
        "version": "4.0-alpha1",
        "policy": {
            "minimum_completeness": (
                minimum_completeness
            ),
            "maximum_contamination": (
                maximum_contamination
            ),
            "strict_minimum_completeness": (
                strict_minimum_completeness
            ),
            "strict_maximum_contamination": (
                strict_maximum_contamination
            ),
            "representatives_only_for_primary_set": (
                representative_primary
            ),
        },
        "total_target_genomes": len(all_records),
        "total_primary_genomes": len(
            primary_records
        ),
        "rank_counts": dict(
            sorted(rank_counts.items())
        ),
        "family_counts": dict(
            sorted(family_counts.items())
        ),
        "primary_family_counts": dict(
            sorted(primary_family_counts.items())
        ),
        "genome_category_counts": dict(
            sorted(category_counts.items())
        ),
        "quality_sources": dict(
            sorted(source_counts.items())
        ),
    }

    summary_json_path = (
        output_dir
        / "genome_selection_summary.json"
    )

    summary_json_path.write_text(
        json.dumps(
            summary_json,
            indent=2,
            ensure_ascii=False,
        )
        + "\n",
        encoding="utf-8",
    )

    print("=" * 72)
    print("[SUCCESS] Genome selection completed.")
    print(f"Target genomes:  {len(all_records):,}")
    print(f"Primary genomes: {len(primary_records):,}")
    print()

    for family in sorted(target_families):
        print(
            f"{family}: "
            f"{family_counts[family]:,} total, "
            f"{primary_family_counts[family]:,} primary"
        )

    print()
    print("Quality ranks:")

    for rank in ["A", "B", "C", "D", "E"]:
        print(
            f"  Rank {rank}: "
            f"{rank_counts[rank]:,}"
        )

    print()
    print(f"[OUTPUT] {all_path}")
    print(f"[OUTPUT] {primary_path}")
    print(f"[OUTPUT] {summary_tsv_path}")
    print(f"[OUTPUT] {summary_json_path}")
    print("=" * 72)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
