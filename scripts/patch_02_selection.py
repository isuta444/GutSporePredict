from pathlib import Path

path = Path("scripts/02_select_genomes.py")
text = path.read_text(encoding="utf-8")


old_type = '''def detect_type_material(row: dict[str, str]) -> tuple[bool, str]:
    value, column = first_existing_value(
        row,
        [
            "gtdb_type_designation",
            "gtdb_type_designation_sources",
            "ncbi_type_material_designation",
            "ncbi_type_material",
        ],
    )

    value_lower = value.lower()

    negative_terms = {
        "",
        "not type material",
        "none",
        "na",
    }

    if value_lower in negative_terms:
        return False, ""

    positive_terms = [
        "type strain",
        "type material",
        "type species",
        "type subspecies",
        "neotype",
        "reference strain",
    ]

    is_type = any(term in value_lower for term in positive_terms)

    return is_type, f"{column}:{value}" if value else ""
'''

new_type = '''def detect_type_material(row: dict[str, str]) -> tuple[bool, str]:
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
'''

old_category = '''def detect_genome_category(row: dict[str, str]) -> str:
    value, _ = first_existing_value(
        row,
        [
            "ncbi_genome_representation",
            "ncbi_genome_category",
            "genome_category",
        ],
    )

    value_lower = value.lower()

    if "metagenome" in value_lower or value_lower == "mag":
        return "MAG"

    if "single cell" in value_lower or value_lower == "sag":
        return "SAG"

    if "isolate" in value_lower:
        return "isolate"

    return value or "unknown"
'''

new_category = '''def detect_genome_category(row: dict[str, str]) -> str:
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
'''

old_human = '''def detect_human_gut_score(
    row: dict[str, str],
    keywords: list[str],
) -> tuple[int, str]:
    candidate_columns = [
        "ncbi_isolation_source",
        "ncbi_host",
        "ncbi_environment",
        "ncbi_env_broad_scale",
        "ncbi_env_local_scale",
        "ncbi_env_medium",
        "ncbi_biome",
    ]

    text_parts = []

    for column in candidate_columns:
        value = clean_text(row.get(column))

        if value:
            text_parts.append(value)

    combined = " | ".join(text_parts)
    lower = combined.lower()

    matched = sorted({
        keyword
        for keyword in keywords
        if keyword.lower() in lower
    })

    if not matched:
        return 0, combined

    strong_terms = {
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
    }

    strong_match = any(term in lower for term in strong_terms)
    human_match = (
        "human" in lower
        or "homo sapiens" in lower
    )

    if strong_match and human_match:
        return 3, combined

    if strong_match:
        return 2, combined

    return 1, combined
'''

new_human = '''def classify_source_origin(
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
'''

replacements = [
    (old_type, new_type),
    (old_category, new_category),
    (old_human, new_human),
]

for old, new in replacements:
    if old not in text:
        raise SystemExit(
            "Expected code block was not found. "
            "The script may already have been modified."
        )

    text = text.replace(old, new)


old_call = '''            human_gut_score, source_text = (
                detect_human_gut_score(
                    row,
                    human_keywords,
                )
            )
'''

new_call = '''            (
                host_origin,
                body_site,
                human_gut_score,
                source_text,
            ) = classify_source_origin(row)
'''

if old_call not in text:
    raise SystemExit("Source-classification call was not found.")

text = text.replace(old_call, new_call)


old_record = '''                "human_gut_score": human_gut_score,
                "source_metadata": source_text,
'''

new_record = '''                "host_origin": host_origin,
                "body_site": body_site,
                "human_gut_score": human_gut_score,
                "source_metadata": source_text,
                "genome_representation": clean_text(
                    row.get("ncbi_genome_representation")
                ),
'''

if old_record not in text:
    raise SystemExit("Output record block was not found.")

text = text.replace(old_record, new_record)


old_fields = '''        "human_gut_score",
        "source_metadata",
        "ncbi_bioproject",
'''

new_fields = '''        "host_origin",
        "body_site",
        "human_gut_score",
        "source_metadata",
        "genome_representation",
        "ncbi_bioproject",
'''

if old_fields not in text:
    raise SystemExit("Output fields block was not found.")

text = text.replace(old_fields, new_fields)

path.write_text(text, encoding="utf-8")

print(f"Patched: {path}")
