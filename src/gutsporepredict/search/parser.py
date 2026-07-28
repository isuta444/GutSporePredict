"""Parsers for tabular sequence-search output."""

from pathlib import Path

from gutsporepredict.search.exceptions import SearchOutputError
from gutsporepredict.search.models import SearchHit

DIAMOND_COLUMN_COUNT = 12


def parse_diamond_output(
    output_file: str | Path,
) -> list[SearchHit]:
    """Parse GutSporePredict DIAMOND tabular output."""

    output_file = Path(output_file)

    if not output_file.exists():
        raise SearchOutputError(
            f"DIAMOND output file does not exist: {output_file}"
        )

    if not output_file.is_file():
        raise SearchOutputError(
            f"DIAMOND output path is not a file: {output_file}"
        )

    hits: list[SearchHit] = []

    with output_file.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()

            if not line or line.startswith("#"):
                continue

            columns = line.split("\t")

            if len(columns) != DIAMOND_COLUMN_COUNT:
                raise SearchOutputError(
                    "Invalid DIAMOND output at "
                    f"{output_file}:{line_number}. "
                    f"Expected {DIAMOND_COLUMN_COUNT} columns, "
                    f"found {len(columns)}."
                )

            (
                query_id,
                target_id,
                identity_text,
                alignment_length_text,
                query_length_text,
                target_length_text,
                _query_start,
                _query_end,
                _target_start,
                _target_end,
                evalue_text,
                bitscore_text,
            ) = columns

            try:
                identity = float(identity_text)
                alignment_length = int(alignment_length_text)
                query_length = int(query_length_text)
                target_length = int(target_length_text)
                evalue = float(evalue_text)
                bitscore = float(bitscore_text)
            except ValueError as exc:
                raise SearchOutputError(
                    "Non-numeric value in DIAMOND output at "
                    f"{output_file}:{line_number}."
                ) from exc

            if query_length <= 0:
                raise SearchOutputError(
                    "Query length must be greater than zero at "
                    f"{output_file}:{line_number}."
                )

            if target_length <= 0:
                raise SearchOutputError(
                    "Target length must be greater than zero at "
                    f"{output_file}:{line_number}."
                )

            query_coverage = (
                alignment_length / query_length * 100.0
            )
            target_coverage = (
                alignment_length / target_length * 100.0
            )

            hits.append(
                SearchHit(
                    query_id=query_id,
                    target_id=target_id,
                    identity=identity,
                    alignment_length=alignment_length,
                    query_length=query_length,
                    target_length=target_length,
                    query_coverage=query_coverage,
                    target_coverage=target_coverage,
                    evalue=evalue,
                    bitscore=bitscore,
                    method="diamond",
                )
            )

    return hits
