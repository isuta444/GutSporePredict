#!/usr/bin/env python3

"""
GutSporePredict v4.0-alpha1
Download the current GTDB bacterial metadata and phylogenetic tree.

Files:
    bac120_metadata.tsv.gz
    bac120.tree.gz
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
import urllib.error
import urllib.request
from datetime import datetime
from pathlib import Path


DEFAULT_BASE_URL = "https://data.gtdb.ecogenomic.org/releases/latest"

FILES = {
    "metadata": "bac120_metadata.tsv.gz",
    "tree": "bac120.tree.gz",
}


def calculate_md5(path: Path, chunk_size: int = 1024 * 1024) -> str:
    """Calculate the MD5 checksum of a file."""
    digest = hashlib.md5()

    with path.open("rb") as handle:
        while True:
            chunk = handle.read(chunk_size)

            if not chunk:
                break

            digest.update(chunk)

    return digest.hexdigest()


def download_file(
    url: str,
    output_path: Path,
    force: bool = False,
) -> dict:
    """Download one file with a temporary partial file."""

    output_path.parent.mkdir(parents=True, exist_ok=True)

    if output_path.exists() and not force:
        print(f"[SKIP] File already exists: {output_path}")

        return {
            "url": url,
            "path": str(output_path),
            "status": "existing",
            "size_bytes": output_path.stat().st_size,
            "md5": calculate_md5(output_path),
        }

    partial_path = Path(str(output_path) + ".part")

    if partial_path.exists():
        partial_path.unlink()

    print(f"[DOWNLOAD] {url}")
    print(f"[OUTPUT]   {output_path}")

    request = urllib.request.Request(
        url,
        headers={
            "User-Agent": "GutSporePredict-v4.0-alpha1",
        },
    )

    try:
        with urllib.request.urlopen(request, timeout=60) as response:
            total_size_header = response.headers.get("Content-Length")
            total_size = (
                int(total_size_header)
                if total_size_header is not None
                else None
            )

            downloaded = 0

            with partial_path.open("wb") as output_handle:
                while True:
                    block = response.read(1024 * 1024)

                    if not block:
                        break

                    output_handle.write(block)
                    downloaded += len(block)

                    if total_size:
                        percentage = downloaded / total_size * 100

                        print(
                            f"\r"
                            f"  {downloaded / 1024 / 1024:,.1f} MB"
                            f" / "
                            f"{total_size / 1024 / 1024:,.1f} MB"
                            f"  ({percentage:.1f}%)",
                            end="",
                            flush=True,
                        )
                    else:
                        print(
                            f"\r"
                            f"  {downloaded / 1024 / 1024:,.1f} MB",
                            end="",
                            flush=True,
                        )

            print()

    except urllib.error.HTTPError as error:
        if partial_path.exists():
            partial_path.unlink()

        raise RuntimeError(
            f"HTTP error {error.code} while downloading:\n{url}"
        ) from error

    except urllib.error.URLError as error:
        if partial_path.exists():
            partial_path.unlink()

        raise RuntimeError(
            f"Network error while downloading:\n{url}\n{error}"
        ) from error

    partial_path.replace(output_path)

    file_size = output_path.stat().st_size
    md5 = calculate_md5(output_path)

    print(f"[COMPLETE] {output_path}")
    print(f"[SIZE]     {file_size / 1024 / 1024:,.1f} MB")
    print(f"[MD5]      {md5}")

    return {
        "url": url,
        "path": str(output_path),
        "status": "downloaded",
        "size_bytes": file_size,
        "md5": md5,
    }


def test_gzip(path: Path) -> None:
    """Use gzip -t to check whether the downloaded gzip file is valid."""

    gzip_executable = shutil.which("gzip")

    if gzip_executable is None:
        print("[WARNING] gzip command was not found; integrity test skipped.")
        return

    import subprocess

    result = subprocess.run(
        [gzip_executable, "-t", str(path)],
        capture_output=True,
        text=True,
    )

    if result.returncode != 0:
        raise RuntimeError(
            f"gzip integrity test failed for:\n{path}\n"
            f"{result.stderr}"
        )

    print(f"[GZIP OK]  {path}")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Download GTDB bacterial metadata and tree."
    )

    parser.add_argument(
        "--output-dir",
        default="database/gtdb",
        help="Output directory. Default: database/gtdb",
    )

    parser.add_argument(
        "--base-url",
        default=DEFAULT_BASE_URL,
        help="GTDB release directory URL.",
    )

    parser.add_argument(
        "--metadata-only",
        action="store_true",
        help="Download metadata without the tree.",
    )

    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing files.",
    )

    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    requested_files = ["metadata"]

    if not args.metadata_only:
        requested_files.append("tree")

    records = []

    print("=" * 70)
    print("GutSporePredict v4.0-alpha1")
    print("GTDB data download")
    print("=" * 70)
    print(f"Base URL:   {args.base_url}")
    print(f"Output dir: {output_dir}")
    print()

    try:
        for file_key in requested_files:
            filename = FILES[file_key]
            url = f"{args.base_url.rstrip('/')}/{filename}"
            output_path = output_dir / filename

            record = download_file(
                url=url,
                output_path=output_path,
                force=args.force,
            )

            test_gzip(output_path)

            record["file_type"] = file_key
            records.append(record)

    except Exception as error:
        print()
        print(f"[ERROR] {error}", file=sys.stderr)
        return 1

    manifest = {
        "pipeline": "GutSporePredict",
        "version": "4.0-alpha1",
        "downloaded_at": datetime.now().astimezone().isoformat(),
        "base_url": args.base_url,
        "files": records,
    }

    manifest_path = output_dir / "download_manifest.json"

    manifest_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    print()
    print("=" * 70)
    print("[SUCCESS] GTDB files are ready.")
    print(f"[MANIFEST] {manifest_path}")
    print("=" * 70)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
