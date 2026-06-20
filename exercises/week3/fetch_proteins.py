#!/usr/bin/env python3
"""Fetch Week 3 starter proteins from UniProt and write proteins.fasta."""

from __future__ import annotations

import csv
import sys
import textwrap
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import urlopen


ROOT = Path(__file__).resolve().parent
ACCESSIONS = ROOT / "protein_accessions.tsv"
OUTPUT_FASTA = ROOT / "proteins.fasta"
UNIPROT_FASTA_URL = "https://rest.uniprot.org/uniprotkb/{accession}.fasta"


def read_accessions(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(
            (line for line in handle if not line.startswith("#") and line.strip()),
            delimiter="\t",
            fieldnames=["family", "accession", "label"],
        )
        return list(reader)


def fetch_sequence(accession: str) -> tuple[str, str]:
    url = UNIPROT_FASTA_URL.format(accession=accession)
    try:
        with urlopen(url, timeout=30) as response:
            fasta = response.read().decode("utf-8")
    except (HTTPError, URLError, TimeoutError) as exc:
        raise RuntimeError(f"Could not fetch {accession} from UniProt: {exc}") from exc

    lines = [line.strip() for line in fasta.splitlines() if line.strip()]
    if not lines or not lines[0].startswith(">"):
        raise RuntimeError(f"UniProt returned an unexpected FASTA response for {accession}")

    description = lines[0][1:]
    sequence = "".join(lines[1:])
    return description, sequence


def write_fasta(records: list[dict[str, str]], output_path: Path) -> None:
    with output_path.open("w") as handle:
        for record in records:
            description, sequence = fetch_sequence(record["accession"])
            header = (
                f">{record['label']}|family={record['family']}|"
                f"accession={record['accession']} {description}"
            )
            handle.write(header + "\n")
            handle.write("\n".join(textwrap.wrap(sequence, width=80)) + "\n")


def main() -> int:
    records = read_accessions(ACCESSIONS)
    if not records:
        print(f"No accessions found in {ACCESSIONS}", file=sys.stderr)
        return 1

    write_fasta(records, OUTPUT_FASTA)
    print(f"Wrote {len(records)} sequences to {OUTPUT_FASTA}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
