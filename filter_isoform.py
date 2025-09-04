#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""
filter_isoforms.py

Filter FASTA entries down to one isoform per gene, supporting both Ensembl and
NCBI RefSeq-style headers.

Strategies:
  - longest (default): keep the isoform with the longest sequence
  - first: keep the first isoform encountered
  - curated: prefer curated accessions (RefSeq NP_ > XP_ > others; Ensembl ENSP > others),
             tie-break by length, then by lexicographic accession.

Gene detection (best-effort, common cases):
  Ensembl peptides/transcripts may include tokens like 'gene:', 'gene_symbol:'.
  RefSeq often includes 'gene=' or 'GN='; otherwise falls back to the first
  pipe/space-delimited token if it looks like a gene symbol.

Usage:
  python filter_isoforms.py -i input.faa -o filtered.faa
  python filter_isoforms.py -i input.faa.gz -o filtered.faa --strategy curated
"""

from __future__ import annotations
import argparse
import gzip
import io
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, Tuple, Optional

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# --- Heuristics ---------------------------------------------------------------

# Common patterns for gene extraction
RE_GENE_TOKENS = [
    re.compile(r"(?:^|\s)gene:([A-Za-z0-9_.-]+)"),
    re.compile(r"(?:^|\s)gene_symbol:([A-Za-z0-9_.-]+)"),
    re.compile(r"(?:^|\s)gene=([A-Za-z0-9_.-]+)"),
    re.compile(r"(?:^|\s)GN=([A-Za-z0-9_.-]+)"),
]

# Accession priority (higher is better) for 'curated' strategy
def accession_priority(acc: str) -> int:
    # RefSeq proteins: NP_ (curated) > XP_ (model) > YP_ (prok/vir) > others
    # Ensembl protein IDs: ENSP... generally curated reference set
    if acc.startswith("NP_"):
        return 100
    if acc.startswith("ENSP"):
        return 90
    if acc.startswith("XP_"):
        return 80
    if acc.startswith("YP_"):
        return 70
    if acc.startswith(("NM_", "XM_")):
        return 60
    if acc.startswith(("ENST", "ENSPTRP", "ENSPPAP", "ENSGGOP")):
        return 50
    return 10

def open_maybe_gzip(path: Path):
    # Transparent open for .gz or plain files
    if str(path).endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8")
    return open(path, "r", encoding="utf-8")

def extract_gene(record: SeqRecord) -> Optional[str]:
    """
    Try to extract a gene identifier (symbol or stable gene ID) from header.
    Returns None if no reasonable gene could be inferred.
    """
    desc = record.description or record.id or ""
    # 1) Look for explicit gene tokens
    for rx in RE_GENE_TOKENS:
        m = rx.search(desc)
        if m:
            return m.group(1)

    # 2) Ensembl FASTA headers often contain 'gene:' or gene_symbol in long desc
    # already covered by RE_GENE_TOKENS.

    # 3) Fallbacks:
    # - If ID looks like SYMBOL|ACCESSION → take SYMBOL
    if "|" in record.id:
        left = record.id.split("|", 1)[0]
        if re.fullmatch(r"[A-Za-z0-9_.-]+", left):
            return left

    # - If first token (before space) looks like a plausible gene symbol
    first = desc.split()[0]
    if re.fullmatch(r"[A-Za-z0-9_.-]+", first) and not first.startswith(("ENSP", "ENST", "ENSG", "XP_", "NP_", "XM_", "NM_")):
        return first

    # As a last resort, if it's clearly an Ensembl gene ID in description (ENSG...)
    m = re.search(r"(ENSG[0-9]+)", desc)
    if m:
        return m.group(1)

    return None

def get_accession(record: SeqRecord) -> str:
    """
    Return an accession-like token from the record ID/description to help ranking.
    """
    # Prefer the raw ID if it looks accession-like
    rid = record.id
    if re.match(r"^(?:[A-Z]{2,4}_\d+(\.\d+)?|ENS[TPG][A-Z]*\d+(\.\d+)?)$", rid):
        return rid
    # Otherwise search for a familiar accession inside description
    m = re.search(r"(NP|XP|YP|NM|XM)_(\d+)(\.\d+)?", record.description)
    if m:
        return m.group(0)
    m = re.search(r"(ENSP|ENST|ENSG)[A-Z]*\d+(\.\d+)?", record.description)
    if m:
        return m.group(0)
    return rid

# --- Filtering strategies -----------------------------------------------------

def choose_isoform(records: list[SeqRecord], strategy: str) -> SeqRecord:
    if strategy == "first":
        return records[0]
    if strategy == "longest":
        return max(records, key=lambda r: (len(r.seq), -records.index(r)))
    if strategy == "curated":
        return max(
            records,
            key=lambda r: (accession_priority(get_accession(r)), len(r.seq), -records.index(r)),
        )
    raise ValueError(f"Unknown strategy: {strategy}")

# --- CLI ----------------------------------------------------------------------

def parse_args():
    ap = argparse.ArgumentParser(
        description="Filter isoforms to one sequence per gene (Ensembl & RefSeq headers).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("-i", "--input", required=True, type=Path, help="Input FASTA (.fa/.faa/.fasta[.gz])")
    ap.add_argument("-o", "--output", required=True, type=Path, help="Output FASTA with one isoform per gene")
    ap.add_argument("--strategy", choices=["longest", "first", "curated"], default="longest",
                    help="Isoform selection strategy")
    ap.add_argument("--write-map", type=Path, default=None,
                    help="Optional TSV: gene_id   kept_record_id")
    return ap.parse_args()

def main():
    args = parse_args()

    # Group records by gene
    gene_to_records: Dict[str, list[SeqRecord]] = defaultdict(list)
    total = 0
    with open_maybe_gzip(args.input) as handle:
        for rec in SeqIO.parse(handle, "fasta"):
            total += 1
            gene = extract_gene(rec)
            if gene:
                gene_to_records[gene].append(rec)

    if not gene_to_records:
        raise SystemExit("No genes could be extracted from headers. Check your FASTA format.")

    # Choose one isoform per gene
    kept: list[SeqRecord] = []
    kept_map: list[Tuple[str, str]] = []
    for gene, recs in gene_to_records.items():
        chosen = choose_isoform(recs, args.strategy)
        # Clean header: keep record.id, store gene as prefix to be explicit
        acc = get_accession(chosen)
        chosen.id = f"{gene}|{acc}"
        chosen.description = ""
        kept.append(chosen)
        kept_map.append((gene, chosen.id))

    # Write outputs
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", encoding="utf-8") as out_f:
        SeqIO.write(kept, out_f, "fasta")

    if args.write_map:
        args.write_map.parent.mkdir(parents=True, exist_ok=True)
        with open(args.write_map, "w", encoding="utf-8") as m:
            for gene, chosen_id in kept_map:
                m.write(f"{gene}\t{chosen_id}\n")

    print(f"[ok] Read {total} records; kept {len(kept)} genes → {args.output}")
    if args.write_map:
        print(f"[ok] Wrote map → {args.write_map}")

if __name__ == "__main__":
    main()
