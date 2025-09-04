#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""
isoform_check.py

Count how many isoforms per gene are present in a FASTA file.

Usage:
    python isoform_check.py -i transcripts.fasta
"""

from __future__ import annotations
from collections import defaultdict
from pathlib import Path
import argparse
import re
from Bio import SeqIO

def parse_args():
    p = argparse.ArgumentParser(description="Check isoform counts per gene from a FASTA file.")
    p.add_argument("-i", "--input", required=True, type=Path, help="Input FASTA file")
    p.add_argument("-n", "--top", type=int, default=10, help="Show top N genes with most isoforms")
    return p.parse_args()

def main():
    args = parse_args()
    isoform_counts = defaultdict(int)

    for record in SeqIO.parse(str(args.input), "fasta"):
        # Match gene ID at start of header, e.g. LOC12345 or a word
        m = re.match(r"(LOC\d+|\w+)", record.description)
        if m:
            gene_id = m.group(1)
            isoform_counts[gene_id] += 1

    total_genes = len(isoform_counts)
    total_isoforms = sum(isoform_counts.values())
    multi_iso = {g:c for g,c in isoform_counts.items() if c > 1}

    print(f"Total genes: {total_genes}")
    print(f"Total isoforms: {total_isoforms}")
    print(f"Genes with multiple isoforms: {len(multi_iso)}")

    print(f"\nTop {args.top} genes with most isoforms:")
    for gene, count in sorted(isoform_counts.items(), key=lambda x: x[1], reverse=True)[:args.top]:
        print(f"{gene}: {count}")

if __name__ == "__main__":
    main()
