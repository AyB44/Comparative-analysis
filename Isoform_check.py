#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""
isoform_check.py

Count how many isoforms per gene are present in a FASTA file.

Usage:
    python isoform_check.py -i transcripts.fasta
"""
# Enable future compatibility for postponed evaluation of type annotations
from __future__ import annotations
# Import necessary modules for data structures, file handling, argument parsing, regex and FASTA parsing
from collections import defaultdict
from pathlib import Path
import argparse
import re
from Bio import SeqIO

# Define a function for command-line arguments for input file and number of top genes to display
def parse_args():
    p = argparse.ArgumentParser(description="Check isoform counts per gene from a FASTA file.")
    p.add_argument("-i", "--input", required=True, type=Path, help="Input FASTA file")
    p.add_argument("-n", "--top", type=int, default=10, help="Show top N genes with most isoforms")
    return p.parse_args()

# Main function to process the FASTA file and compute isoform statistics
def main():
    args = parse_args()
    isoform_counts = defaultdict(int)

    for record in SeqIO.parse(str(args.input), "fasta"):
        # Match gene ID at start of header and change according to it
        m = re.match(r"(LOC\d+|\w+)", record.description)
        # For Ensembl headers use:
        # m = re.search(r"gene:(ENSG\d+)", record.description)
        if m:
            gene_id = m.group(1)
            isoform_counts[gene_id] += 1
            
# Calculate total unique genes, total isoforms and gene with more than one isoforms
    total_genes = len(isoform_counts)
    total_isoforms = sum(isoform_counts.values())
    multi_iso = {g:c for g,c in isoform_counts.items() if c > 1}
    
# Print summary statistics
    print(f"Total genes: {total_genes}")
    print(f"Total isoforms: {total_isoforms}")
    print(f"Genes with multiple isoforms: {len(multi_iso)}")
# Print the top N genes with the highest number of isoforms
    print(f"\nTop {args.top} genes with most isoforms:")
    for gene, count in sorted(isoform_counts.items(), key=lambda x: x[1], reverse=True)[:args.top]:
        print(f"{gene}: {count}")

if __name__ == "__main__":
    main()
