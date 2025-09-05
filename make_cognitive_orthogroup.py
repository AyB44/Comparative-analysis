# Filter orthogroups based on cognitive genes

# --- cognitive_orthogroups_filter.py ---
import re
cognitive_gene_file = "cognitive_genes.txt" # File with the list of all cognitive genes extracted using eggNOG-mapper
orthogroups_file = "Orthogroups.tsv"  # Output file from Orthofinder
outfile = "cognitive.tsv"  # Output containing congitive orthogroups

# Load cognitive genes into a set
with open(cognitive_gene_file, "r") as f:
    cognitive_genes = set()
    for line in f:
        line = line.strip()
        if not line:
            continue
        # Extract gene/protein identifier from various formats
        # extract ENSG, ENST, ENSP, XP_, or full fallback as needed
        match = re.search(r"(ENSG\w+|ENST\w+|ENSP\w+|XP_\w+)", line)
        if match:
            cognitive_genes.add(match.group(1))
        else:
            cognitive_genes.add(line)  # fallback if no pattern matched

print(f"Loaded {len(cognitive_genes)} cognitive gene identifiers.")

def extract_all_ids(line):
    """
    Extract gene/protein identifiers from an orthogroup line.
    Handles mixed formatting (commas, pipes, version numbers).
    """
    ids = set()
    parts = line.strip().split('\t')[1:]  # skip Orthogroup ID
    for field in parts:
        for item in field.split(','):
            item = item.strip()
            if not item:
                continue
            # Extract common identifiers
            matches = re.findall(r"(ENSG\w+|ENST\w+|ENSP\w+|XP_\w+)", item)
            ids.update(matches)
    return ids

with open(orthogroups_file, "r") as infile, open(outfile, "w") as out:
    header = infile.readline()
    out.write(header)  # Write header
    for line in infile:
        gene_ids = extract_all_ids(line)
        if cognitive_genes & gene_ids:
            out.write(line)

print(f" cognitive.tsv written: filtered orthogroups with cognitive genes.")
