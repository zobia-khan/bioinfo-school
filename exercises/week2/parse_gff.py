#!/usr/bin/env python3
import sys
from pathlib import Path
from collections import defaultdict

def parse_gff(filepath):
    """
    Parses a GFF3 file and groups CDS features by their parent gene/transcript.
    
    Returns:
        A dictionary mapping gene/transcript IDs to a list of CDS features (start, end, length).
    """
    genes_cds = defaultdict(list)
    
    with open(filepath, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
                
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
                
            feature_type = parts[2]
            if feature_type != "CDS":
                continue
                
            # Parse coordinates (1-based, inclusive)
            try:
                start = int(parts[3])
                end = int(parts[4])
            except ValueError:
                continue
                
            length = end - start + 1
            
            # Parse attributes
            attrs = {}
            for item in parts[8].split(";"):
                if "=" in item:
                    key, val = item.split("=", 1)
                    attrs[key.strip()] = val.strip()
            
            # Find a suitable gene/transcript identifier:
            # 1. Parent attribute (standard for GFF3 CDS)
            # 2. gene_id or gene_name or Name
            # 3. ID (if it's a standalone CDS)
            gene_id = attrs.get("Parent") or attrs.get("gene_id") or attrs.get("gene_name") or attrs.get("Name") or attrs.get("ID") or "unknown_gene"
            
            genes_cds[gene_id].append({
                "start": start,
                "end": end,
                "length": length,
                "chrom": parts[0],
                "strand": parts[6]
            })
            
    return genes_cds

def main():
    if len(sys.argv) < 2:
        print("Usage: python parse_gff.py <path_to_gff3>")
        sys.exit(1)
        
    gff_path = sys.argv[1]
    if not Path(gff_path).exists():
        print(f"Error: File '{gff_path}' not found.")
        sys.exit(1)
        
    genes_cds = parse_gff(gff_path)
    
    print(f"{'Gene/Transcript':<25} {'Num CDS Exons':<15} {'Max Exon Length':<18} {'Total CDS Length':<18}")
    print("-" * 78)
    
    for gene_id, cds_list in sorted(genes_cds.items()):
        total_len = sum(cds["length"] for cds in cds_list)
        max_exon_len = max(cds["length"] for cds in cds_list)
        num_exons = len(cds_list)
        print(f"{gene_id:<25} {num_exons:<15} {max_exon_len:<18} {total_len:<18}")

if __name__ == "__main__":
    main()
