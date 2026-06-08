import sys
from pathlib import Path

# Standard translation codon table
CODON_TABLE = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
}

COMPLEMENT_MAP = str.maketrans('ATCGatcgNn', 'TAGCtagcNn')

def reverse_complement(seq):
    """Returns the reverse complement of a nucleotide sequence."""
    return seq.translate(COMPLEMENT_MAP)[::-1]

def translate(seq):
    """Translates a nucleotide sequence to protein using standard genetic code."""
    seq_upper = seq.upper()
    protein = []
    for i in range(0, len(seq_upper) - 2, 3):
        codon = seq_upper[i:i+3]
        protein.append(CODON_TABLE.get(codon, "X"))
    return "".join(protein)

def parse_fasta(filepath):
    """Parses a FASTA file and returns a dictionary of records."""
    records = {}
    current_name = None
    current_seq = []
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_name:
                    records[current_name] = "".join(current_seq)
                current_name = line[1:].split()[0]  # Get only the sequence ID
                current_seq = []
            else:
                current_seq.append(line)
        if current_name:
            records[current_name] = "".join(current_seq)
    return records

def parse_gff(filepath):
    """Parses a GFF3 file and returns CDS features (1-based, inclusive coordinates)."""
    cds_features = []
    with open(filepath, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue
            
            # Parse attributes (e.g. ID=cds_alpha;Name=alpha_orf)
            attrs = {}
            for item in parts[8].split(";"):
                if "=" in item:
                    key, val = item.split("=", 1)
                    attrs[key.strip()] = val.strip()
                    
            name = attrs.get("Name") or attrs.get("ID") or "unknown"
            
            cds_features.append({
                "chrom": parts[0],
                "start": int(parts[3]),
                "end": int(parts[4]),
                "strand": parts[6],
                "name": name
            })
    return cds_features

def main():
    genome_file = "genome.fa"
    gff_file = "annotations.gff3"
    
    if not Path(genome_file).exists() or not Path(gff_file).exists():
        print(f"Error: Make sure {genome_file} and {gff_file} are in the current directory.")
        sys.exit(1)
        
    # Load genome FASTA
    genome = parse_fasta(genome_file)
    
    # Load GFF annotations
    cds_features = parse_gff(gff_file)
    
    for cds in cds_features:
        chrom = cds["chrom"]
        if chrom not in genome:
            print(f"Warning: Chromosome {chrom} not found in genome.")
            continue
            
        seq = genome[chrom]
        
        # GFF coordinates are 1-based, inclusive.
        # Python slicing is 0-based, exclusive end.
        start_idx = cds["start"] - 1
        end_idx = cds["end"]
        
        nt_seq = seq[start_idx:end_idx]
        
        # Reverse complement if on the negative strand
        if cds["strand"] == "-":
            nt_seq = reverse_complement(nt_seq)
            
        # Translate to protein using the standard genetic code
        protein_seq = translate(nt_seq)
        
        print(f"{cds['name']}\t{nt_seq}\t{protein_seq}")

if __name__ == "__main__":
    main()
