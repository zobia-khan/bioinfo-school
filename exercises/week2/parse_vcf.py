#!/usr/bin/env python3
import sys
from pathlib import Path
from collections import defaultdict

def parse_info_field(info_str):
    """Parses the INFO string in a VCF line into a key-value dictionary."""
    info_dict = {}
    if info_str == ".":
        return info_dict
    for item in info_str.split(";"):
        if "=" in item:
            key, val = item.split("=", 1)
            info_dict[key.strip()] = val.strip()
        else:
            info_dict[item.strip()] = True
    return info_dict

def calculate_gt_frequency(gt_format, format_parts, sample_parts):
    """
    Calculates the alternative allele frequency from the sample genotype fields.
    Returns (ac, an, frequency) or (None, None, None) if GT is missing.
    """
    if not gt_format or "GT" not in gt_format.split(":"):
        return None, None, None
        
    gt_index = gt_format.split(":").index("GT")
    
    total_alleles = 0
    alt_alleles = 0
    
    for sample in sample_parts:
        if sample == "." or ":" not in sample and sample == "./.":
            continue
            
        sample_subparts = sample.split(":")
        if len(sample_subparts) <= gt_index:
            continue
            
        gt = sample_subparts[gt_index]
        # Genotype can be separated by '/' (unphased) or '|' (phased)
        delimiters = ["/", "|"]
        sep = None
        for d in delimiters:
            if d in gt:
                sep = d
                break
                
        if sep:
            alleles = gt.split(sep)
        else:
            alleles = [gt]  # Haploid organism/sex chromosome
            
        for allele in alleles:
            if allele == ".":
                continue  # Missing genotype
            total_alleles += 1
            # 0 is REF, any integer > 0 is ALT
            try:
                allele_idx = int(allele)
                if allele_idx > 0:
                    alt_alleles += 1
            except ValueError:
                continue
                
    if total_alleles == 0:
        return 0, 0, 0.0
    return alt_alleles, total_alleles, alt_alleles / total_alleles

def parse_vcf(filepath):
    """
    Parses a VCF file and calculates:
    1. Allele frequencies based on INFO/AF tag.
    2. Allele frequencies calculated from Genotypes (GT).
    """
    chrom_info_af = defaultdict(list)
    chrom_gt_af = defaultdict(list)
    
    with open(filepath, "r") as f:
        header_cols = []
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                header_cols = line.strip().split("\t")
                continue
                
            parts = line.strip().split("\t")
            if len(parts) < 8:
                continue
                
            chrom = parts[0]
            info_str = parts[7]
            info_dict = parse_info_field(info_str)
            
            # 1. Parse AF from INFO
            if "AF" in info_dict:
                af_val = info_dict["AF"]
                # AF can be comma-separated for multi-allelic sites
                try:
                    for val in af_val.split(","):
                        chrom_info_af[chrom].append(float(val))
                except ValueError:
                    pass
            
            # 2. Parse from Genotypes (if sample columns are present)
            if len(parts) > 9 and len(header_cols) > 9:
                format_field = parts[8]
                samples = parts[9:]
                ac, an, gt_af = calculate_gt_frequency(format_field, format_field, samples)
                if gt_af is not None:
                    chrom_gt_af[chrom].append(gt_af)
                    
    return chrom_info_af, chrom_gt_af

def main():
    if len(sys.argv) < 2:
        print("Usage: python parse_vcf.py <path_to_vcf>")
        sys.exit(1)
        
    vcf_path = sys.argv[1]
    if not Path(vcf_path).exists():
        print(f"Error: File '{vcf_path}' not found.")
        sys.exit(1)
        
    info_af, gt_af = parse_vcf(vcf_path)
    
    all_chroms = sorted(list(set(info_af.keys()) | set(gt_af.keys())))
    
    print(f"{'Chromosome':<15} | {'INFO/AF (Mean)':<20} | {'Genotype/AF (Mean)':<20}")
    print("-" * 65)
    
    for chrom in all_chroms:
        # Calculate mean INFO/AF
        info_list = info_af.get(chrom, [])
        info_mean_str = f"{sum(info_list)/len(info_list):.4f} (n={len(info_list)})" if info_list else "N/A"
        
        # Calculate mean Genotype/AF
        gt_list = gt_af.get(chrom, [])
        gt_mean_str = f"{sum(gt_list)/len(gt_list):.4f} (n={len(gt_list)})" if gt_list else "N/A"
        
        print(f"{chrom:<15} | {info_mean_str:<20} | {gt_mean_str:<20}")

if __name__ == "__main__":
    main()
