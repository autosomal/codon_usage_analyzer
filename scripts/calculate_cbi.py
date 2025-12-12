#!/usr/bin/env python3
"""
Calculate Codon Bias Index (CBI) for genes
Based on Bennetzen & Hall (1982) Nucleic Acids Research
"""

import pandas as pd
import numpy as np
import sys
import os
from collections import defaultdict

def load_codon_usage(codon_table_path):
    """
    Load codon usage table from file
    Expected columns: Gene, Codon, Count, AminoAcid
    """
    df = pd.read_csv(codon_table_path, sep='\t')
    return df

def calculate_cbi_for_gene(gene_data):
    """
    Calculate CBI for a single gene
    CBI = sum(|fi - fmax|) / (2 * (Naa - 1) / Naa)
    where fi is frequency of codon i, fmax is max frequency in synonymous family,
    Naa is number of amino acids in the gene
    """
    # Group by amino acid to get synonymous codon families
    aa_groups = gene_data.groupby('AminoAcid')
    
    # Calculate CBI components
    total_cbi = 0.0
    total_codons_in_syn_families = 0
    
    for aa, group in aa_groups:
        # Skip Met, Trp, and stop codons (only one codon per amino acid)
        if len(group) <= 1:
            continue
            
        # Calculate relative frequencies for this amino acid
        codon_counts = group['Count'].values
        total_aa_count = codon_counts.sum()
        
        if total_aa_count == 0:
            continue
            
        codon_freqs = codon_counts / total_aa_count
        max_freq = codon_freqs.max()
        
        # Calculate sum of absolute differences from max frequency
        abs_diff_sum = np.sum(np.abs(codon_freqs - max_freq))
        
        # Weight by number of codons in this amino acid family
        total_cbi += abs_diff_sum
        total_codons_in_syn_families += len(group)  # number of synonymous codons
    
    # Calculate normalization factor
    # For CBI, the formula is slightly different depending on implementation
    # Here we implement the original formula: CBI = sum(|fi - fmax|) / (2*(Nc-1)/Nc)
    # where Nc is the number of codons in the gene
    if total_codons_in_syn_families <= 1:
        return 0.0
    
    # Calculate CBI
    # Original CBI formula: CBI = sum(|fi - fmax|) / (2*(Nc-1)/Nc)
    # But typically for CBI: CBI = sum(fi_max - fi) / (N_synonymous_families)
    # where fi_max is the max frequency in each synonymous family
    
    # Let's recalculate using the proper CBI formula
    cbi_numerator = 0.0
    synonymous_family_count = 0
    
    for aa, group in aa_groups:
        if len(group) <= 1:  # Skip amino acids with only one codon
            continue
        
        codon_counts = group['Count'].values
        total_aa_count = codon_counts.sum()
        
        if total_aa_count == 0:
            continue
        
        codon_freqs = codon_counts / total_aa_count
        max_freq = codon_freqs.max()
        
        # Sum of (max_freq - each_freq) for this amino acid family
        cbi_contribution = len(group) * (max_freq - np.mean(codon_freqs))
        cbi_numerator += (max_freq - codon_freqs).sum()
        synonymous_family_count += 1
    
    if synonymous_family_count == 0:
        return 0.0
    
    # Proper CBI calculation
    cbi_value = cbi_numerator / synonymous_family_count
    
    return cbi_value

def calculate_cbi_for_all_genes(codon_table_path):
    """
    Calculate CBI for all genes in the codon usage table
    """
    # Load data
    codon_df = load_codon_usage(codon_table_path)
    
    # Group by gene
    gene_groups = codon_df.groupby('Gene')
    
    # Calculate CBI for each gene
    cbi_results = []
    for gene_id, gene_data in gene_groups:
        cbi_value = calculate_cbi_for_gene(gene_data)
        cbi_results.append({'Gene': gene_id, 'CBI': cbi_value})
    
    return pd.DataFrame(cbi_results)

def main():
    # Get input and output file paths from Snakemake
    codon_table_path = snakemake.input.codon_table
    output_path = snakemake.output.cbi_table
    
    # Calculate CBI values
    cbi_df = calculate_cbi_for_all_genes(codon_table_path)
    
    # Save results
    cbi_df.to_csv(output_path, sep='\t', index=False)
    
    print(f"CBI values calculated and saved to {output_path}")

if __name__ == "__main__":
    # For Snakemake integration:
    main()