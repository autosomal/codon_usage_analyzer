#!/usr/bin/env python3
"""
Calculate Frequency of Optimal Codons (FOP) for genes
Based on Ikemura (1981) Journal of Molecular Biology
"""

import pandas as pd
import numpy as np
import sys
import os

def load_codon_usage(codon_table_path):
    """
    Load codon usage table from file
    Expected columns: Gene, Codon, Count, AminoAcid
    """
    df = pd.read_csv(codon_table_path, sep='\t')
    return df

def load_preferred_codons(preferred_codons_path):
    """
    Load list of preferred codons from file
    Expected format: AminoAcid, PreferredCodon
    """
    df = pd.read_csv(preferred_codons_path, sep='\t')
    # Create a dictionary mapping amino acid to its preferred codon(s)
    preferred_dict = {}
    for _, row in df.iterrows():
        aa = row['AminoAcid']
        codon = row['Codon']
        if aa not in preferred_dict:
            preferred_dict[aa] = []
        preferred_dict[aa].append(codon)
    return preferred_dict

def calculate_fop_for_gene(gene_data, preferred_codons_dict):
    """
    Calculate FOP for a single gene
    FOP = (No. of optimal codons in gene) / (Total number of codons in gene)
    """
    # Count total codons in the gene (excluding Met, Trp, and Stop which have only one codon)
    total_codons = 0
    optimal_codons = 0
    
    # Group by amino acid to handle synonymous codon families
    aa_groups = gene_data.groupby('AminoAcid')
    
    for aa, group in aa_groups:
        # Skip amino acids with only one codon (Met, Trp, and Stop)
        if len(group) <= 1:
            continue
        
        # Count total codons for this amino acid
        aa_total_codons = group['Count'].sum()
        total_codons += aa_total_codons
        
        # Check if this amino acid has preferred codons defined
        if aa in preferred_codons_dict:
            preferred_codon_list = preferred_codons_dict[aa]
            
            # Count how many of the codons in this gene are preferred
            for _, codon_row in group.iterrows():
                if codon_row['Codon'] in preferred_codon_list:
                    optimal_codons += codon_row['Count']
    
    # Calculate FOP
    if total_codons == 0:
        return 0.0
    
    fop_value = optimal_codons / total_codons
    return fop_value

def calculate_fop_for_all_genes(codon_table_path, preferred_codons_path):
    """
    Calculate FOP for all genes in the codon usage table
    """
    # Load data
    codon_df = load_codon_usage(codon_table_path)
    preferred_codons_dict = load_preferred_codons(preferred_codons_path)
    
    # Group by gene
    gene_groups = codon_df.groupby('Gene')
    
    # Calculate FOP for each gene
    fop_results = []
    for gene_id, gene_data in gene_groups:
        fop_value = calculate_fop_for_gene(gene_data, preferred_codons_dict)
        fop_results.append({'Gene': gene_id, 'FOP': fop_value})
    
    return pd.DataFrame(fop_results)

def main():
    # Get input and output file paths from Snakemake
    codon_table_path = snakemake.input.codon_table
    preferred_codons_path = snakemake.input.preferred_codons
    output_path = snakemake.output.fop_table
    
    # Calculate FOP values
    fop_df = calculate_fop_for_all_genes(codon_table_path, preferred_codons_path)
    
    # Save results
    fop_df.to_csv(output_path, sep='\t', index=False)
    
    print(f"FOP values calculated and saved to {output_path}")

if __name__ == "__main__":
    # For Snakemake integration:
    main()