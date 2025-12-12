#!/usr/bin/env python3
"""
Calculate codon usage from parsed CodonW table
"""

import pandas as pd
import numpy as np
import sys
import os

def load_parsed_codonw(parsed_file_path):
    """
    Load parsed CodonW output file
    Expected format: gene-specific codon usage data
    """
    df = pd.read_csv(parsed_file_path, sep='\t')
    return df

def calculate_codon_usage_from_table(parsed_file_path):
    """
    Calculate comprehensive codon usage statistics from parsed table
    """
    # Load the parsed data
    df = load_parsed_codonw(parsed_file_path)
    
    # Assuming the parsed file has columns: Gene, Codon, Count, AminoAcid
    # If the structure is different, we'll need to adjust accordingly
    
    # Create codon usage table
    codon_usage_data = []
    
    # Group by gene and codon to get counts
    gene_codon_groups = df.groupby(['Gene', 'Codon'])
    
    for (gene, codon), group in gene_codon_groups:
        # Calculate total for this gene
        gene_total = df[df['Gene'] == gene]['Count'].sum()
        
        # Get count for this specific codon
        codon_count = group['Count'].sum()
        
        # Calculate frequency
        frequency = codon_count / gene_total if gene_total > 0 else 0
        
        # Get amino acid (assuming it's in the data)
        amino_acid = group['AminoAcid'].iloc[0] if 'AminoAcid' in group.columns else 'Unknown'
        
        codon_usage_data.append({
            'Gene': gene,
            'Codon': codon,
            'AminoAcid': amino_acid,
            'Count': codon_count,
            'Frequency': frequency
        })
    
    codon_usage_df = pd.DataFrame(codon_usage_data)
    
    # Also calculate amino acid usage
    aa_usage_data = []
    gene_aa_groups = codon_usage_df.groupby(['Gene', 'AminoAcid'])
    
    for (gene, aa), group in gene_aa_groups:
        total_aa_count = group['Count'].sum()
        aa_usage_data.append({
            'Gene': gene,
            'AminoAcid': aa,
            'Count': total_aa_count
        })
    
    amino_acid_df = pd.DataFrame(aa_usage_data)
    
    return codon_usage_df, amino_acid_df

def main():
    # Get input and output file paths from Snakemake
    parsed_file_path = snakemake.input.parsed_file
    codon_table_path = snakemake.output.codon_table
    amino_acid_table_path = snakemake.output.amino_acid_table
    
    # Calculate codon usage
    codon_usage_df, amino_acid_df = calculate_codon_usage_from_table(parsed_file_path)
    
    # Save results
    codon_usage_df.to_csv(codon_table_path, sep='\t', index=False)
    amino_acid_df.to_csv(amino_acid_table_path, sep='\t', index=False)
    
    print(f"Codon usage table saved to {codon_table_path}")
    print(f"Amino acid usage table saved to {amino_acid_table_path}")

if __name__ == "__main__":
    # For Snakemake integration:
    main()