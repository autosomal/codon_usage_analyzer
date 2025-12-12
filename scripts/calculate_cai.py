#!/usr/bin/env python3
"""
Calculate Codon Adaptation Index (CAI) for genes
Based on Sharp & Li (1987) Nucleic Acids Research
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

def load_rscu_values(rscu_table_path):
    """
    Load RSCU values to determine the most frequent codon for each amino acid
    """
    df = pd.read_csv(rscu_table_path, sep='\t')
    return df

def calculate_relative_adaptiveness(codon_table_path, rscu_table_path):
    """
    Calculate relative adaptiveness values for each codon
    This is done by finding the most frequently used codon for each amino acid
    and scaling all other codons relative to it
    """
    # Load codon usage and RSCU data
    codon_df = load_codon_usage(codon_table_path)
    rscu_df = load_rscu_values(rscu_table_path)
    
    # Group by amino acid and find the codon with highest RSCU (most frequent)
    aa_groups = rscu_df.groupby('AminoAcid')
    reference_codons = {}
    
    for aa, group in aa_groups:
        if len(group) <= 1:  # Skip amino acids with only one codon
            continue
        # Find codon with highest RSCU value for this amino acid
        max_rscu_idx = group['RSCU'].idxmax()
        reference_codons[aa] = group.loc[max_rscu_idx, 'Codon']
    
    # Calculate relative adaptiveness (w) for each codon
    codon_adaptiveness = {}
    aa_groups_codon = codon_df.groupby('AminoAcid')
    
    for aa, group in aa_groups_codon:
        if aa not in reference_codons:
            # For amino acids with only one codon, set adaptiveness to 1.0
            if len(group) == 1:
                codon_adaptiveness[group.iloc[0]['Codon']] = 1.0
            continue
            
        # Get reference codon for this amino acid
        ref_codon = reference_codons[aa]
        
        # Calculate adaptiveness for all codons of this amino acid
        aa_codons = group.set_index('Codon')
        
        for codon in aa_codons.index:
            if codon == ref_codon:
                # Reference codon has adaptiveness of 1.0
                codon_adaptiveness[codon] = 1.0
            else:
                # For other codons, we need to calculate relative frequency
                # In a more complete implementation, this would be based on actual usage
                # For now, we'll calculate it based on the codon's relative frequency in the dataset
                ref_freq = aa_codons.loc[ref_codon, 'Count'] if ref_codon in aa_codons.index else 1
                codon_freq = aa_codons.loc[codon, 'Count'] if codon in aa_codons.index else 0
                if ref_freq > 0:
                    codon_adaptiveness[codon] = codon_freq / ref_freq if codon_freq <= ref_freq else 1.0
                else:
                    codon_adaptiveness[codon] = 0.0
    
    return codon_adaptiveness

def calculate_cai_for_gene(gene_data, codon_adaptiveness):
    """
    Calculate CAI for a single gene using the geometric mean method
    CAI = (product of wi for all codons)^(1/n)
    where wi is the relative adaptiveness of the i-th codon and n is the number of codons
    """
    # Filter to keep only codons with defined adaptiveness values
    gene_codons = gene_data[gene_data['Codon'].isin(codon_adaptiveness.keys())]
    
    if len(gene_codons) == 0:
        return 0.0
    
    # Get adaptiveness values for codons in this gene
    w_values = []
    for _, row in gene_codons.iterrows():
        codon = row['Codon']
        count = row['Count']
        if codon in codon_adaptiveness and codon_adaptiveness[codon] > 0:
            # Add the adaptiveness value 'count' times (for multiple occurrences)
            w_values.extend([codon_adaptiveness[codon]] * int(count))
    
    if len(w_values) == 0:
        return 0.0
    
    # Calculate geometric mean of w values
    # To avoid numerical issues with very small numbers, use log transformation
    log_w_sum = sum(np.log(w) for w in w_values if w > 0)
    n = len(w_values)
    
    if n == 0:
        return 0.0
    
    # CAI = exp(log(sum of w_i)/n) = geometric mean
    cai_value = np.exp(log_w_sum / n)
    
    return cai_value

def calculate_cai_for_all_genes(codon_table_path, rscu_table_path):
    """
    Calculate CAI for all genes in the codon usage table
    """
    # Calculate codon adaptiveness values
    codon_adaptiveness = calculate_relative_adaptiveness(codon_table_path, rscu_table_path)
    
    # Load codon usage data
    codon_df = load_codon_usage(codon_table_path)
    
    # Group by gene
    gene_groups = codon_df.groupby('Gene')
    
    # Calculate CAI for each gene
    cai_results = []
    for gene_id, gene_data in gene_groups:
        cai_value = calculate_cai_for_gene(gene_data, codon_adaptiveness)
        cai_results.append({'Gene': gene_id, 'CAI': cai_value})
    
    return pd.DataFrame(cai_results)

def main():
    # Get input and output file paths from Snakemake
    codon_table_path = snakemake.input.codon_table
    rscu_table_path = snakemake.input.rscu_table
    output_path = snakemake.output.cai_table
    
    # Calculate CAI values
    cai_df = calculate_cai_for_all_genes(codon_table_path, rscu_table_path)
    
    # Save results
    cai_df.to_csv(output_path, sep='\t', index=False)
    
    print(f"CAI values calculated and saved to {output_path}")

if __name__ == "__main__":
    # For Snakemake integration:
    main()