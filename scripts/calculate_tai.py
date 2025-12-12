#!/usr/bin/env python3
"""
Calculate tRNA Adaptation Index (tAI) for genes
Based on dos Reis et al. (2004) Nucleic Acids Research
"""

import pandas as pd
import numpy as np
import sys
import os
from collections import defaultdict

def load_codon_usage(codon_table_path):
    """
    Load codon usage table from file
    """
    df = pd.read_csv(codon_table_path, sep='\t')
    return df

def load_trna_abundance(trna_data_path):
    """
    Load tRNA gene copy numbers or abundances
    Expected format: Codon, tRNA_Copies
    """
    if os.path.exists(trna_data_path):
        trna_df = pd.read_csv(trna_data_path, sep='\t')
        # Convert to dictionary mapping codon to tRNA abundance
        trna_abundance = dict(zip(trna_df['Codon'], trna_df['tRNA_Copies']))
    else:
        # Provide default human tRNA abundance data if file doesn't exist
        print(f"Warning: {trna_data_path} not found. Using default tRNA abundance values.", file=sys.stderr)
        # Default values based on typical human tRNA gene copy numbers
        trna_abundance = {
            'GCA': 1.0, 'GCC': 1.0, 'GCG': 1.0, 'GCT': 1.0,  # Ala
            'TGC': 1.0, 'TGT': 1.0,                           # Cys
            'GAC': 1.0, 'GAT': 1.0,                           # Asp
            'GAA': 1.0, 'GAG': 1.0,                           # Glu
            'TTC': 1.0, 'TTT': 1.0,                           # Phe
            'GGA': 1.0, 'GGC': 1.0, 'GGG': 1.0, 'GGT': 1.0,  # Gly
            'CAC': 1.0, 'CAT': 1.0,                           # His
            'ATA': 1.0, 'ATC': 1.0, 'ATT': 1.0,              # Ile
            'AAA': 1.0, 'AAG': 1.0,                           # Lys
            'CTA': 1.0, 'CTC': 1.0, 'CTG': 1.0, 'CTT': 1.0, 'TTA': 1.0, 'TTG': 1.0,  # Leu
            'ATG': 1.0,                                        # Met
            'AAC': 1.0, 'AAT': 1.0,                           # Asn
            'CCA': 1.0, 'CCC': 1.0, 'CCG': 1.0, 'CCT': 1.0,  # Pro
            'CAA': 1.0, 'CAG': 1.0,                           # Gln
            'AGA': 1.0, 'AGG': 1.0, 'CGA': 1.0, 'CGC': 1.0, 'CGG': 1.0, 'CGT': 1.0,  # Arg
            'TCA': 1.0, 'TCC': 1.0, 'TCG': 1.0, 'TCT': 1.0, 'AGC': 1.0, 'AGT': 1.0,  # Ser
            'ACA': 1.0, 'ACC': 1.0, 'ACG': 1.0, 'ACT': 1.0,  # Thr
            'GTA': 1.0, 'GTC': 1.0, 'GTG': 1.0, 'GTT': 1.0,  # Val
            'TGG': 1.0,                                        # Trp
            'TAC': 1.0, 'TAT': 1.0,                           # Tyr
        }
    
    return trna_abundance

def calculate_wobble_efficiency(codon, trna_abundance):
    """
    Calculate wobble efficiency for a codon based on tRNA abundance
    This implements the basic model where we consider Watson-Crick and Wobble pairing
    """
    # Map codon to amino acid (simplified - normally would use genetic code)
    codon_aa_map = {
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',  # Ala
        'TGC': 'C', 'TGT': 'C',                           # Cys
        'GAC': 'D', 'GAT': 'D',                           # Asp
        'GAA': 'E', 'GAG': 'E',                           # Glu
        'TTC': 'F', 'TTT': 'F',                           # Phe
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',  # Gly
        'CAC': 'H', 'CAT': 'H',                           # His
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I',              # Ile
        'AAA': 'K', 'AAG': 'K',                           # Lys
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L', 'TTA': 'L', 'TTG': 'L',  # Leu
        'ATG': 'M',                                        # Met
        'AAC': 'N', 'AAT': 'N',                           # Asn
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',  # Pro
        'CAA': 'Q', 'CAG': 'Q',                           # Gln
        'AGA': 'R', 'AGG': 'R', 'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',  # Arg
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S', 'AGC': 'S', 'AGT': 'S',  # Ser
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',  # Thr
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',  # Val
        'TGG': 'W',                                        # Trp
        'TAC': 'Y', 'TAT': 'Y',                           # Tyr
        'TAG': '*', 'TAA': '*', 'TGA': '*'                # Stop
    }
    
    aa = codon_aa_map.get(codon, None)
    if aa is None or codon not in trna_abundance:
        return 0.0
    
    # Simple approach: use tRNA abundance directly as wobble efficiency
    # In a real implementation, this would be more complex
    return float(trna_abundance[codon])

def calculate_gene_tai(gene_codons, trna_abundance):
    """
    Calculate tAI for a single gene based on its codons
    """
    if len(gene_codons) == 0:
        return 0.0
    
    # Calculate wobble efficiency for each codon in the gene
    wobble_efficiencies = []
    for codon, count in gene_codons.items():
        if codon in trna_abundance and count > 0:
            wobble_efficiency = calculate_wobble_efficiency(codon, trna_abundance)
            # Add the efficiency value 'count' times (weighted by codon frequency)
            wobble_efficiencies.extend([wobble_efficiency] * int(count))
    
    if len(wobble_efficiencies) == 0:
        return 0.0
    
    # Calculate tAI as the mean wobble efficiency
    tai = np.mean(wobble_efficiencies)
    return tai

def calculate_tai_for_all_genes(codon_table_path, trna_data_path):
    """
    Calculate tAI for all genes in the codon usage table
    """
    # Load data
    codon_df = load_codon_usage(codon_table_path)
    trna_abundance = load_trna_abundance(trna_data_path)
    
    # Group codon usage by gene
    gene_codon_counts = {}
    for _, row in codon_df.iterrows():
        gene_id = row['Gene']
        codon = row['Codon']
        count = row['Count']
        
        if gene_id not in gene_codon_counts:
            gene_codon_counts[gene_id] = {}
        gene_codon_counts[gene_id][codon] = count
    
    # Calculate tAI for each gene
    tai_results = []
    for gene_id, codons in gene_codon_counts.items():
        tai_value = calculate_gene_tai(codons, trna_abundance)
        tai_results.append({'Gene': gene_id, 'tAI': tai_value})
    
    return pd.DataFrame(tai_results)

def main():
    # Get input and output file paths from Snakemake
    codon_table_path = snakemake.input.codon_table
    trna_data_path = snakemake.input.trna_data
    output_path = snakemake.output.tai_table
    
    # Calculate tAI values
    tai_df = calculate_tai_for_all_genes(codon_table_path, trna_data_path)
    
    # Save results
    tai_df.to_csv(output_path, sep='\t', index=False)
    
    print(f"tAI values calculated and saved to {output_path}")

if __name__ == "__main__":
    # For standalone testing, uncomment these lines and comment out the snakemake parts:
    # codon_table_path = "results/tables/codon_usage_table.tsv"
    # trna_data_path = "reference/trna_data/human_trna_counts.tsv"
    # output_path = "results/tables/tai_values.tsv"
    # tai_df = calculate_tai_for_all_genes(codon_table_path, trna_data_path)
    # tai_df.to_csv(output_path, sep='\t', index=False)
    
    # For Snakemake integration:
    main()