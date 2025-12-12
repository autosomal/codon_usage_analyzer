#!/usr/bin/env python3
"""
Codon pair analysis to identify codon pair usage bias
"""

import pandas as pd
import numpy as np
import sys
import os
from Bio import SeqIO
from collections import defaultdict, Counter

def extract_codons(sequence):
    """
    Extract codons from a DNA sequence
    """
    codons = []
    # Ensure sequence length is divisible by 3
    seq_len = len(sequence) - (len(sequence) % 3)
    for i in range(0, seq_len, 3):
        codon = sequence[i:i+3].upper()
        if len(codon) == 3 and all(base in 'ATCG' for base in codon):
            codons.append(codon)
    return codons

def analyze_codon_pairs_from_fasta(fasta_file):
    """
    Analyze codon pairs from FASTA file
    """
    gene_codon_pairs = {}
    dipeptide_counts = Counter()
    
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq)
        codons = extract_codons(sequence)
        
        # Get codon pairs (adjacent codons)
        codon_pairs = []
        for i in range(len(codons) - 1):
            pair = (codons[i], codons[i+1])
            codon_pairs.append(pair)
            
            # Also track dipeptides (could reveal bias)
            # This would require a genetic code mapping
            # For now, just store the codon pair
            
        gene_codon_pairs[record.id] = codon_pairs
        
        # Count dipeptide occurrences (conceptual - would need genetic code)
        for codon1, codon2 in codon_pairs:
            # In a real implementation, we'd translate to amino acids
            dipeptide = f"{codon1}-{codon2}"  # Placeholder
            dipeptide_counts[dipeptide] += 1
    
    return gene_codon_pairs, dipeptide_counts

def calculate_codon_pair_usage(gene_codon_pairs):
    """
    Calculate overall codon pair usage across all genes
    """
    pair_counts = Counter()
    total_pairs = 0
    
    for gene_id, pairs in gene_codon_pairs.items():
        for pair in pairs:
            pair_key = f"{pair[0]}-{pair[1]}"
            pair_counts[pair_key] += 1
            total_pairs += 1
    
    # Convert to DataFrame
    pairs_data = []
    for pair, count in pair_counts.items():
        freq = count / total_pairs if total_pairs > 0 else 0
        codons = pair.split('-')
        pairs_data.append({
            'Codon1': codons[0],
            'Codon2': codons[1],
            'Pair': pair,
            'Count': count,
            'Frequency': freq
        })
    
    return pd.DataFrame(pairs_data)

def calculate_codon_pair_bias(codon_pair_df):
    """
    Calculate basic codon pair bias metrics
    """
    # Calculate ratio of observed vs expected if independent
    # This is a simplified version - full analysis would be more complex
    
    # Count individual codons appearing in first and second position
    codon1_counts = codon_pair_df.groupby('Codon1')['Count'].sum()
    codon2_counts = codon_pair_df.groupby('Codon2')['Count'].sum()
    
    total_codon_pairs = codon_pair_df['Count'].sum()
    
    # Calculate expected counts if codons were independent
    bias_data = []
    for _, row in codon_pair_df.iterrows():
        observed = row['Count']
        expected = (codon1_counts[row['Codon1']] * codon2_counts[row['Codon2']]) / total_codon_pairs
        bias_score = observed / expected if expected > 0 else 0
        
        bias_data.append({
            'Pair': row['Pair'],
            'Observed': observed,
            'Expected': expected,
            'BiasScore': bias_score,
            'IsEnriched': bias_score > 1.2,  # Arbitrary threshold
            'IsDepleted': bias_score < 0.8   # Arbitrary threshold
        })
    
    return pd.DataFrame(bias_data)

def main():
    # Get input and output file paths from Snakemake
    fasta_file = snakemake.input.fasta
    codon_pairs_output = snakemake.output.codon_pairs
    dipeptide_output = snakemake.output.dipeptide_bias
    
    # Analyze codon pairs
    gene_codon_pairs, dipeptide_counts = analyze_codon_pairs_from_fasta(fasta_file)
    
    # Calculate codon pair usage
    codon_pair_df = calculate_codon_pair_usage(gene_codon_pairs)
    
    # Calculate codon pair bias
    bias_df = calculate_codon_pair_bias(codon_pair_df)
    
    # Save codon pair usage
    codon_pair_df.to_csv(codon_pairs_output, sep='\t', index=False)
    
    # Save dipeptide bias (simplified)
    dipeptide_df = pd.DataFrame([
        {'Dipeptide': dp, 'Count': count} 
        for dp, count in dipeptide_counts.items()
    ])
    dipeptide_df.to_csv(dipeptide_output, sep='\t', index=False)
    
    print(f"Codon pair analysis completed:")
    print(f"- Codon pair usage saved to {codon_pairs_output}")
    print(f"- Dipeptide bias saved to {dipeptide_output}")

if __name__ == "__main__":
    # For Snakemake integration:
    main()