#!/usr/bin/env python3
"""
Calculate Effective Number of Codons (ENC) for genes
Based on Wright (1990) Gene
"""

import pandas as pd
import numpy as np
import sys
import os
from scipy.stats import chi2

def load_codon_usage(codon_table_path):
    """
    Load codon usage table from file
    Expected columns: Gene, Codon, Count, AminoAcid
    """
    df = pd.read_csv(codon_table_path, sep='\t')
    return df

def calculate_fi_for_aa(codon_counts_for_aa):
    """
    Calculate Fi (frequency of each codon) for amino acids with multiple codons
    Fi = (observed count of codon i - 1) / (total count of codons for this amino acid - number of synonymous codons)
    Actually, Fi is the frequency of codon i in the synonymous family
    """
    if len(codon_counts_for_aa) <= 1:
        # For amino acids with only one codon, no calculation needed
        return []
    
    total_count = sum(codon_counts_for_aa.values())
    if total_count == 0:
        return []
    
    # Calculate frequency of each codon
    frequencies = [count / total_count for count in codon_counts_for_aa.values()]
    
    return frequencies

def calculate_enc_for_synonymous_family(frequencies):
    """
    Calculate ENC for a single synonymous codon family
    ENC = 1 / sum(fi^2) where fi is frequency of codon i
    """
    if len(frequencies) == 0:
        return 0
    
    # Calculate sum of squared frequencies
    sum_sq_freq = sum(f**2 for f in frequencies)
    
    if sum_sq_freq == 0:
        return len(frequencies)  # Return number of codons if all have zero frequency
    
    # ENC for this family
    enc_family = 1 / sum_sq_freq
    return enc_family

def calculate_enc_for_gene(gene_data):
    """
    Calculate ENC for a single gene
    Wright's formula: ENC = 2 + (9/(1/S2) + 1/(1/S3) + 5/(1/S4) + 3/(1/S6))
    where S2, S3, S4, S6 are synonymous codon families with 2, 3, 4, and 6 codons respectively
    """
    # Group codons by amino acid
    aa_groups = gene_data.groupby('AminoAcid')
    
    # Calculate F statistics for different degenerate codon sets
    # 2-fold degenerate (ending in T/Y or A/T, etc.)
    two_fold_codons = {}
    three_fold_codons = {}
    four_fold_codons = {}
    six_fold_codons = {}
    
    # Define amino acids by number of codons
    aa_2_fold = ['Phe', 'Tyr', 'His', 'Asn', 'Asp', 'Cys', 'Gln', 'Lys', 'Glu']  # 2 codons each
    aa_3_fold = ['Ile']  # 3 codons (ATT, ATC, ATA)
    aa_4_fold = ['Leu', 'Val', 'Ser', 'Pro', 'Thr', 'Ala', 'Gly']  # multiple codons, some with 4-fold degeneracy
    aa_6_fold = ['Arg', 'Leu', 'Ser']  # these have 6 codons each
    
    # For Wright's ENC calculation, we need to calculate F2, F3, F4, F6 (mean frequency of codons)
    # F = (sum of ni*fi^2 - 1) / (n - 1) where ni is count of codon i, fi is frequency, n is total codons
    
    enc_values = []
    
    for aa, group in aa_groups:
        # Skip amino acids with only one codon (Met, Trp, Stop)
        if len(group) <= 1:
            continue
            
        # Calculate frequencies for this amino acid's codons
        codon_counts = dict(zip(group['Codon'], group['Count']))
        frequencies = calculate_fi_for_aa(codon_counts)
        
        if frequencies:
            enc_aa = calculate_enc_for_synonymous_family(frequencies)
            enc_values.append(enc_aa)
    
    # Wright's original method for overall ENC
    # This is a simplified version - the full method considers different codon degeneracies
    # A more accurate implementation would follow Wright's formula more precisely
    
    # Calculate F statistics for different codon degeneracies
    # We'll implement Wright's approach more directly:
    
    # Get total counts for each amino acid
    aa_totals = gene_data.groupby('AminoAcid')['Count'].sum()
    
    # Calculate the F statistics for each amino acid group (2, 3, 4, 6 fold degenerate)
    f_stats = []
    for aa, group in aa_groups:
        if len(group) <= 1:  # Skip single-codon amino acids
            continue
            
        counts = group['Count'].values
        total = sum(counts)
        
        if total == 0:
            continue
            
        # Calculate frequency for each codon
        freqs = counts / total
        
        # Calculate F (frequency-based statistic) for this amino acid
        # F = sum(fi^2) where fi is frequency of codon i
        f_sum = sum(f**2 for f in freqs)
        
        # Expected F if all codons equally used = 1/ni (where ni is number of synonymous codons)
        expected_f = 1.0 / len(group)
        
        # Calculate F statistic: (observed F - expected F) / (max possible F - expected F)
        # Max possible F is 1 (when one codon used 100% of time)
        if (1 - expected_f) != 0:
            f_stat = (f_sum - expected_f) / (1 - expected_f)
        else:
            f_stat = 0
            
        f_stats.append((f_sum, len(group), total))
    
    # Wright's ENC formula
    # ENC = 2 + (9/S2) + (1/S3) + (5/S4) + (3/S6) where Si = mean F for i-fold degenerate sites
    # Actually, it's: ENC = 2 + (9/F2) + (1/F3) + (5/F4) + (3/F6) where Fi is mean F for i-fold system
    
    # Group F values by number of codons per amino acid
    f_by_codons = {}
    for f_sum, n_codons, total in f_stats:
        if n_codons not in f_by_codons:
            f_by_codons[n_codons] = []
        # Calculate F for this codon family
        if total > n_codons:  # Need at least 2 codons to calculate F
            # F = (sum of ni*pi^2 - 1) / (n - 1) where pi = freq of codon i, n = total codons
            # Actually for ENC: F = sum of (ni * (pi)^2) where ni is count of codon i, pi is freq
            # Then F = (sum of (ni * pi^2) - 1) / (n - 1) for expected heterozygosity
            # Simpler: F = sum of (ni * pi^2) where pi is relative frequency in the codon family
            f_by_codons[n_codons].append(f_sum)
    
    # Calculate mean F for each codon degeneracy
    mean_f_2 = np.mean(f_by_codons.get(2, [])) if f_by_codons.get(2) else None
    mean_f_3 = np.mean(f_by_codons.get(3, [])) if f_by_codons.get(3) else None
    mean_f_4 = np.mean(f_by_codons.get(4, [])) if f_by_codons.get(4) else None
    mean_f_6 = np.mean(f_by_codons.get(6, [])) if f_by_codons.get(6) else None
    
    # Apply Wright's formula
    enc = 2  # Start with 2 (for 2-codon families)
    
    if mean_f_2 and mean_f_2 != 0:
        enc += 9 / (mean_f_2 * len(f_by_codons.get(2, [])) / len(f_by_codons.get(2, []))) if len(f_by_codons.get(2, [])) > 0 else 9 / mean_f_2
    else:
        enc += 9  # If no 2-fold codons or F2 is 0, max contribution
    
    if mean_f_3 and mean_f_3 != 0:
        enc += 1 / mean_f_3
    else:
        enc += 1  # If no 3-fold codons or F3 is 0, max contribution
    
    if mean_f_4 and mean_f_4 != 0:
        enc += 5 / mean_f_4
    else:
        enc += 5  # If no 4-fold codons or F4 is 0, max contribution
    
    if mean_f_6 and mean_f_6 != 0:
        enc += 3 / mean_f_6
    else:
        enc += 3  # If no 6-fold codons or F6 is 0, max contribution
    
    # Cap ENC at 61 (maximum possible for 61 sense codons)
    enc = min(enc, 61)
    
    # A simpler implementation of ENC calculation
    # ENC = sum of (ni * Fi) for each codon family i
    # Where ni is number of codons in family i, Fi = 1/sum(pi^2) for codon frequencies pi
    
    # Actually, let's implement the correct formula:
    # ENC = sum over all codon families of: ni / sum(pi^2) 
    # where ni is number of codons in family i, pi is frequency of codon i in family
    
    # Simpler approach following the definition more directly:
    total_enc_contrib = 0
    family_count = 0
    
    for aa, group in aa_groups:
        if len(group) <= 1:  # Skip single-codon amino acids
            continue
            
        counts = group['Count'].values
        total = sum(counts)
        
        if total == 0:
            continue
            
        # Calculate frequency for each codon in this family
        freqs = counts / total
        
        # Calculate sum of squared frequencies
        sum_sq_freq = sum(f**2 for f in freqs)
        
        if sum_sq_freq > 0:
            # ENC contribution for this amino acid = 1 / sum of squared frequencies
            enc_contrib = 1 / sum_sq_freq
            total_enc_contrib += enc_contrib
            family_count += 1
    
    # Average ENC across amino acid families (this is a simplified approach)
    if family_count > 0:
        enc_simple = total_enc_contrib / family_count
    else:
        enc_simple = 0
    
    # Return the more accurate calculation based on Wright's method
    # Let's implement a cleaner version
    return calculate_enc_simple(gene_data)

def calculate_enc_simple(gene_data):
    """
    Simplified ENC calculation following Wright (1990)
    """
    # Group by amino acid
    aa_groups = gene_data.groupby('AminoAcid')
    
    # Calculate ENC using Wright's formula
    # ENC = 2 + (9/F2) + (1/F3) + (5/F4) + (3/F6)
    # where F2, F3, F4, F6 are the mean F statistics for 2, 3, 4, and 6-codon amino acids
    
    # Define amino acids by number of codons
    aa_2_codon = {'Phe', 'Cys', 'Tyr', 'His', 'Gln', 'Asn', 'Lys', 'Asp', 'Glu'}  # 2 codons each
    aa_3_codon = {'Ile'}  # 3 codons
    aa_4_codon = {'Leu', 'Val', 'Ser', 'Pro', 'Thr', 'Ala', 'Gly'}  # 4 codons (some have more due to alternative codons)
    aa_6_codon = {'Arg'}  # 6 codons
    
    # Calculate F statistics for each amino acid
    f2_values = []
    f3_values = []
    f4_values = []
    f6_values = []
    
    for aa, group in aa_groups:
        if len(group) <= 1:  # Skip Met, Trp, and stop codons
            continue
            
        counts = group['Count'].values
        total = sum(counts)
        
        if total == 0:
            continue
            
        # Calculate frequencies
        freqs = counts / total
        
        # Calculate F (homozygosity) = sum of squared frequencies
        f_stat = sum(freq**2 for freq in freqs)
        
        # Assign to appropriate group
        n_codons = len(group)
        if n_codons == 2:
            f2_values.append(f_stat)
        elif n_codons == 3:
            f3_values.append(f_stat)
        elif n_codons == 4:
            f4_values.append(f_stat)
        elif n_codons == 6:
            f6_values.append(f_stat)
    
    # Calculate mean F values
    mean_f2 = np.mean(f2_values) if f2_values else None
    mean_f3 = np.mean(f3_values) if f3_values else None
    mean_f4 = np.mean(f4_values) if f4_values else None
    mean_f6 = np.mean(f6_values) if f6_values else None
    
    # Calculate ENC using Wright's formula
    enc = 2  # Start with 2-codon contribution
    
    if mean_f2 and mean_f2 != 0:
        # For 2-codon families: 9/F2, but weighted by number of such families
        enc += 9 / mean_f2
    else:
        enc += 9  # Maximum contribution if no data or F2=0
    
    if mean_f3 and mean_f3 != 0:
        enc += 1 / mean_f3
    elif mean_f3 is not None:  # If we have 3-codon families but F3 is 0
        enc += 1
    
    if mean_f4 and mean_f4 != 0:
        enc += 5 / mean_f4
    elif mean_f4 is not None:  # If we have 4-codon families but F4 is 0
        enc += 5
    
    if mean_f6 and mean_f6 != 0:
        enc += 3 / mean_f6
    elif mean_f6 is not None:  # If we have 6-codon families but F6 is 0
        enc += 3
    
    # Cap at theoretical maximum
    enc = min(enc, 61.0)
    
    return enc

def calculate_enc_for_all_genes(codon_table_path):
    """
    Calculate ENC for all genes in the codon usage table
    """
    # Load data
    codon_df = load_codon_usage(codon_table_path)
    
    # Group by gene
    gene_groups = codon_df.groupby('Gene')
    
    # Calculate ENC for each gene
    enc_results = []
    for gene_id, gene_data in gene_groups:
        enc_value = calculate_enc_simple(gene_data)
        enc_results.append({'Gene': gene_id, 'ENC': enc_value})
    
    return pd.DataFrame(enc_results)

def main():
    # Get input and output file paths from Snakemake
    codon_table_path = snakemake.input.codon_table
    output_path = snakemake.output.enc_table
    
    # Calculate ENC values
    enc_df = calculate_enc_for_all_genes(codon_table_path)
    
    # Save results
    enc_df.to_csv(output_path, sep='\t', index=False)
    
    print(f"ENC values calculated and saved to {output_path}")

if __name__ == "__main__":
    # For Snakemake integration:
    main()