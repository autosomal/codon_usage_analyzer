#!/usr/bin/env python3
"""
Identify preferred and avoided codons based on RSCU values
"""

import pandas as pd
from collections import defaultdict

def main():
    # Get input/output from Snakemake
    rscu_file = snakemake.input[0]
    preferred_output = snakemake.output.preferred_codons
    avoided_output = snakemake.output.avoided_codons
    
    # Get parameters
    preferred_threshold = snakemake.params.preferred_threshold
    avoided_threshold = snakemake.params.avoided_threshold
    
    # Read RSCU table
    rscu_df = pd.read_csv(rscu_file, sep='\t')
    
    # Group codons by amino acid
    aa_groups = rscu_df.groupby('amino_acid')
    
    # Identify preferred and avoided codons for each amino acid
    preferred_codons = []
    avoided_codons = []
    
    for aa, group in aa_groups:
        # Skip stop codons and unknown amino acids
        if aa == '*' or aa == 'X':
            continue
            
        # Get codons with sufficient counts
        valid_codons = group[group['count'] > 0].copy()
        
        if len(valid_codons) < 2:
            # Not enough codons to determine preference
            continue
            
        # Sort by RSCU
        sorted_codons = valid_codons.sort_values('rscu', ascending=False)
        
        # Identify preferred codons (highest RSCU)
        max_rscu = sorted_codons['rscu'].max()
        preferred = sorted_codons[sorted_codons['rscu'] >= preferred_threshold * max_rscu]
        
        # Identify avoided codons (lowest RSCU)
        min_rscu = sorted_codons['rscu'].min()
        avoided = sorted_codons[sorted_codons['rscu'] <= avoided_threshold * max_rscu]
        
        # Add to results
        for _, codon in preferred.iterrows():
            preferred_codons.append({
                'amino_acid': aa,
                'codon': codon['codon'],
                'rscu': codon['rscu'],
                'count': codon['count'],
                'frequency': codon['frequency'],
                'preference_rank': len(preferred_codons) % len(preferred) + 1
            })
        
        for _, codon in avoided.iterrows():
            avoided_codons.append({
                'amino_acid': aa,
                'codon': codon['codon'],
                'rscu': codon['rscu'],
                'count': codon['count'],
                'frequency': codon['frequency'],
                'avoidance_rank': len(avoided_codons) % len(avoided) + 1
            })
    
    # Also identify globally preferred/avoided based on absolute thresholds
    global_preferred = rscu_df[(rscu_df['rscu'] > preferred_threshold) & 
                              (rscu_df['amino_acid'] != '*') & 
                              (rscu_df['amino_acid'] != 'X')]
    
    global_avoided = rscu_df[(rscu_df['rscu'] < avoided_threshold) & 
                            (rscu_df['amino_acid'] != '*') & 
                            (rscu_df['amino_acid'] != 'X') &
                            (rscu_df['count'] > 0)]
    
    # Combine results
    preferred_df = pd.DataFrame(preferred_codons)
    global_preferred_df = global_preferred[['amino_acid', 'codon', 'rscu', 'count', 'frequency']].copy()
    global_preferred_df['type'] = 'global'
    
    avoided_df = pd.DataFrame(avoided_codons)
    global_avoided_df = global_avoided[['amino_acid', 'codon', 'rscu', 'count', 'frequency']].copy()
    global_avoided_df['type'] = 'global'
    
    # Save to files
    preferred_df.to_csv(preferred_output, sep='\t', index=False)
    avoided_df.to_csv(avoided_output, sep='\t', index=False)
    
    print(f"Preferred codons identified: {len(preferred_df)}")
    print(f"Avoided codons identified: {len(avoided_df)}")
    print(f"Global preferred codons: {len(global_preferred_df)}")
    print(f"Global avoided codons: {len(global_avoided_df)}")
    print(f"Preferred codons saved to: {preferred_output}")
    print(f"Avoided codons saved to: {avoided_output}")

if __name__ == "__main__":
    main()