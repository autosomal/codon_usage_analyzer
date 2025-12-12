#!/usr/bin/env python3
"""
Calculate Relative Synonymous Codon Usage (RSCU)
"""

import pandas as pd
from collections import defaultdict
from Bio.Data import CodonTable

def get_genetic_code(table_id=1):
    """Get genetic code table from BioPython"""
    try:
        return CodonTable.unambiguous_dna_by_id[table_id]
    except KeyError:
        print(f"Warning: Genetic code table {table_id} not found. Using standard table (1).")
        return CodonTable.unambiguous_dna_by_id[1]

def calculate_rscu(codon_counts, genetic_code):
    """Calculate Relative Synonymous Codon Usage (RSCU)"""
    # Group codons by amino acid
    aa_codons = defaultdict(list)
    for codon, aa in genetic_code.forward_table.items():
        aa_codons[aa].append(codon)
    
    # Add stop codons
    aa_codons['*'] = genetic_code.stop_codons
    
    # Calculate RSCU
    rscu = {}
    total_codons = sum(codon_counts.values())
    
    for aa, codons in aa_codons.items():
        # Get counts for all codons of this amino acid
        aa_counts = {codon: codon_counts.get(codon, 0) for codon in codons}
        total_aa = sum(aa_counts.values())
        
        if total_aa == 0:
            # No usage for this amino acid
            for codon in codons:
                rscu[codon] = 0.0
            continue
        
        # Number of synonymous codons
        n_codons = len([c for c in codons if codon_counts.get(c, 0) > 0])
        if n_codons == 0:
            n_codons = len(codons)
        
        # Calculate RSCU for each codon
        for codon, count in aa_counts.items():
            if count == 0:
                rscu[codon] = 0.0
            else:
                rscu[codon] = count / (total_aa / n_codons)
    
    return rscu

def main():
    # Get input/output from Snakemake
    codon_table_file = snakemake.input[0]
    output_file = snakemake.output[0]
    
    # Get genetic code configuration
    genetic_code_config = snakemake.params.genetic_code
    
    # Get genetic code
    if 'table' in genetic_code_config:
        genetic_code = get_genetic_code(genetic_code_config['table'])
    else:
        genetic_code = get_genetic_code(1)  # Default to standard table
    
    # Read codon usage table
    codon_df = pd.read_csv(codon_table_file, sep='\t')
    
    # Create codon counts dictionary
    codon_counts = dict(zip(codon_df['codon'], codon_df['count']))
    
    # Calculate RSCU
    rscu_values = calculate_rscu(codon_counts, genetic_code)
    
    # Create RSCU DataFrame
    rscu_data = []
    for codon, rscu in rscu_values.items():
        aa = genetic_code.forward_table.get(codon, '*' if codon in genetic_code.stop_codons else 'X')
        count = codon_counts.get(codon, 0)
        rscu_data.append({
            'codon': codon,
            'amino_acid': aa,
            'count': count,
            'rscu': rscu,
            'rscu_category': 'preferred' if rscu > 1.2 else 'avoided' if rscu < 0.8 else 'neutral'
        })
    
    rscu_df = pd.DataFrame(rscu_data)
    
    # Add additional metrics
    rscu_df['frequency'] = rscu_df['count'] / rscu_df['count'].sum() if rscu_df['count'].sum() > 0 else 0
    
    # Save to output file
    rscu_df.to_csv(output_file, sep='\t', index=False)
    
    print(f"RSCU calculation completed successfully")
    print(f"Total codons analyzed: {rscu_df['count'].sum()}")
    print(f"Preferred codons: {len(rscu_df[rscu_df['rscu_category'] == 'preferred'])}")
    print(f"Avoided codons: {len(rscu_df[rscu_df['rscu_category'] == 'avoided'])}")
    print(f"Neutral codons: {len(rscu_df[rscu_df['rscu_category'] == 'neutral'])}")
    print(f"Output saved to: {output_file}")

if __name__ == "__main__":
    main()