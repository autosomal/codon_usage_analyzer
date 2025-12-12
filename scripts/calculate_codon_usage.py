#!/usr/bin/env python3
"""
Calculate codon usage from FASTA file
"""

import os
import sys
import re
import pandas as pd
from collections import defaultdict
from Bio import SeqIO
from Bio.Data import CodonTable

def get_genetic_code(table_id=1):
    """Get genetic code table from BioPython"""
    try:
        return CodonTable.unambiguous_dna_by_id[table_id]
    except KeyError:
        print(f"Warning: Genetic code table {table_id} not found. Using standard table (1).")
        return CodonTable.unambiguous_dna_by_id[1]

def reverse_complement(seq):
    """Calculate reverse complement of DNA sequence"""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join([complement[base] for base in reversed(seq)])

def translate_codon(codon, genetic_code):
    """Translate codon to amino acid using given genetic code"""
    codon = codon.upper()
    if codon in genetic_code.stop_codons:
        return '*'
    elif codon in genetic_code.forward_table:
        return genetic_code.forward_table[codon]
    else:
        return 'X'  # Unknown amino acid

def calculate_codon_usage(sequence, genetic_code):
    """Calculate codon usage for a single sequence"""
    codon_counts = defaultdict(int)
    sequence = sequence.upper().replace('U', 'T')  # Convert RNA to DNA if needed
    
    # Process each codon in the sequence
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i+3]
        if len(codon) == 3 and 'N' not in codon:
            codon_counts[codon] += 1
    
    return codon_counts

def process_fasta_file(fasta_file, genetic_code):
    """Process FASTA file and calculate codon usage"""
    total_codons = defaultdict(int)
    gene_counts = []
    
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq)
        codon_counts = calculate_codon_usage(sequence, genetic_code)
        
        # Update total counts
        for codon, count in codon_counts.items():
            total_codons[codon] += count
        
        # Save gene-level counts
        gene_counts.append({
            'gene_id': record.id,
            'length': len(sequence),
            'codon_counts': codon_counts
        })
    
    return total_codons, gene_counts

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
    fasta_file = snakemake.input[0]
    codon_table_output = snakemake.output.codon_table
    amino_acid_output = snakemake.output.amino_acid_table
    
    # Get parameters
    genetic_code_config = snakemake.params.genetic_code
    calculate_rscu_flag = snakemake.params.calculate_rscu
    calculate_cai_flag = snakemake.params.calculate_cai
    calculate_enc_flag = snakemake.params.calculate_enc
    
    # Get genetic code
    if 'table' in genetic_code_config:
        genetic_code = get_genetic_code(genetic_code_config['table'])
    elif 'custom_table' in genetic_code_config:
        # Create custom genetic code
        custom_table = genetic_code_config['custom_table']
        # This would need more complex handling
        genetic_code = get_genetic_code(1)  # Default to standard table for now
    else:
        genetic_code = get_genetic_code(1)  # Default to standard table
    
    print(f"Using genetic code table: {genetic_code}")
    
    # Process FASTA file
    print(f"Processing FASTA file: {fasta_file}")
    total_codons, gene_counts = process_fasta_file(fasta_file, genetic_code)
    
    total_count = sum(total_codons.values())
    print(f"Total codons analyzed: {total_count}")
    
    # Create codon usage table
    codon_data = []
    for codon, count in sorted(total_codons.items()):
        aa = translate_codon(codon, genetic_code)
        frequency = count / total_count if total_count > 0 else 0
        
        codon_data.append({
            'codon': codon,
            'amino_acid': aa,
            'count': count,
            'frequency': frequency,
            'percentage': frequency * 100
        })
    
    # Calculate RSCU if requested
    if calculate_rscu_flag and total_count > 0:
        rscu_values = calculate_rscu(total_codons, genetic_code)
        for entry in codon_data:
            entry['rscu'] = rscu_values.get(entry['codon'], 0.0)
    
    # Create DataFrame and save
    codon_df = pd.DataFrame(codon_data)
    codon_df.to_csv(codon_table_output, sep='\t', index=False)
    
    # Create amino acid usage table
    aa_data = defaultdict(lambda: {'count': 0, 'frequency': 0.0, 'percentage': 0.0})
    
    for entry in codon_data:
        aa = entry['amino_acid']
        aa_data[aa]['count'] += entry['count']
        aa_data[aa]['frequency'] += entry['frequency']
        aa_data[aa]['percentage'] += entry['percentage']
    
    aa_df = pd.DataFrame([{'amino_acid': aa, **data} for aa, data in aa_data.items()])
    aa_df.to_csv(amino_acid_output, sep='\t', index=False)
    
    print(f"Codon usage table saved to: {codon_table_output}")
    print(f"Amino acid usage table saved to: {amino_acid_output}")
    print(f"Total codons: {total_count}")
    print(f"Unique codons: {len(codon_data)}")
    print(f"Unique amino acids: {len(aa_data)}")

if __name__ == "__main__":
    main()