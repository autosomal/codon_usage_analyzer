#!/usr/bin/env python3
"""
Parse CodonW output file into standardized format
"""

import os
import re
import pandas as pd
from collections import defaultdict

def parse_codonw_block(line):
    """Parse a single line from CodonW block format"""
    # Split line into parts, handling multiple spaces
    parts = re.split(r'\s+', line.strip())
    
    # Extract codon, amino acid, count, and frequency
    codon_data = []
    i = 0
    while i < len(parts):
        # Check if this part is an amino acid (3 letters)
        if len(parts[i]) == 3 and parts[i].isupper() and not parts[i].isdigit():
            aa = parts[i]
            i += 1
            # Next part should be the codon
            if i < len(parts) and len(parts[i]) == 3 and parts[i].isupper():
                codon = parts[i]
                i += 1
                # Next parts should be count and frequency
                if i + 1 < len(parts) and parts[i].isdigit() and '.' in parts[i+1]:
                    count = int(parts[i])
                    freq = float(parts[i+1])
                    i += 2
                    codon_data.append({
                        'codon': codon,
                        'amino_acid': aa,
                        'count': count,
                        'frequency': freq
                    })
                else:
                    # Skip if not a valid count/frequency pair
                    i += 1
            else:
                # Skip if not a valid codon
                i += 1
        else:
            # Skip non-amino acid parts
            i += 1
    
    return codon_data

def main():
    # Get input/output from Snakemake
    input_file = snakemake.input[0]
    output_file = snakemake.output[0]
    
    # Read input file
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    # Parse each line
    all_codons = []
    
    for line_num, line in enumerate(lines, 1):
        line = line.strip()
        if not line:
            continue
            
        # Skip header lines if present
        if line.startswith(('Phe', 'Leu', 'Ile', 'Val', 'Ser', 'Pro', 'Thr', 'Ala',
                           'Tyr', 'His', 'Gln', 'Asn', 'Lys', 'Asp', 'Glu', 'Cys',
                           'Trp', 'Arg', 'Ser', 'Arg', 'Gly')):
            codon_data = parse_codonw_block(line)
            all_codons.extend(codon_data)
    
    # Create DataFrame
    if all_codons:
        df = pd.DataFrame(all_codons)
        
        # Calculate additional metrics
        df['relative_frequency'] = df['frequency'] / df['frequency'].sum()
        
        # Save to output file
        df.to_csv(output_file, sep='\t', index=False)
        
        print(f"Successfully parsed CodonW file")
        print(f"Total codons parsed: {len(df)}")
        print(f"Output saved to: {output_file}")
    else:
        print("No codon data found in input file")
        # Create empty file
        open(output_file, 'w').close()

if __name__ == "__main__":
    main()