#!/usr/bin/env python3
"""
Cross-Species Codon Usage Comparison
Compares codon usage patterns across multiple species
"""

import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram
import argparse
import sys
import os
import json

def load_codon_data(file_path, species_name):
    """
    Load codon usage data for a specific species
    """
    try:
        data = pd.read_csv(file_path, sep='\t')
        required_cols = ['gene_id', 'codon', 'frequency']
        if not all(col in data.columns for col in required_cols):
            raise ValueError(f"Required columns {required_cols} not found in {file_path}")
        
        # Add species column
        data['species'] = species_name
        return data
    except Exception as e:
        raise ValueError(f"Error loading codon data for {species_name}: {str(e)}")

def calculate_codon_frequencies_by_species(codon_data):
    """
    Calculate average codon frequencies per species
    """
    # Group by species and codon to get average frequencies
    avg_freq = codon_data.groupby(['species', 'codon'])['frequency'].mean().reset_index()
    
    # Pivot to wide format: rows=species, columns=codon, values=frequency
    freq_matrix = avg_freq.pivot(index='species', columns='codon', values='frequency')
    
    # Fill missing values with 0
    freq_matrix = freq_matrix.fillna(0)
    
    return freq_matrix

def calculate_distance_matrix(freq_matrix):
    """
    Calculate distance matrix between species based on codon usage
    """
    # Calculate Euclidean distance between species
    distances = pdist(freq_matrix.values, metric='euclidean')
    distance_matrix = squareform(distances)
    
    # Convert to DataFrame with proper labels
    species_names = freq_matrix.index.tolist()
    distance_df = pd.DataFrame(distance_matrix, index=species_names, columns=species_names)
    
    return distance_df

def calculate_correlation_matrix(freq_matrix):
    """
    Calculate correlation matrix between species based on codon usage
    """
    # Calculate correlation matrix
    correlation_df = freq_matrix.T.corr(method='pearson')
    return correlation_df

def calculate_codon_bias_metrics(freq_matrix):
    """
    Calculate various codon bias metrics for each species
    """
    metrics = []
    
    for species in freq_matrix.index:
        species_freqs = freq_matrix.loc[species]
        
        # Calculate Effective Number of Codons (simplified)
        # This is a basic approximation - a full ENC calculation would be more complex
        non_zero_freqs = species_freqs[species_freqs > 0]
        if len(non_zero_freqs) > 0:
            enc_approx = len(non_zero_freqs) / (species_freqs**2).sum() if (species_freqs**2).sum() > 0 else len(non_zero_freqs)
        else:
            enc_approx = 0
            
        # Calculate codon usage evenness (Shannon entropy-based)
        freqs_no_zero = non_zero_freqs[non_zero_freqs > 0] / non_zero_freqs.sum()
        entropy = -sum(f * np.log(f) for f in freqs_no_zero if f > 0)
        max_entropy = np.log(len(non_zero_freqs)) if len(non_zero_freqs) > 0 else 1
        evenness = entropy / max_entropy if max_entropy > 0 else 0
        
        metrics.append({
            'species': species,
            'mean_codon_frequency': species_freqs.mean(),
            'std_codon_frequency': species_freqs.std(),
            'num_codons_measured': len(non_zero_freqs),
            'approximate_enc': enc_approx,
            'codon_usage_evenness': evenness
        })
    
    return pd.DataFrame(metrics)

def main():
    parser = argparse.ArgumentParser(description='Compare codon usage across multiple species')
    parser.add_argument('input_files', nargs='+', help='Input codon usage tables (TSV) with format species:file')
    parser.add_argument('output_comparison', help='Output comparison table (TSV)')
    parser.add_argument('output_distance', help='Output distance matrix (TSV)')
    
    args = parser.parse_args()
    
    # Parse input files in format species:file
    species_files = []
    for arg in args.input_files:
        if ':' in arg:
            species, file_path = arg.split(':', 1)
            species_files.append((species, file_path))
        else:
            # If no species specified, extract from filename
            base_name = os.path.basename(arg)
            species = base_name.split('.')[0]  # Remove extension
            species_files.append((species, arg))
    
    # Validate input files
    for species, file_path in species_files:
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"Input file not found: {file_path} for species {species}")
    
    # Load all codon data
    print("Loading codon usage data for all species...")
    all_data = []
    for species, file_path in species_files:
        print(f"  Loading data for {species} from {file_path}")
        species_data = load_codon_data(file_path, species)
        all_data.append(species_data)
    
    # Combine all data
    combined_data = pd.concat(all_data, ignore_index=True)
    
    # Calculate codon frequency matrix
    print("Calculating codon frequency matrix...")
    freq_matrix = calculate_codon_frequencies_by_species(combined_data)
    
    # Calculate distance matrix
    print("Calculating distance matrix...")
    distance_matrix = calculate_distance_matrix(freq_matrix)
    
    # Calculate correlation matrix
    print("Calculating correlation matrix...")
    correlation_matrix = calculate_correlation_matrix(freq_matrix)
    
    # Calculate codon bias metrics
    print("Calculating codon bias metrics...")
    bias_metrics = calculate_codon_bias_metrics(freq_matrix)
    
    # Create comprehensive comparison table
    comparison_data = []
    for species in freq_matrix.index:
        species_row = {'species': species}
        
        # Add basic metrics
        metrics_row = bias_metrics[bias_metrics['species'] == species].iloc[0] if len(bias_metrics[bias_metrics['species'] == species]) > 0 else None
        if metrics_row is not None:
            species_row.update({
                'mean_codon_frequency': metrics_row['mean_codon_frequency'],
                'std_codon_frequency': metrics_row['std_codon_frequency'],
                'num_codons_measured': metrics_row['num_codons_measured'],
                'approximate_enc': metrics_row['approximate_enc'],
                'codon_usage_evenness': metrics_row['codon_usage_evenness']
            })
        
        # Add similarity scores to other species
        for other_species in freq_matrix.index:
            if species != other_species:
                similarity = correlation_matrix.loc[species, other_species]
                species_row[f'similarity_to_{other_species}'] = similarity
        
        comparison_data.append(species_row)
    
    comparison_df = pd.DataFrame(comparison_data)
    
    # Save results
    print(f"Saving comparison results to {args.output_comparison}...")
    comparison_df.to_csv(args.output_comparison, sep='\t', index=False)
    
    print(f"Saving distance matrix to {args.output_distance}...")
    distance_matrix.to_csv(args.output_distance, sep='\t')
    
    # Print summary
    print(f"\nCross-Species Comparison Summary:")
    print(f"Number of species compared: {len(freq_matrix.index)}")
    print(f"Number of codons analyzed: {len(freq_matrix.columns)}")
    print(f"Species analyzed: {', '.join(freq_matrix.index.tolist())}")
    
    print(f"\nTop species similarities:")
    # Get top similarities (excluding self-comparisons)
    upper_triangle = np.triu(correlation_matrix.values, k=1)
    upper_triangle[upper_triangle == 0] = np.nan  # Replace zeros with NaN for off-diagonal
    flat_similarities = upper_triangle.flatten()
    flat_similarities = flat_similarities[~np.isnan(flat_similarities)]
    
    if len(flat_similarities) > 0:
        top_similarities = sorted(flat_similarities, reverse=True)[:5]
        print("Top 5 species pair similarities:")
        for i, sim in enumerate(top_similarities):
            print(f"  {i+1}. {sim:.4f}")
    
    print(f"\nResults saved to:")
    print(f"  Comparison table: {args.output_comparison}")
    print(f"  Distance matrix: {args.output_distance}")

if __name__ == "__main__":
    main()