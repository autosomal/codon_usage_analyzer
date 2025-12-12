#!/usr/bin/env python3
"""
Generate Comprehensive Analysis Report
Creates a detailed HTML report summarizing all codon usage analysis results
"""

import pandas as pd
import numpy as np
import os
import sys
from datetime import datetime
import argparse
import json

def load_analysis_results(output_dir):
    """
    Load all analysis results from the output directory
    """
    results = {}
    
    # Load main codon usage tables
    codon_usage_file = os.path.join(output_dir, "tables", "codon_usage_table.tsv")
    if os.path.exists(codon_usage_file):
        results['codon_usage'] = pd.read_csv(codon_usage_file, sep='\t')
    
    cai_file = os.path.join(output_dir, "tables", "cai_values.tsv")
    if os.path.exists(cai_file):
        results['cai_values'] = pd.read_csv(cai_file, sep='\t')
    
    enc_file = os.path.join(output_dir, "tables", "enc_values.tsv")
    if os.path.exists(enc_file):
        results['enc_values'] = pd.read_csv(enc_file, sep='\t')
    
    rscu_file = os.path.join(output_dir, "tables", "rscu_values.tsv")
    if os.path.exists(rscu_file):
        results['rscu_values'] = pd.read_csv(rscu_file, sep='\t')
    
    # Load other metrics if they exist
    tai_file = os.path.join(output_dir, "tables", "tai_values.tsv")
    if os.path.exists(tai_file):
        results['tai_values'] = pd.read_csv(tai_file, sep='\t')
    
    correlation_file = os.path.join(output_dir, "tables", "codon_expression_correlation.tsv")
    if os.path.exists(correlation_file):
        results['correlation'] = pd.read_csv(correlation_file, sep='\t')
    
    return results

def generate_statistics_summary(results):
    """
    Generate summary statistics from analysis results
    """
    summary = {}
    
    if 'cai_values' in results:
        cai_data = results['cai_values']
        summary['cai_stats'] = {
            'count': len(cai_data),
            'mean': cai_data['cai_value'].mean() if 'cai_value' in cai_data.columns else 0,
            'std': cai_data['cai_value'].std() if 'cai_value' in cai_data.columns else 0,
            'min': cai_data['cai_value'].min() if 'cai_value' in cai_data.columns else 0,
            'max': cai_data['cai_value'].max() if 'cai_value' in cai_data.columns else 0
        }
    
    if 'enc_values' in results:
        enc_data = results['enc_values']
        summary['enc_stats'] = {
            'count': len(enc_data),
            'mean': enc_data['enc_value'].mean() if 'enc_value' in enc_data.columns else 0,
            'std': enc_data['enc_value'].std() if 'enc_value' in enc_data.columns else 0,
            'min': enc_data['enc_value'].min() if 'enc_value' in enc_data.columns else 0,
            'max': enc_data['enc_value'].max() if 'enc_value' in enc_data.columns else 0
        }
    
    if 'rscu_values' in results:
        rscu_data = results['rscu_values']
        summary['rscu_stats'] = {
            'unique_codons': rscu_data['codon'].nunique() if 'codon' in rscu_data.columns else 0,
            'mean_rscu': rscu_data['rscu_value'].mean() if 'rscu_value' in rscu_data.columns else 0,
            'std_rscu': rscu_data['rscu_value'].std() if 'rscu_value' in rscu_data.columns else 0
        }
    
    if 'tai_values' in results:
        tai_data = results['tai_values']
        summary['tai_stats'] = {
            'count': len(tai_data),
            'mean': tai_data['tai_value'].mean() if 'tai_value' in tai_data.columns else 0,
            'std': tai_data['tai_value'].std() if 'tai_value' in tai_data.columns else 0,
            'min': tai_data['tai_value'].min() if 'tai_value' in tai_data.columns else 0,
            'max': tai_data['tai_value'].max() if 'tai_value' in tai_data.columns else 0
        }
    
    return summary

def generate_html_report(results, summary, output_file, config_info):
    """
    Generate comprehensive HTML report
    """
    html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>Comprehensive Codon Usage Analysis Report</title>
    <meta charset="UTF-8">
    <style>
        body {{
            font-family: Arial, sans-serif;
            margin: 40px;
            background-color: #f9f9f9;
            color: #333;
        }}
        .header {{
            background-color: #2c3e50;
            color: white;
            padding: 20px;
            border-radius: 5px;
            margin-bottom: 20px;
        }}
        .section {{
            background-color: white;
            padding: 20px;
            margin-bottom: 20px;
            border-radius: 5px;
            box-shadow: 0 2px 5px rgba(0,0,0,0.1);
        }}
        .stats-table {{
            width: 100%;
            border-collapse: collapse;
            margin-top: 10px;
        }}
        .stats-table th, .stats-table td {{
            border: 1px solid #ddd;
            padding: 8px;
            text-align: left;
        }}
        .stats-table th {{
            background-color: #f2f2f2;
        }}
        .plot-container {{
            text-align: center;
            margin: 20px 0;
        }}
        .plot {{
            max-width: 100%;
            height: auto;
        }}
        h1, h2, h3 {{
            color: #2c3e50;
        }}
    </style>
</head>
<body>
    <div class="header">
        <h1>Comprehensive Codon Usage Analysis Report</h1>
        <p><strong>Analysis Date:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        <p><strong>Analysis Pipeline:</strong> Codon Usage Analyzer v1.0</p>
    </div>
    
    <div class="section">
        <h2>Analysis Configuration</h2>
        <table class="stats-table">
            <tr><th>Parameter</th><th>Value</th></tr>
            <tr><td>Input Type</td><td>{config_info.get('input_type', 'N/A')}</td></tr>
            <tr><td>Genetic Code</td><td>{config_info.get('genetic_code', 'N/A')}</td></tr>
            <tr><td>Calculate RSCU</td><td>{config_info.get('calculate_rscu', 'N/A')}</td></tr>
            <tr><td>Calculate CAI</td><td>{config_info.get('calculate_cai', 'N/A')}</td></tr>
            <tr><td>Calculate ENC</td><td>{config_info.get('calculate_enc', 'N/A')}</td></tr>
            <tr><td>Calculate tAI</td><td>{config_info.get('calculate_tai', 'N/A')}</td></tr>
        </table>
    </div>
    
    <div class="section">
        <h2>Summary Statistics</h2>
        <table class="stats-table">
            <tr><th>Metric</th><th>Count/Value</th><th>Mean</th><th>Std Dev</th><th>Min</th><th>Max</th></tr>
    """
    
    # Add CAI statistics
    if 'cai_stats' in summary:
        stats = summary['cai_stats']
        html_content += f"""
            <tr>
                <td>CAI (Codon Adaptation Index)</td>
                <td>{stats['count']}</td>
                <td>{stats['mean']:.4f}</td>
                <td>{stats['std']:.4f}</td>
                <td>{stats['min']:.4f}</td>
                <td>{stats['max']:.4f}</td>
            </tr>
        """
    
    # Add ENC statistics
    if 'enc_stats' in summary:
        stats = summary['enc_stats']
        html_content += f"""
            <tr>
                <td>ENC (Effective Number of Codons)</td>
                <td>{stats['count']}</td>
                <td>{stats['mean']:.4f}</td>
                <td>{stats['std']:.4f}</td>
                <td>{stats['min']:.4f}</td>
                <td>{stats['max']:.4f}</td>
            </tr>
        """
    
    # Add RSCU statistics
    if 'rscu_stats' in summary:
        stats = summary['rscu_stats']
        html_content += f"""
            <tr>
                <td>RSCU (Relative Synonymous Codon Usage)</td>
                <td>{stats['unique_codons']} codons</td>
                <td>{stats['mean_rscu']:.4f}</td>
                <td>{stats['std_rscu']:.4f}</td>
                <td>N/A</td>
                <td>N/A</td>
            </tr>
        """
    
    # Add tAI statistics
    if 'tai_stats' in summary:
        stats = summary['tai_stats']
        html_content += f"""
            <tr>
                <td>tAI (tRNA Adaptation Index)</td>
                <td>{stats['count']}</td>
                <td>{stats['mean']:.4f}</td>
                <td>{stats['std']:.4f}</td>
                <td>{stats['min']:.4f}</td>
                <td>{stats['max']:.4f}</td>
            </tr>
        """
    
    html_content += """
        </table>
    </div>
    
    <div class="section">
        <h2>Analysis Overview</h2>
        <p>This comprehensive codon usage analysis report summarizes the results of multiple codon bias metrics calculated from the input sequences.</p>
        
        <h3>Key Findings:</h3>
        <ul>
            <li><strong>CAI (Codon Adaptation Index)</strong>: Measures how closely gene sequences match the most abundant codons in the organism. Values range from 0 (no bias) to 1 (maximum bias).</li>
            <li><strong>ENC (Effective Number of Codons)</strong>: Quantifies codon usage bias. Values range from 20 (maximum bias) to 61 (no bias).</li>
            <li><strong>RSCU (Relative Synonymous Codon Usage)</strong>: Shows the relative usage of each codon compared to what would be expected under equal usage of synonymous codons. Values >1.0 indicate preferred codons.</li>
            <li><strong>tAI (tRNA Adaptation Index)</strong>: Estimates translational efficiency based on tRNA gene copy numbers.</li>
        </ul>
    </div>
    
    <div class="section">
        <h2>Visualizations</h2>
        <p>The following plots provide visual representations of the codon usage patterns:</p>
        
        <div class="plot-container">
            <h3>Codon Usage Heatmap</h3>
            <img src="../plots/codon_usage_heatmap.pdf" alt="Codon Usage Heatmap" class="plot" onerror="this.onerror=null; this.src='placeholder.png';">
        </div>
        
        <div class="plot-container">
            <h3>CAI Distribution</h3>
            <img src="../plots/cai_distribution.pdf" alt="CAI Distribution" class="plot" onerror="this.onerror=null; this.src='placeholder.png';">
        </div>
        
        <div class="plot-container">
            <h3>ENC Plot (GC3 vs ENC)</h3>
            <img src="../plots/enc_plot.pdf" alt="ENC Plot" class="plot" onerror="this.onerror=null; this.src='placeholder.png';">
        </div>
        
        <div class="plot-container">
            <h3>Neutrality Plot (GC12 vs GC3)</h3>
            <img src="../plots/neutrality_plot.pdf" alt="Neutrality Plot" class="plot" onerror="this.onerror=null; this.src='placeholder.png';">
        </div>
        
        <div class="plot-container">
            <h3>PR2 Plot</h3>
            <img src="../plots/pr2_plot.pdf" alt="PR2 Plot" class="plot" onerror="this.onerror=null; this.src='placeholder.png';">
        </div>
    </div>
    
    <div class="section">
        <h2>Interpretation</h2>
        <p>Codon usage bias reflects the evolutionary forces that have shaped the genome, including mutational biases, translational selection, and genetic drift. The relative contributions of these forces can be inferred from the patterns observed in the various metrics:</p>
        
        <ul>
            <li><strong>Strong codon bias</strong> (low ENC, high CAI) suggests translational selection for efficient protein synthesis.</li>
            <li><strong>Correlation between GC12 and GC3</strong> in neutrality plots suggests mutation pressure as the primary force.</li>
            <li><strong>Deviations from PR2 parity</strong> indicate directional mutational pressure.</li>
        </ul>
    </div>
    
    <div class="section">
        <h2>Methods</h2>
        <p>All analyses were performed using the Codon Usage Analyzer pipeline. Standard formulas were used for calculating each metric:</p>
        <ul>
            <li>CAI: Sharp & Li (1987)</li>
            <li>ENC: Wright (1990)</li>
            <li>RSCU: Sharp et al. (1986)</li>
            <li>tAI: dos Reis et al. (2004)</li>
        </ul>
    </div>
    
    <div class="section">
        <h2>References</h2>
        <ul>
            <li>Sharp PM, Li WH. The codon Adaptation Index-a measure of directional synonymous codon usage bias, and its potential applications. Nucleic Acids Res. 1987;15(3):1281-1295.</li>
            <li>Wright F. The 'effective number of codons' used in a gene. Gene. 1990;87(1):23-29.</li>
            <li>Sharp PM, Tuohy TM, Mosurski KR. Codon usage in yeast: cluster analysis clearly differentiates highly and lowly expressed genes. Nucleic Acids Res. 1986;14(13):5125-5143.</li>
            <li>dos Reis M, Savva R, Wernisch L. Solving the riddle of codon usage preferences: a test for translational selection. Nucleic Acids Res. 2004;32(17):5036-5044.</li>
        </ul>
    </div>
</body>
</html>
    """
    
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(html_content)

def main():
    parser = argparse.ArgumentParser(description='Generate comprehensive analysis report')
    parser.add_argument('output_dir', help='Output directory containing analysis results')
    parser.add_argument('config_file', help='Configuration file used for analysis')
    parser.add_argument('output_report', help='Output HTML report file')
    
    args = parser.parse_args()
    
    # Load configuration
    with open(args.config_file, 'r') as f:
        import yaml
        config = yaml.safe_load(f)
    
    config_info = {
        'input_type': config.get('input', {}).get('type', 'N/A'),
        'genetic_code': config.get('genetic_code', {}).get('table', 'N/A'),
        'calculate_rscu': config.get('analysis', {}).get('calculate_rscu', False),
        'calculate_cai': config.get('analysis', {}).get('calculate_cai', False),
        'calculate_enc': config.get('analysis', {}).get('calculate_enc', False),
        'calculate_tai': config.get('analysis', {}).get('calculate_tai', False),
    }
    
    # Load analysis results
    print("Loading analysis results...")
    results = load_analysis_results(args.output_dir)
    
    # Generate summary statistics
    print("Generating summary statistics...")
    summary = generate_statistics_summary(results)
    
    # Generate HTML report
    print(f"Generating comprehensive report: {args.output_report}")
    generate_html_report(results, summary, args.output_report, config_info)
    
    print("Comprehensive analysis report generated successfully!")

if __name__ == "__main__":
    main()