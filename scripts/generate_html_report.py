#!/usr/bin/env python3
"""
Generate final HTML report for codon usage analysis
"""

import os
import sys
import json
import pandas as pd
from jinja2 import Environment, FileSystemLoader

def main():
    # Get input/output from Snakemake
    summary_report = snakemake.input.summary_report
    summary_stats = snakemake.input.summary_stats
    plots = snakemake.input.plots
    tables = snakemake.input.tables
    output_file = snakemake.output[0]
    
    # Get parameters
    config = snakemake.params.config
    timestamp = snakemake.params.timestamp
    
    # Read summary statistics
    stats_data = {}
    with open(summary_stats, 'r') as f:
        for line in f:
            line = line.strip()
            if ':' in line:
                key, value = line.split(':', 1)
                stats_data[key.strip()] = value.strip()
    
    # Read summary report
    summary_df = pd.read_csv(summary_report, sep='\t')
    
    # Organize plots by type
    plot_types = {
        'heatmap': [],
        'correlation': [],
        'distribution': [],
        'pca': [],
        'other': []
    }
    
    for plot in plots:
        plot_name = os.path.basename(plot)
        if 'heatmap' in plot_name:
            plot_types['heatmap'].append(plot_name)
        elif 'correlation' in plot_name:
            plot_types['correlation'].append(plot_name)
        elif 'distribution' in plot_name or 'violin' in plot_name:
            plot_types['distribution'].append(plot_name)
        elif 'pca' in plot_name:
            plot_types['pca'].append(plot_name)
        else:
            plot_types['other'].append(plot_name)
    
    # Organize tables by type
    table_types = {
        'codon_usage': [],
        'rscu': [],
        'cai': [],
        'enc': [],
        'preferred': [],
        'avoided': [],
        'other': []
    }
    
    for table in tables:
        table_name = os.path.basename(table)
        if 'codon_usage' in table_name:
            table_types['codon_usage'].append(table_name)
        elif 'rscu' in table_name:
            table_types['rscu'].append(table_name)
        elif 'cai' in table_name:
            table_types['cai'].append(table_name)
        elif 'enc' in table_name:
            table_types['enc'].append(table_name)
        elif 'preferred' in table_name:
            table_types['preferred'].append(table_name)
        elif 'avoided' in table_name:
            table_types['avoided'].append(table_name)
        else:
            table_types['other'].append(table_name)
    
    # Create HTML template
    template = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Codon Usage Analysis Report</title>
    <style>
        body { font-family: Arial, sans-serif; line-height: 1.6; margin: 0; padding: 20px; }
        .container { max-width: 1200px; margin: 0 auto; }
        h1, h2, h3 { color: #2c3e50; }
        .header { background-color: #f8f9fa; padding: 20px; border-radius: 8px; margin-bottom: 30px; }
        .summary { background-color: #e8f4f8; padding: 20px; border-radius: 8px; margin-bottom: 30px; }
        .section { margin-bottom: 40px; }
        .plot-container { margin-bottom: 30px; }
        .plot-title { font-weight: bold; margin-bottom: 10px; }
        .table-container { overflow-x: auto; margin-bottom: 30px; }
        table { border-collapse: collapse; width: 100%; }
        th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
        th { background-color: #f2f2f2; }
        tr:nth-child(even) { background-color: #f9f9f9; }
        .footer { margin-top: 50px; padding-top: 20px; border-top: 1px solid #ddd; text-align: center; color: #666; }
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>Codon Usage Analysis Report</h1>
            <p>Analysis Date: {{ timestamp }}</p>
            <p>Organism: {{ config.get('organism', 'Unknown') }}</p>
            <p>Input Type: {{ config.get('input', {}).get('type', 'Unknown') }}</p>
        </div>

        <div class="summary">
            <h2>Summary Statistics</h2>
            <ul>
                {% for key, value in stats_data.items() %}
                <li><strong>{{ key }}:</strong> {{ value }}</li>
                {% endfor %}
            </ul>
        </div>

        <div class="section">
            <h2>Codon Usage Overview</h2>
            <p>This section provides an overview of codon usage patterns in the analyzed dataset.</p>
            
            <div class="plot-container">
                <h3>Codon Usage Heatmap</h3>
                {% for plot in plot_types.heatmap %}
                <img src="{{ plot }}" alt="Codon Usage Heatmap" style="max-width: 100%; height: auto;">
                {% endfor %}
            </div>
        </div>

        <div class="section">
            <h2>Codon Bias Analysis</h2>
            <p>This section presents the results of codon bias analysis, including RSCU values and preferred/avoided codons.</p>
            
            <div class="table-container">
                <h3>RSCU Values</h3>
                <table>
                    <tr>
                        <th>Codon</th>
                        <th>Amino Acid</th>
                        <th>RSCU</th>
                        <th>Count</th>
                        <th>Frequency</th>
                    </tr>
                    {% for _, row in summary_df.iterrows() %}
                    <tr>
                        <td>{{ row.codon }}</td>
                        <td>{{ row.amino_acid }}</td>
                        <td>{{ row.rscu|round(2) }}</td>
                        <td>{{ row.count }}</td>
                        <td>{{ row.frequency|round(4) }}</td>
                    </tr>
                    {% endfor %}
                </table>
            </div>
        </div>

        <div class="section">
            <h2>Visualizations</h2>
            
            <div class="plot-container">
                <h3>Codon Usage Correlation</h3>
                {% for plot in plot_types.correlation %}
                <img src="{{ plot }}" alt="Codon Usage Correlation" style="max-width: 100%; height: auto;">
                {% endfor %}
            </div>
            
            <div class="plot-container">
                <h3>CAI Distribution</h3>
                {% for plot in plot_types.distribution %}
                <img src="{{ plot }}" alt="CAI Distribution" style="max-width: 100%; height: auto;">
                {% endfor %}
            </div>
            
            <div class="plot-container">
                <h3>PCA Analysis</h3>
                {% for plot in plot_types.pca %}
                <img src="{{ plot }}" alt="PCA Analysis" style="max-width: 100%; height: auto;">
                {% endfor %}
            </div>
        </div>

        <div class="section">
            <h2>Data Tables</h2>
            <p>All analysis results are available in the following tables:</p>
            <ul>
                {% for table in table_types.codon_usage %}
                <li><a href="{{ table }}">{{ table }}</a> - Codon usage statistics</li>
                {% endfor %}
                {% for table in table_types.rscu %}
                <li><a href="{{ table }}">{{ table }}</a> - RSCU values</li>
                {% endfor %}
                {% for table in table_types.cai %}
                <li><a href="{{ table }}">{{ table }}</a> - CAI values</li>
                {% endfor %}
                {% for table in table_types.enc %}
                <li><a href="{{ table }}">{{ table }}</a> - ENC values</li>
                {% endfor %}
                {% for table in table_types.preferred %}
                <li><a href="{{ table }}">{{ table }}</a> - Preferred codons</li>
                {% endfor %}
                {% for table in table_types.avoided %}
                <li><a href="{{ table }}">{{ table }}</a> - Avoided codons</li>
                {% endfor %}
            </ul>
        </div>

        <div class="footer">
            <p>Generated by Codon Usage Analyzer workflow</p>
            <p>Report generated on: {{ timestamp }}</p>
        </div>
    </div>
</body>
</html>
    """
    
    # Render template
    env = Environment()
    html_content = env.from_string(template).render(
        timestamp=timestamp,
        config=config,
        stats_data=stats_data,
        summary_df=summary_df,
        plot_types=plot_types,
        table_types=table_types
    )
    
    # Save HTML report
    with open(output_file, 'w') as f:
        f.write(html_content)
    
    print(f"Final HTML report generated successfully")
    print(f"Report saved to: {output_file}")
    print(f"Total plots included: {len(plots)}")
    print(f"Total tables included: {len(tables)}")

if __name__ == "__main__":
    main()