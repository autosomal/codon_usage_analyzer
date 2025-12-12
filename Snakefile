"""
Codon Usage Analyzer Workflow
Snakemake workflow for comprehensive codon usage analysis
"""

import os
import yaml
from datetime import datetime

# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------

# Load configuration file
configfile: "config.yaml"

# Set global variables
now = datetime.now()
config["timestamp"] = now.strftime("%Y%m%d_%H%M%S")
config["date"] = now.strftime("%Y-%m-%d")

# Create output directories
OUTPUT_DIRS = [
    config["output"]["base_dir"],
    config["output"]["plots_dir"],
    config["output"]["tables_dir"],
    config["output"]["reports_dir"],
    os.path.join(config["output"]["base_dir"], "preprocessing"),
    os.path.join(config["output"]["base_dir"], "quality_control"),
    os.path.join(config["output"]["base_dir"], "intermediate"),
]

# ------------------------------------------------------------------------------
# Include rules from separate files
# ------------------------------------------------------------------------------

include: "rules/preprocessing.smk"
include: "rules/codon_calculation.smk"
include: "rules/bias_analysis.smk"
include: "rules/tai_analysis.smk"
include: "rules/visualization.smk"
include: "rules/reporting.smk"

# ------------------------------------------------------------------------------
# Utility functions
# ------------------------------------------------------------------------------

def get_input_files():
    """Get input files based on input type"""
    input_type = config["input"]["type"]
    
    if input_type == "fastq":
        return expand(config["input"]["fastq_files"], 
                     sample=config["input"].get("samples", []),
                     read=config["input"].get("reads", []))
    elif input_type == "fasta":
        return config["input"]["files"]
    elif input_type == "cds":
        return config["input"]["cds_files"]
    elif input_type == "codonw":
        return config["input"]["codonw_files"]
    else:
        raise ValueError(f"Unsupported input type: {input_type}")

# ------------------------------------------------------------------------------
# Create directories
# ------------------------------------------------------------------------------

rule all:
    """Build all final targets"""
    input:
        # Main output files
        os.path.join(config["output"]["tables_dir"], "codon_usage_table.tsv"),
        os.path.join(config["output"]["tables_dir"], "rscu_values.tsv"),
        os.path.join(config["output"]["tables_dir"], "cai_values.tsv"),
        os.path.join(config["output"]["tables_dir"], "enc_values.tsv"),
        os.path.join(config["output"]["tables_dir"], "preferred_codons.tsv"),
        
        # Advanced analysis output files
        os.path.join(config["output"]["tables_dir"], "tai_values.tsv") if config["analysis"]["calculate_tai"] else None,
        os.path.join(config["output"]["tables_dir"], "cbi_values.tsv") if config["analysis"]["calculate_cbi"] else None,
        os.path.join(config["output"]["tables_dir"], "fop_values.tsv") if config["analysis"]["calculate_fop"] else None,
        os.path.join(config["output"]["tables_dir"], "codon_pair_usage.tsv") if config["analysis"]["codon_pair_analysis"] else None,
        os.path.join(config["output"]["tables_dir"], "codon_expression_correlation.tsv") if config["analysis"]["correlation_analysis"] else None,
        os.path.join(config["output"]["tables_dir"], "cross_species_comparison.tsv") if config["analysis"].get("comparison_species") else None,
        
        # Visualization outputs
        os.path.join(config["output"]["plots_dir"], "codon_usage_heatmap.pdf"),
        os.path.join(config["output"]["plots_dir"], "rscu_correlation.pdf"),
        os.path.join(config["output"]["plots_dir"], "cai_distribution.pdf"),
        os.path.join(config["output"]["plots_dir"], "neutrality_plot.pdf"),
        os.path.join(config["output"]["plots_dir"], "pr2_plot.pdf"),
        os.path.join(config["output"]["plots_dir"], "codon_usage_effectiveness.pdf"),
        
        # Advanced visualization outputs
        os.path.join(config["output"]["plots_dir"], "tai_distribution.pdf") if config["analysis"]["calculate_tai"] else None,
        os.path.join(config["output"]["plots_dir"], "cbi_distribution.pdf") if config["analysis"]["calculate_cbi"] else None,
        os.path.join(config["output"]["plots_dir"], "fop_distribution.pdf") if config["analysis"]["calculate_fop"] else None,
        os.path.join(config["output"]["plots_dir"], "bias_metrics_correlation.pdf") if config["analysis"]["calculate_tai"] or config["analysis"]["calculate_cbi"] or config["analysis"]["calculate_fop"] else None,
        
        # Final report
        os.path.join(config["output"]["reports_dir"], "final_report.html")

rule create_directories:
    """Create all necessary directories"""
    output:
        directories = [directory for directory in OUTPUT_DIRS]
    shell:
        """
        for dir in {output.directories}; do
            mkdir -p "$dir"
        done
        """

# ------------------------------------------------------------------------------
# Workflow visualization
# ------------------------------------------------------------------------------

rule plot_workflow:
    """Plot the workflow DAG"""
    output:
        "workflow_dag.pdf"
    shell:
        """
        snakemake --dag | dot -Tpdf > {output}
        """

rule plot_rulegraph:
    """Plot the rule graph"""
    output:
        "rule_graph.pdf"
    shell:
        """
        snakemake --rulegraph | dot -Tpdf > {output}
        """