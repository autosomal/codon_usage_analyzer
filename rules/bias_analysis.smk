"""
Codon bias analysis rules
"""

# ------------------------------------------------------------------------------
# Identify preferred and avoided codons
# ------------------------------------------------------------------------------

rule identify_preferred_codons:
    """Identify preferred and avoided codons based on RSCU values"""
    input:
        rscu_table = os.path.join(config["output"]["tables_dir"], "rscu_values.tsv")
    output:
        preferred_codons = os.path.join(config["output"]["tables_dir"], "preferred_codons.tsv"),
        avoided_codons = os.path.join(config["output"]["tables_dir"], "avoided_codons.tsv")
    params:
        preferred_threshold = 1.2,
        avoided_threshold = 0.8
    log:
        os.path.join(config["output"]["base_dir"], "logs", "identify_preferred_codons.log")
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/identify_preferred_codons.py"

# ------------------------------------------------------------------------------
# Codon pair bias analysis
# ------------------------------------------------------------------------------

rule codon_pair_analysis:
    """Analyze codon pair bias"""
    input:
        fasta = config["input"]["files"]
    output:
        pair_table = os.path.join(config["output"]["tables_dir"], "codon_pair_usage.tsv"),
        pair_bias = os.path.join(config["output"]["tables_dir"], "codon_pair_bias.tsv")
    params:
        genetic_code = config["genetic_code"],
        min_count = 10
    log:
        os.path.join(config["output"]["base_dir"], "logs", "codon_pair_analysis.log")
    conda:
        "../envs/base.yaml"
    threads:
        config["resources"]["threads"]
    script:
        "../scripts/codon_pair_analysis.py"

# ------------------------------------------------------------------------------
# Correlation analysis between codon usage and gene expression
# ------------------------------------------------------------------------------

rule correlation_analysis:
    """Analyze correlation between codon usage and gene expression"""
    input:
        codon_table = os.path.join(config["output"]["tables_dir"], "codon_usage_table.tsv"),
        expression_data = lambda wildcards: config["input"].get("expression_data", "data/expression_data.tsv") if config["input"].get("expression_data") else "data/expression_data.tsv"
    output:
        correlation_table = os.path.join(config["output"]["tables_dir"], "codon_expression_correlation.tsv")
    params:
        genetic_code = config["genetic_code"]
    log:
        os.path.join(config["output"]["base_dir"], "logs", "correlation_analysis.log")
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/correlation_analysis.py"

# ------------------------------------------------------------------------------
# Cross-species codon usage comparison
# ------------------------------------------------------------------------------

rule cross_species_comparison:
    """Compare codon usage across multiple species"""
    input:
        # This would typically take multiple codon usage tables
        # For now, we'll use a placeholder - in a real workflow this would be expand() with species
        codon_tables = lambda wildcards: expand(os.path.join(config["output"]["tables_dir"], "{species}_codon_usage.tsv"),
                            species=config["analysis"].get("comparison_species", [])) if config["analysis"].get("comparison_species") else [os.path.join(config["output"]["tables_dir"], "codon_usage_table.tsv")]
    output:
        comparison_table = os.path.join(config["output"]["tables_dir"], "cross_species_comparison.tsv"),
        distance_matrix = os.path.join(config["output"]["tables_dir"], "codon_usage_distance.tsv")
    log:
        os.path.join(config["output"]["base_dir"], "logs", "cross_species_comparison.log")
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/cross_species_comparison.py"