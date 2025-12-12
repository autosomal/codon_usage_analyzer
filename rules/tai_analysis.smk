"""
tRNA Adaptation Index (tAI) analysis rules
"""

# ------------------------------------------------------------------------------
# Calculate tAI (tRNA Adaptation Index)
# ------------------------------------------------------------------------------

rule calculate_tai:
    """Calculate tRNA Adaptation Index (tAI) for genes"""
    input:
        codon_table = os.path.join(config["output"]["tables_dir"], "codon_usage_table.tsv"),
        trna_data = config["analysis"]["trna_data"]
    output:
        tai_table = os.path.join(config["output"]["tables_dir"], "tai_values.tsv")
    params:
        genetic_code = config["genetic_code"]
    log:
        os.path.join(config["output"]["base_dir"], "logs", "calculate_tai.log")
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/calculate_tai.py"

# ------------------------------------------------------------------------------
# Calculate CBI (Codon Bias Index)
# ------------------------------------------------------------------------------

rule calculate_cbi:
    """Calculate Codon Bias Index (CBI)"""
    input:
        codon_table = os.path.join(config["output"]["tables_dir"], "codon_usage_table.tsv")
    output:
        cbi_table = os.path.join(config["output"]["tables_dir"], "cbi_values.tsv")
    params:
        genetic_code = config["genetic_code"]
    log:
        os.path.join(config["output"]["base_dir"], "logs", "calculate_cbi.log")
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/calculate_cbi.py"

# ------------------------------------------------------------------------------
# Calculate FOP (Frequency of Optimal Codons)
# ------------------------------------------------------------------------------

rule calculate_fop:
    """Calculate Frequency of Optimal Codons (FOP)"""
    input:
        codon_table = os.path.join(config["output"]["tables_dir"], "codon_usage_table.tsv"),
        preferred_codons = os.path.join(config["output"]["tables_dir"], "preferred_codons.tsv")
    output:
        fop_table = os.path.join(config["output"]["tables_dir"], "fop_values.tsv")
    params:
        genetic_code = config["genetic_code"]
    log:
        os.path.join(config["output"]["base_dir"], "logs", "calculate_fop.log")
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/calculate_fop.py"

# ------------------------------------------------------------------------------
# Codon Pair Analysis
# ------------------------------------------------------------------------------

rule codon_pair_analysis:
    """Analyze codon pair usage bias"""
    input:
        fasta = config["input"]["files"]
    output:
        codon_pairs = os.path.join(config["output"]["tables_dir"], "codon_pair_usage.tsv"),
        dipeptide_bias = os.path.join(config["output"]["tables_dir"], "dipeptide_bias.tsv")
    params:
        genetic_code = config["genetic_code"],
        calculate_cpbia = config["analysis"]["codon_pair_analysis"]
    log:
        os.path.join(config["output"]["base_dir"], "logs", "codon_pair_analysis.log")
    conda:
        "../envs/bioinformatics.yaml"
    script:
        "../scripts/codon_pair_analysis.py"