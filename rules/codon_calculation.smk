"""
Codon usage calculation rules
"""

# ------------------------------------------------------------------------------
# Calculate codon usage from FASTA file
# ------------------------------------------------------------------------------

rule calculate_codon_usage_fasta:
    """Calculate codon usage from FASTA file"""
    input:
        fasta = config["input"]["files"]
    output:
        codon_table = os.path.join(config["output"]["tables_dir"], "codon_usage_table.tsv"),
        amino_acid_table = os.path.join(config["output"]["tables_dir"], "amino_acid_usage.tsv")
    params:
        genetic_code = config["genetic_code"],
        calculate_rscu = config["analysis"]["calculate_rscu"],
        calculate_cai = config["analysis"]["calculate_cai"],
        calculate_enc = config["analysis"]["calculate_enc"]
    log:
        os.path.join(config["output"]["base_dir"], "logs", "calculate_codon_usage.log")
    conda:
        "../envs/base.yaml"
    threads:
        config["resources"]["threads"]
    script:
        "../scripts/calculate_codon_usage.py"

# ------------------------------------------------------------------------------
# Calculate codon usage from parsed CodonW file
# ------------------------------------------------------------------------------

rule calculate_codon_usage_codonw:
    """Calculate codon usage from parsed CodonW file"""
    input:
        parsed_file = os.path.join(config["output"]["base_dir"], "preprocessing", "codonw_parsed.tsv")
    output:
        codon_table = os.path.join(config["output"]["tables_dir"], "codon_usage_table.tsv"),
        amino_acid_table = os.path.join(config["output"]["tables_dir"], "amino_acid_usage.tsv")
    params:
        genetic_code = config["genetic_code"],
        calculate_rscu = config["analysis"]["calculate_rscu"],
        calculate_cai = config["analysis"]["calculate_cai"],
        calculate_enc = config["analysis"]["calculate_enc"]
    log:
        os.path.join(config["output"]["base_dir"], "logs", "calculate_codon_usage_codonw.log")
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/calculate_codon_usage_from_table.py"

# ------------------------------------------------------------------------------
# Calculate RSCU (Relative Synonymous Codon Usage)
# ------------------------------------------------------------------------------

rule calculate_rscu:
    """Calculate Relative Synonymous Codon Usage (RSCU)"""
    input:
        codon_table = os.path.join(config["output"]["tables_dir"], "codon_usage_table.tsv")
    output:
        rscu_table = os.path.join(config["output"]["tables_dir"], "rscu_values.tsv")
    params:
        genetic_code = config["genetic_code"]
    log:
        os.path.join(config["output"]["base_dir"], "logs", "calculate_rscu.log")
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/calculate_rscu.py"

# ------------------------------------------------------------------------------
# Calculate CAI (Codon Adaptation Index)
# ------------------------------------------------------------------------------

rule calculate_cai:
    """Calculate Codon Adaptation Index (CAI)"""
    input:
        codon_table = os.path.join(config["output"]["tables_dir"], "codon_usage_table.tsv"),
        rscu_table = os.path.join(config["output"]["tables_dir"], "rscu_values.tsv")
    output:
        cai_table = os.path.join(config["output"]["tables_dir"], "cai_values.tsv")
    params:
        genetic_code = config["genetic_code"]
    log:
        os.path.join(config["output"]["base_dir"], "logs", "calculate_cai.log")
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/calculate_cai.py"

# ------------------------------------------------------------------------------
# Calculate ENC (Effective Number of Codons)
# ------------------------------------------------------------------------------

rule calculate_enc:
    """Calculate Effective Number of Codons (ENC)"""
    input:
        codon_table = os.path.join(config["output"]["tables_dir"], "codon_usage_table.tsv")
    output:
        enc_table = os.path.join(config["output"]["tables_dir"], "enc_values.tsv")
    params:
        genetic_code = config["genetic_code"]
    log:
        os.path.join(config["output"]["base_dir"], "logs", "calculate_enc.log")
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/calculate_enc.py"