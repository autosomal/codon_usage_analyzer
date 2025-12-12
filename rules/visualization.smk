"""
Visualization rules for codon usage analysis
"""

# ------------------------------------------------------------------------------
# Generate codon usage heatmap
# ------------------------------------------------------------------------------

rule generate_codon_heatmap:
    """Generate heatmap of codon usage patterns"""
    input:
        codon_table = os.path.join(config["output"]["tables_dir"], "codon_usage_table.tsv"),
        rscu_table = os.path.join(config["output"]["tables_dir"], "rscu_values.tsv")
    output:
        heatmap = os.path.join(config["output"]["plots_dir"], "codon_usage_heatmap.pdf"),
        heatmap_png = os.path.join(config["output"]["plots_dir"], "codon_usage_heatmap.png")
    params:
        genetic_code = config["genetic_code"],
        plot_format = config["analysis"].get("plot_formats", ["pdf", "png"])
    log:
        os.path.join(config["output"]["base_dir"], "logs", "generate_codon_heatmap.log")
    conda:
        "../envs/r_analysis.yaml"
    script:
        "../scripts/plot_codon_heatmap.R"

# ------------------------------------------------------------------------------
# Generate RSCU correlation plots
# ------------------------------------------------------------------------------

rule generate_rscu_correlation:
    """Generate correlation plots for RSCU values"""
    input:
        rscu_table = os.path.join(config["output"]["tables_dir"], "rscu_values.tsv")
    output:
        correlation_plot = os.path.join(config["output"]["plots_dir"], "rscu_correlation.pdf"),
        correlation_matrix = os.path.join(config["output"]["plots_dir"], "rscu_correlation_matrix.pdf")
    params:
        genetic_code = config["genetic_code"]
    log:
        os.path.join(config["output"]["base_dir"], "logs", "generate_rscu_correlation.log")
    conda:
        "../envs/r_analysis.yaml"
    script:
        "../scripts/plot_rscu_correlation.R"

# ------------------------------------------------------------------------------
# Generate CAI distribution plot
# ------------------------------------------------------------------------------

rule generate_cai_distribution:
    """Generate distribution plot of CAI values"""
    input:
        cai_table = os.path.join(config["output"]["tables_dir"], "cai_values.tsv")
    output:
        cai_plot = os.path.join(config["output"]["plots_dir"], "cai_distribution.pdf"),
        cai_violin = os.path.join(config["output"]["plots_dir"], "cai_violin_plot.pdf")
    log:
        os.path.join(config["output"]["base_dir"], "logs", "generate_cai_distribution.log")
    conda:
        "../envs/r_analysis.yaml"
    script:
        "../scripts/plot_cai_distribution.R"

# ------------------------------------------------------------------------------
# Generate ENC-plot (ENC vs GC3)
# ------------------------------------------------------------------------------

rule generate_enc_plot:
    """Generate ENC-plot (ENC vs GC3 content)"""
    input:
        enc_table = os.path.join(config["output"]["tables_dir"], "enc_values.tsv"),
        codon_table = os.path.join(config["output"]["tables_dir"], "codon_usage_table.tsv")
    output:
        enc_plot = os.path.join(config["output"]["plots_dir"], "enc_plot.pdf")
    log:
        os.path.join(config["output"]["base_dir"], "logs", "generate_enc_plot.log")
    conda:
        "../envs/r_analysis.yaml"
    script:
        "../scripts/plot_enc_plot.R"

# ------------------------------------------------------------------------------
# Generate PCA plot of codon usage
# ------------------------------------------------------------------------------

rule generate_pca_plot:
    """Generate PCA plot of codon usage patterns"""
    input:
        codon_table = os.path.join(config["output"]["tables_dir"], "codon_usage_table.tsv")
    output:
        pca_plot = os.path.join(config["output"]["plots_dir"], "codon_usage_pca.pdf")
    params:
        group_by = config["analysis"].get("pca_group_by", "none")
    log:
        os.path.join(config["output"]["base_dir"], "logs", "generate_pca_plot.log")
    conda:
        "../envs/r_analysis.yaml"
    script:
        "../scripts/plot_codon_pca.R"