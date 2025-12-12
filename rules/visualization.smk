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

# ------------------------------------------------------------------------------
# Generate tAI distribution plot
# ------------------------------------------------------------------------------

rule generate_tai_distribution:
    """Generate distribution plot of tAI values"""
    input:
        tai_table = os.path.join(config["output"]["tables_dir"], "tai_values.tsv")
    output:
        tai_plot = os.path.join(config["output"]["plots_dir"], "tai_distribution.pdf"),
        tai_histogram = os.path.join(config["output"]["plots_dir"], "tai_histogram.png")
    log:
        os.path.join(config["output"]["base_dir"], "logs", "generate_tai_distribution.log")
    conda:
        "../envs/r_analysis.yaml"
    script:
        "../scripts/plot_tai_distribution.R"

# ------------------------------------------------------------------------------
# Generate CBI distribution plot
# ------------------------------------------------------------------------------

rule generate_cbi_distribution:
    """Generate distribution plot of CBI values"""
    input:
        cbi_table = os.path.join(config["output"]["tables_dir"], "cbi_values.tsv")
    output:
        cbi_plot = os.path.join(config["output"]["plots_dir"], "cbi_distribution.pdf"),
        cbi_histogram = os.path.join(config["output"]["plots_dir"], "cbi_histogram.png")
    log:
        os.path.join(config["output"]["base_dir"], "logs", "generate_cbi_distribution.log")
    conda:
        "../envs/r_analysis.yaml"
    script:
        "../scripts/plot_cbi_distribution.R"

# ------------------------------------------------------------------------------
# Generate FOP distribution plot
# ------------------------------------------------------------------------------

rule generate_fop_distribution:
    """Generate distribution plot of FOP values"""
    input:
        fop_table = os.path.join(config["output"]["tables_dir"], "fop_values.tsv")
    output:
        fop_plot = os.path.join(config["output"]["plots_dir"], "fop_distribution.pdf"),
        fop_histogram = os.path.join(config["output"]["plots_dir"], "fop_histogram.png")
    log:
        os.path.join(config["output"]["base_dir"], "logs", "generate_fop_distribution.log")
    conda:
        "../envs/r_analysis.yaml"
    script:
        "../scripts/plot_fop_distribution.R"

# ------------------------------------------------------------------------------
# Generate correlation plots between codon bias metrics
# ------------------------------------------------------------------------------

rule generate_bias_correlation:
    """Generate correlation plots between different codon bias metrics"""
    input:
        cai_table = os.path.join(config["output"]["tables_dir"], "cai_values.tsv"),
        enc_table = os.path.join(config["output"]["tables_dir"], "enc_values.tsv"),
        tai_table = os.path.join(config["output"]["tables_dir"], "tai_values.tsv"),
        cbi_table = os.path.join(config["output"]["tables_dir"], "cbi_values.tsv"),
        fop_table = os.path.join(config["output"]["tables_dir"], "fop_values.tsv")
    output:
        correlation_plot = os.path.join(config["output"]["plots_dir"], "bias_metrics_correlation.pdf"),
        correlation_matrix = os.path.join(config["output"]["plots_dir"], "bias_metrics_correlation_matrix.tsv")
    log:
        os.path.join(config["output"]["base_dir"], "logs", "generate_bias_correlation.log")
    conda:
        "../envs/r_analysis.yaml"
    script:
        "../scripts/plot_bias_correlation.R"

# ------------------------------------------------------------------------------
# Generate neutrality analysis plot (GC12 vs GC3)
# ------------------------------------------------------------------------------

rule generate_neutrality_plot:
    """Generate neutrality analysis plot (GC12 vs GC3)"""
    input:
        codon_table = os.path.join(config["output"]["tables_dir"], "codon_usage_table.tsv")
    output:
        neutrality_plot = os.path.join(config["output"]["plots_dir"], "neutrality_plot.pdf")
    log:
        os.path.join(config["output"]["base_dir"], "logs", "generate_neutrality_plot.log")
    conda:
        "../envs/r_analysis.yaml"
    script:
        "../scripts/plot_neutrality_analysis.R"

# ------------------------------------------------------------------------------
# Generate PR2 (Parity Rule 2) analysis plot
# ------------------------------------------------------------------------------

rule generate_pr2_plot:
    """Generate PR2 (Parity Rule 2) analysis plot"""
    input:
        codon_table = os.path.join(config["output"]["tables_dir"], "codon_usage_table.tsv")
    output:
        pr2_plot = os.path.join(config["output"]["plots_dir"], "pr2_plot.pdf")
    log:
        os.path.join(config["output"]["base_dir"], "logs", "generate_pr2_plot.log")
    conda:
        "../envs/r_analysis.yaml"
    script:
        "../scripts/plot_pr2_analysis.R"

# ------------------------------------------------------------------------------
# Generate codon usage effectiveness analysis plot
# ------------------------------------------------------------------------------

rule generate_codon_effectiveness_plot:
    """Generate codon usage effectiveness analysis plot"""
    input:
        codon_table = os.path.join(config["output"]["tables_dir"], "codon_usage_table.tsv"),
        cai_table = os.path.join(config["output"]["tables_dir"], "cai_values.tsv")
    output:
        effectiveness_plot = os.path.join(config["output"]["plots_dir"], "codon_usage_effectiveness.pdf")
    log:
        os.path.join(config["output"]["base_dir"], "logs", "generate_codon_effectiveness_plot.log")
    conda:
        "../envs/r_analysis.yaml"
    script:
        "../scripts/plot_codon_usage_effectiveness.R"