"""
Reporting rules for codon usage analysis
"""

# ------------------------------------------------------------------------------
# Generate MultiQC report
# ------------------------------------------------------------------------------

rule multiqc_report:
    """Generate MultiQC report for quality control"""
    input:
        fastqc = expand(os.path.join(config["output"]["base_dir"], "quality_control", "fastqc", "*.html"), glob=True),
        star_logs = expand(os.path.join(config["output"]["base_dir"], "preprocessing", "*_star_log.txt"), glob=True),
        stringtie_logs = expand(os.path.join(config["output"]["base_dir"], "preprocessing", "stringtie_*.log"), glob=True)
    output:
        directory(os.path.join(config["output"]["reports_dir"], "multiqc"))
    log:
        os.path.join(config["output"]["base_dir"], "logs", "multiqc.log")
    conda:
        "../envs/multiqc.yaml"
    shell:
        """
        multiqc {config["output"]["base_dir"]} \\
            -o {output} \\
            --title "Codon Usage Analyzer QC Report" 2>&1 | tee {log}
        """

# ------------------------------------------------------------------------------
# Generate codon usage summary report
# ------------------------------------------------------------------------------

rule generate_summary_report:
    """Generate summary report of codon usage analysis"""
    input:
        codon_table = os.path.join(config["output"]["tables_dir"], "codon_usage_table.tsv"),
        rscu_table = os.path.join(config["output"]["tables_dir"], "rscu_values.tsv"),
        cai_table = os.path.join(config["output"]["tables_dir"], "cai_values.tsv"),
        enc_table = os.path.join(config["output"]["tables_dir"], "enc_values.tsv"),
        preferred_codons = os.path.join(config["output"]["tables_dir"], "preferred_codons.tsv"),
        plots = expand(os.path.join(config["output"]["plots_dir"], "*.pdf"), glob=True),
        multiqc_report = os.path.join(config["output"]["reports_dir"], "multiqc", "multiqc_report.html")
    output:
        summary_report = os.path.join(config["output"]["reports_dir"], "summary_report.tsv"),
        summary_stats = os.path.join(config["output"]["reports_dir"], "summary_statistics.txt")
    params:
        config = config
    log:
        os.path.join(config["output"]["base_dir"], "logs", "generate_summary_report.log")
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/generate_summary_report.py"

# ------------------------------------------------------------------------------
# Generate final HTML report
# ------------------------------------------------------------------------------

rule generate_final_report:
    """Generate final HTML report with all analysis results"""
    input:
        summary_report = os.path.join(config["output"]["reports_dir"], "summary_report.tsv"),
        summary_stats = os.path.join(config["output"]["reports_dir"], "summary_statistics.txt"),
        codon_table = os.path.join(config["output"]["tables_dir"], "codon_usage_table.tsv"),
        cai_table = os.path.join(config["output"]["tables_dir"], "cai_values.tsv"),
        enc_table = os.path.join(config["output"]["tables_dir"], "enc_values.tsv"),
        rscu_table = os.path.join(config["output"]["tables_dir"], "rscu_values.tsv"),
        tai_table = os.path.join(config["output"]["tables_dir"], "tai_values.tsv") if config["analysis"]["calculate_tai"] else os.path.join(config["output"]["tables_dir"], "cai_values.tsv"),  # fallback
        plots = expand(os.path.join(config["output"]["plots_dir"], "*.pdf"), glob=True),
        tables = expand(os.path.join(config["output"]["tables_dir"], "*.tsv"), glob=True),
        multiqc_report = os.path.join(config["output"]["reports_dir"], "multiqc", "multiqc_report.html")
    output:
        final_report = os.path.join(config["output"]["reports_dir"], "final_report.html")
    params:
        config_file = "config.yaml",
        output_dir = config["output"]["tables_dir"]
    log:
        os.path.join(config["output"]["base_dir"], "logs", "generate_final_report.log")
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/generate_comprehensive_report.py"