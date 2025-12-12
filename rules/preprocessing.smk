"""
Preprocessing rules for codon usage analysis workflow
"""

# ------------------------------------------------------------------------------
# Quality control for FASTQ files
# ------------------------------------------------------------------------------

rule fastqc:
    """Run FastQC on FASTQ files"""
    input:
        fastq_files = expand(config["input"]["fastq_files"], 
                           sample=config["input"].get("samples", []),
                           read=config["input"].get("reads", []))
    output:
        html = expand(os.path.join(config["output"]["base_dir"], "quality_control", "fastqc", "{sample}_R{read}_fastqc.html"),
                     sample=config["input"].get("samples", []),
                     read=config["input"].get("reads", [])),
        zip = expand(os.path.join(config["output"]["base_dir"], "quality_control", "fastqc", "{sample}_R{read}_fastqc.zip"),
                    sample=config["input"].get("samples", []),
                    read=config["input"].get("reads", []))
    log:
        os.path.join(config["output"]["base_dir"], "quality_control", "fastqc.log")
    conda:
        "../envs/bioinformatics.yaml"
    threads:
        config["resources"]["threads"]
    shell:
        """
        mkdir -p {output.html[0]%/*}
        fastqc {input.fastq_files} \\
            -o {output.html[0]%/*} \\
            -t {threads} 2>&1 | tee {log}
        """

# ------------------------------------------------------------------------------
# Alignment with STAR
# ------------------------------------------------------------------------------

rule star_alignment:
    """Align reads to reference genome using STAR"""
    input:
        genome_index = directory(os.path.join(config["output"]["base_dir"], "preprocessing", "star_index")),
        fastq_files = expand(config["input"]["fastq_files"], 
                           sample=wildcards.sample,
                           read=config["input"].get("reads", []))
    output:
        bam = os.path.join(config["output"]["base_dir"], "preprocessing", "{sample}_aligned.bam"),
        bai = os.path.join(config["output"]["base_dir"], "preprocessing", "{sample}_aligned.bam.bai"),
        log = os.path.join(config["output"]["base_dir"], "preprocessing", "{sample}_star_log.txt")
    log:
        os.path.join(config["output"]["base_dir"], "preprocessing", "star_{sample}.log")
    conda:
        "../envs/bioinformatics.yaml"
    threads:
        config["resources"]["threads"]
    resources:
        mem_mb = int(config["resources"]["memory"].rstrip("G")) * 1024
    shell:
        """
        STAR --genomeDir {input.genome_index} \\
            --readFilesIn {input.fastq_files} \\
            --readFilesCommand zcat \\
            --outFileNamePrefix {output.bam%.bam}_ \\
            --outSAMtype BAM SortedByCoordinate \\
            --runThreadN {threads} \\
            --genomeLoad NoSharedMemory 2>&1 | tee {log}
        
        samtools index {output.bam}
        """

# ------------------------------------------------------------------------------
# STAR genome index
# ------------------------------------------------------------------------------

rule star_index:
    """Build STAR genome index"""
    input:
        fasta = config["input"]["reference"],
        gtf = config["input"]["gtf"]
    output:
        directory(os.path.join(config["output"]["base_dir"], "preprocessing", "star_index"))
    log:
        os.path.join(config["output"]["base_dir"], "preprocessing", "star_index.log")
    conda:
        "../envs/bioinformatics.yaml"
    threads:
        config["resources"]["threads"]
    resources:
        mem_mb = int(config["resources"]["memory"].rstrip("G")) * 1024
    shell:
        """
        STAR --runMode genomeGenerate \\
            --genomeDir {output} \\
            --genomeFastaFiles {input.fasta} \\
            --sjdbGTFfile {input.gtf} \\
            --sjdbOverhang 100 \\
            --runThreadN {threads} 2>&1 | tee {log}
        """

# ------------------------------------------------------------------------------
# Transcript assembly with StringTie
# ------------------------------------------------------------------------------

rule stringtie_assembly:
    """Assemble transcripts using StringTie"""
    input:
        bam = os.path.join(config["output"]["base_dir"], "preprocessing", "{sample}_aligned.bam"),
        gtf = config["input"]["gtf"]
    output:
        gtf = os.path.join(config["output"]["base_dir"], "preprocessing", "{sample}_transcripts.gtf"),
        fasta = os.path.join(config["output"]["base_dir"], "preprocessing", "{sample}_transcripts.fa")
    log:
        os.path.join(config["output"]["base_dir"], "preprocessing", "stringtie_{sample}.log")
    conda:
        "../envs/bioinformatics.yaml"
    threads:
        config["resources"]["threads"]
    shell:
        """
        stringtie {input.bam} \\
            -G {input.gtf} \\
            -o {output.gtf} \\
            -A {output.gtf%.gtf}_gene_abundances.tsv \\
            -C {output.gtf%.gtf}_cov_refs.gtf \\
            -l {wildcards.sample} \\
            -p {threads} 2>&1 | tee {log}
        
        # Extract transcript sequences
        gffread {output.gtf} \\
            -g {config["input"]["reference"]} \\
            -w {output.fasta} 2>&1 | tee -a {log}
        """

# ------------------------------------------------------------------------------
# Extract CDS sequences from GTF
# ------------------------------------------------------------------------------

rule extract_cds:
    """Extract coding sequences from GTF annotation"""
    input:
        fasta = config["input"]["reference"],
        gtf = config["input"]["gtf"]
    output:
        cds_fasta = os.path.join(config["output"]["base_dir"], "preprocessing", "cds_sequences.fa"),
        cds_gtf = os.path.join(config["output"]["base_dir"], "preprocessing", "cds_annotation.gtf")
    log:
        os.path.join(config["output"]["base_dir"], "preprocessing", "extract_cds.log")
    conda:
        "../envs/bioinformatics.yaml"
    script:
        "../scripts/extract_cds.py"

# ------------------------------------------------------------------------------
# Parse CodonW output file
# ------------------------------------------------------------------------------

rule parse_codonw:
    """Parse CodonW output file into standardized format"""
    input:
        codonw_file = config["input"]["codonw_files"]
    output:
        parsed_table = os.path.join(config["output"]["base_dir"], "preprocessing", "codonw_parsed.tsv")
    log:
        os.path.join(config["output"]["base_dir"], "preprocessing", "parse_codonw.log")
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/parse_codonw.py"