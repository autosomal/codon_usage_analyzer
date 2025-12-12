# Codon Usage Analyzer

A comprehensive Snakemake workflow for analyzing codon usage patterns from genomic and transcriptomic data.

## Overview

This project provides a scalable, reproducible workflow for analyzing codon usage bias and patterns. The workflow can process:
- Raw sequencing data (FASTQ)
- Assembled transcriptomes (FASTA)
- CDS sequences (FASTA)
- Pre-computed codon usage tables (from CodonW or other tools)

The workflow integrates multiple analysis steps into a unified pipeline with standardized output formats.

## Features

### **Multi-Input Support**
- **Raw sequencing data**: FASTQ files (paired-end or single-end)
- **Assembled transcripts**: FASTA files
- **CDS sequences**: FASTA files with coding sequences
- **Pre-computed tables**: CodonW output files or custom tables

### **Analysis Modules**
- **Codon Usage Calculation**: Compute RSCU, CAI, ENC, and other metrics
- **Codon Bias Analysis**: Identify preferred and avoided codons
- **tRNA Adaptation Index**: Calculate tAI for translation efficiency
- **Correlation Analysis**: Relate codon usage to gene expression
- **Visualization**: Generate publication-ready plots and heatmaps

### **Cross-Species Support**
- Works with any organism with genomic or transcriptomic data
- Built-in genetic code tables for different organisms
- Customizable codon tables for non-standard genetic codes

## Quick Start

### **1. Clone the repository**
```bash
git clone https://github.com/autosomal/codon_usage_analyzer.git
cd codon_usage_analyzer
```

### **2. Create conda environment**
```bash
conda env create -f environment.yml
conda activate codon_analyzer
```

### **3. Configure the workflow**
Edit the `config.yaml` file to set your parameters:
```yaml
# Input data configuration
input:
  type: "fasta"  # Options: fastq, fasta, cds, codonw
  files: "data/transcripts.fasta"
  
  # For FASTQ input
  # fastq_files: ["data/sample1_R1.fastq.gz", "data/sample1_R2.fastq.gz"]
  # reference: "reference/genome.fa"
  # gtf: "reference/annotation.gtf"

# Genetic code configuration
genetic_code:
  table: 1  # Standard genetic code
  # For custom genetic code, provide a dictionary:
  # custom_table:
  #   TTT: "F"
  #   TTC: "F"
  #   ...

# Analysis parameters
analysis:
  calculate_rscu: true
  calculate_cai: true
  calculate_enc: true
  calculate_tai: true
  correlation_analysis: true

# Output configuration
output:
  base_dir: "results"
  plots_dir: "results/plots"
  tables_dir: "results/tables"
  reports_dir: "results/reports"

# Resource configuration
resources:
  threads: 16
  memory: "32G"
```

### **4. Run the workflow**

#### **Local execution**
```bash
snakemake --cores 16 --use-conda
```

#### **Cluster execution (SLURM example)**
```bash
snakemake --cluster "sbatch --cpus-per-task={threads} --mem={resources.mem_mb}MB" \
          --jobs 10 --use-conda
```

## Project Structure

```
codon_usage_analyzer/
├── README.md                  # Project documentation
├── environment.yml            # Conda environment configuration
├── config.yaml                # Workflow configuration
├── Snakefile                  # Main Snakemake workflow
├── rules/                     # Snakemake rules
│   ├── preprocessing.smk      # Data preprocessing rules
│   ├── codon_calculation.smk  # Codon usage calculation rules
│   ├── bias_analysis.smk      # Codon bias analysis rules
│   ├── tai_analysis.smk       # tAI calculation rules
│   ├── visualization.smk      # Visualization rules
│   └── reporting.smk          # Report generation rules
├── scripts/                   # Analysis scripts
│   ├── *.py                   # Python analysis scripts
│   └── *.R                    # R visualization scripts
├── envs/                      # Conda environment files
│   ├── base.yaml              # Base environment
│   ├── bioinformatics.yaml    # Bioinformatics tools environment
│   └── r_analysis.yaml        # R analysis environment
├── reference/                 # Reference data directory
│   ├── genetic_codes/         # Genetic code tables
│   └── trna_data/             # tRNA gene data
├── data/                      # Input data directory
└── results/                   # Output results directory
```

## Workflow Details

### **Preprocessing**
1. **Quality Control**: FastQC for sequencing data
2. **Alignment**: STAR or HISAT2 for read alignment
3. **Assembly**: StringTie for transcript assembly
4. **CDS Extraction**: Extract coding sequences from annotations

### **Codon Usage Calculation**
- **Relative Synonymous Codon Usage (RSCU)**
- **Codon Adaptation Index (CAI)**
- **Effective Number of Codons (ENC)**
- **Codon Bias Index (CBI)**
- **Frequency of Optimal Codons (Fop)**

### **Codon Bias Analysis**
- **Preferred/Avoided Codons**: Identify optimal codons
- **Codon Pair Bias**: Analyze codon pair usage
- **Correlation with Gene Expression**: Relate codon usage to expression levels
- **Cross-Species Comparison**: Compare codon usage across organisms

### **tRNA Adaptation Index (tAI)**
- **tRNA Gene Copy Number**: Use known tRNA gene counts
- **tRNA Adaptation Calculation**: Compute tAI for each gene
- **Translation Efficiency Prediction**: Predict translation rates

### **Visualization**
- **Codon Usage Heatmaps**: Visualize codon preferences
- **Correlation Plots**: Relate codon usage to other metrics
- **PCA Analysis**: Dimensionality reduction of codon usage
- **Phylogenetic Trees**: Based on codon usage patterns

## Output Files

### **Main Outputs**
- `results/tables/codon_usage_table.tsv`: Comprehensive codon usage statistics
- `results/tables/rscu_values.tsv`: RSCU values for all codons
- `results/tables/cai_values.tsv`: CAI values for all genes
- `results/tables/enc_values.tsv`: ENC values for all genes
- `results/tables/preferred_codons.tsv`: Identified preferred codons
- `results/plots/codon_usage_heatmap.pdf`: Heatmap of codon usage
- `results/reports/final_report.html`: Comprehensive analysis report

### **Example Output**

**Codon Usage Table**:
```
Codon  AminoAcid  Count  Frequency  RSCU  CAI  ENC
TTT    F          1234   0.056      0.89  0.76  45.2
TTC    F          1567   0.071      1.11  0.89  45.2
TTA    L          890    0.040      0.32  0.23  52.1
...
```

## Configuration Guide

### **Input Configuration**
```yaml
input:
  # Input type: fastq, fasta, cds, or codonw
  type: "fasta"
  
  # For FASTQ input
  # fastq_files: ["data/sample1_R1.fastq.gz", "data/sample1_R2.fastq.gz"]
  # reference: "reference/genome.fa"
  # gtf: "reference/annotation.gtf"
  
  # For FASTA input  
  files: "data/transcripts.fasta"
  
  # For CDS input
  # cds_files: "data/cds_sequences.fasta"
  
  # For CodonW input
  # codonw_files: "data/codonw_output.blk"
```

### **Genetic Code Configuration**
```yaml
genetic_code:
  # Use standard genetic code table (1-25)
  table: 1
  
  # For custom genetic code
  # custom_table:
  #   TTT: "F"
  #   TTC: "F"
  #   TTA: "L"
  #   TTG: "L"
  #   ...
```

### **Analysis Parameters**
```yaml
analysis:
  # Basic codon usage metrics
  calculate_rscu: true
  calculate_cai: true
  calculate_enc: true
  calculate_cbi: true
  calculate_fop: true
  
  # Advanced analysis
  codon_pair_analysis: true
  correlation_analysis: true
  phylogenetic_analysis: false
  
  # tAI calculation
  calculate_tai: true
  trna_data: "reference/trna_data/human_trna_counts.tsv"
  
  # Visualization
  generate_heatmaps: true
  generate_correlation_plots: true
  generate_pca: true
```

## Advanced Usage

### **Adding Custom Genetic Codes**
1. Create a custom genetic code file in `reference/genetic_codes/`
2. Update the config file to use your custom table:
```yaml
genetic_code:
  custom_table: "reference/genetic_codes/my_custom_code.tsv"
```

### **Using tRNA Gene Data**
1. Prepare a tRNA gene count file with columns: Codon, tRNA_Copies
2. Update the config file:
```yaml
analysis:
  calculate_tai: true
  trna_data: "reference/trna_data/my_organism_trna.tsv"
```

### **Batch Processing**
For processing multiple samples, use wildcards in your config:
```yaml
input:
  type: "fastq"
  fastq_files: "data/{sample}_R{read}.fastq.gz"
  samples: ["sample1", "sample2", "sample3"]
  reads: [1, 2]
```

## Troubleshooting

### **Common Issues**

#### **Missing Dependencies**
```bash
# Recreate the conda environment
conda env create -f environment.yml
conda activate codon_analyzer

# Check if all prerequisites are installed
snakemake --list-conda-envs
```

#### **Memory Issues**
- Reduce the number of threads in `config.yaml`
- Increase the memory limit
- Use cluster execution with appropriate resource allocation

#### **Genetic Code Issues**
- Ensure you're using the correct genetic code table for your organism
- For custom genetic codes, verify the codon to amino acid mappings

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

- Original authors of codon usage analysis methods
- Contributors to the bioinformatics community
- Snakemake development team for the workflow management system

## Contact

For questions or issues, please:
  - Check the [GitHub Issues](https://github.com/autosomal/codon_usage_analyzer/issues) page
  - Contact the maintainer at phylogenetics@outlook.com
- Refer to the documentation for detailed parameter explanations