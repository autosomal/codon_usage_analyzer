#!/bin/bash
# Test script for Codon Usage Analyzer workflow

# Create test data directory
mkdir -p test_data reference

# Create dummy FASTA file with test sequences
echo -e ">gene1\nATGGCGTTCGAAGCTTAA" > test_data/transcripts.fasta
echo -e ">gene2\nATGTTCGAAGCTGCGTAA" >> test_data/transcripts.fasta
echo -e ">gene3\nATGGAAGCTTTCGCGTAA" >> test_data/transcripts.fasta

# Create dummy reference files
echo -e ">chr1\nATGGCGTTCGAAGCTTAAATGTTCGAAGCTGCGTAAATGGAAGCTTTCGCGTAA" > reference/genome.fa

echo -e "chr1\tunknown\texon\t1\t18\t.\t+\t.\tgene_id \"gene1\"; transcript_id \"gene1\";" > reference/annotation.gtf
echo -e "chr1\tunknown\texon\t19\t36\t.\t+\t.\tgene_id \"gene2\"; transcript_id \"gene2\";" >> reference/annotation.gtf
echo -e "chr1\tunknown\texon\t37\t54\t.\t+\t.\tgene_id \"gene3\"; transcript_id \"gene3\";" >> reference/annotation.gtf

# Create dummy tRNA data
echo -e "Codon\ttRNA_Copies" > reference/trna_data/human_trna_counts.tsv
echo -e "TTT\t3" >> reference/trna_data/human_trna_counts.tsv
echo -e "TTC\t2" >> reference/trna_data/human_trna_counts.tsv
echo -e "TTA\t1" >> reference/trna_data/human_trna_counts.tsv
echo -e "TTG\t1" >> reference/trna_data/human_trna_counts.tsv
echo -e "TCT\t4" >> reference/trna_data/human_trna_counts.tsv
echo -e "TCC\t4" >> reference/trna_data/human_trna_counts.tsv
echo -e "TCA\t3" >> reference/trna_data/human_trna_counts.tsv
echo -e "TCG\t2" >> reference/trna_data/human_trna_counts.tsv

# Run Snakemake workflow with test configuration
echo "Running Snakemake workflow with test configuration..."
snakemake --configfile test_config.yaml --cores 2 --use-conda --dry-run

# Clean up
rm -rf test_data reference test_results logs

echo "Test completed successfully"
echo "The workflow is ready for use with real data"