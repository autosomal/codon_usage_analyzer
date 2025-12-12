# Plot correlation between codon bias metrics
# 
# This script creates correlation plots between different codon bias metrics:
# CAI, ENC, tAI, CBI, FOP

library(ggplot2)
library(dplyr)
library(readr)
library(corrplot)
library(reshape2)

# Function to create correlation plots between bias metrics
plot_bias_correlation <- function(
  cai_file, enc_file, tai_file, cbi_file, fop_file, 
  output_pdf, output_matrix
) {
  # Read all bias metric files
  cai_data <- read_tsv(cai_file)
  enc_data <- read_tsv(enc_file)
  tai_data <- read_tsv(tai_file)
  cbi_data <- read_tsv(cbi_file)
  fop_data <- read_tsv(fop_file)
  
  # Merge all data by Gene
  combined_data <- cai_data %>%
    inner_join(enc_data, by = "Gene") %>%
    inner_join(tai_data, by = "Gene") %>%
    inner_join(cbi_data, by = "Gene") %>%
    inner_join(fop_data, by = "Gene")
  
  # Select only the metric columns for correlation
  metrics_data <- combined_data %>% select(CAI, ENC, tAI, CBI, FOP)
  
  # Calculate correlation matrix
  cor_matrix <- cor(metrics_data, use = "complete.obs")
  
  # Save correlation matrix as TSV
  write.table(cor_matrix, output_matrix, sep = "\t", row.names = TRUE, col.names = NA)
  
  # Create correlation plot
  pdf(output_pdf, width = 10, height = 8)
  
  # Create the correlation plot
  corrplot(cor_matrix, 
           method = "color",
           type = "upper",
           order = "hclust",
           tl.cex = 0.8,
           tl.col = "black",
           tl.srt = 45,
           addCoef.col = "black",
           number.cex = 0.7,
           title = "Correlation Between Codon Bias Metrics",
           mar = c(0,0,1,0))
  
  dev.off()
  
  cat("Correlation analysis completed:\n")
  cat("- Correlation plot saved to: ", output_pdf, "\n")
  cat("- Correlation matrix saved to: ", output_matrix, "\n")
  
  # Print correlation values to console
  cat("\nCorrelation Matrix:\n")
  print(cor_matrix)
}

# Main execution
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  # For testing purposes, use snakemake object if available
  if (exists("snakemake")) {
    cai_file <- snakemake@input[['cai_table']]
    enc_file <- snakemake@input[['enc_table']]
    tai_file <- snakemake@input[['tai_table']]
    cbi_file <- snakemake@input[['cbi_table']]
    fop_file <- snakemake@input[['fop_table']]
    output_pdf <- snakemake@output[['correlation_plot']]
    output_matrix <- snakemake@output[['correlation_matrix']]
  } else {
    stop("No input files provided and not running in Snakemake context")
  }
} else {
  cai_file <- args[1]
  enc_file <- args[2]
  tai_file <- args[3]
  cbi_file <- args[4]
  fop_file <- args[5]
  output_pdf <- args[6]
  output_matrix <- args[7]
}

plot_bias_correlation(cai_file, enc_file, tai_file, cbi_file, fop_file, 
                      output_pdf, output_matrix)