# Plot CBI distribution
# 
# This script creates distribution plots for CBI (Codon Bias Index) values

library(ggplot2)
library(dplyr)
library(readr)

# Function to create CBI distribution plot
plot_cbi_distribution <- function(cbi_file, output_pdf, output_png) {
  # Read CBI data
  cbi_data <- read_tsv(cbi_file)
  
  # Create PDF plot
  p_pdf <- ggplot(cbi_data, aes(x = CBI)) +
    geom_histogram(bins = 50, fill = "coral", color = "black", alpha = 0.7) +
    labs(
      title = "Distribution of Codon Bias Index (CBI) Values",
      x = "CBI Value",
      y = "Frequency"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
  
  # Save PDF
  ggsave(output_pdf, plot = p_pdf, width = 8, height = 6, device = "pdf")
  
  # Create PNG plot
  p_png <- ggplot(cbi_data, aes(x = CBI)) +
    geom_histogram(bins = 50, fill = "orange", color = "black", alpha = 0.7) +
    labs(
      title = "Distribution of Codon Bias Index (CBI) Values",
      x = "CBI Value",
      y = "Frequency"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
  
  # Save PNG
  ggsave(output_png, plot = p_png, width = 8, height = 6, device = "png", dpi = 300)
  
  cat("CBI distribution plots saved to:\n")
  cat("- PDF: ", output_pdf, "\n")
  cat("- PNG: ", output_png, "\n")
}

# Main execution
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  # For testing purposes, use snakemake object if available
  if (exists("snakemake")) {
    cbi_file <- snakemake@input[['cbi_table']]
    output_pdf <- snakemake@output[['cbi_plot']]
    output_png <- snakemake@output[['cbi_histogram']]
  } else {
    stop("No input file provided and not running in Snakemake context")
  }
} else {
  cbi_file <- args[1]
  output_pdf <- args[2]
  output_png <- args[3]
}

plot_cbi_distribution(cbi_file, output_pdf, output_png)