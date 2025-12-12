# Plot FOP distribution
# 
# This script creates distribution plots for FOP (Frequency of Optimal Codons) values

library(ggplot2)
library(dplyr)
library(readr)

# Function to create FOP distribution plot
plot_fop_distribution <- function(fop_file, output_pdf, output_png) {
  # Read FOP data
  fop_data <- read_tsv(fop_file)
  
  # Create PDF plot
  p_pdf <- ggplot(fop_data, aes(x = FOP)) +
    geom_histogram(bins = 50, fill = "lightblue", color = "black", alpha = 0.7) +
    labs(
      title = "Distribution of Frequency of Optimal Codons (FOP) Values",
      x = "FOP Value",
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
  p_png <- ggplot(fop_data, aes(x = FOP)) +
    geom_histogram(bins = 50, fill = "lightgreen", color = "black", alpha = 0.7) +
    labs(
      title = "Distribution of Frequency of Optimal Codons (FOP) Values",
      x = "FOP Value",
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
  
  cat("FOP distribution plots saved to:\n")
  cat("- PDF: ", output_pdf, "\n")
  cat("- PNG: ", output_png, "\n")
}

# Main execution
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  # For testing purposes, use snakemake object if available
  if (exists("snakemake")) {
    fop_file <- snakemake@input[['fop_table']]
    output_pdf <- snakemake@output[['fop_plot']]
    output_png <- snakemake@output[['fop_histogram']]
  } else {
    stop("No input file provided and not running in Snakemake context")
  }
} else {
  fop_file <- args[1]
  output_pdf <- args[2]
  output_png <- args[3]
}

plot_fop_distribution(fop_file, output_pdf, output_png)