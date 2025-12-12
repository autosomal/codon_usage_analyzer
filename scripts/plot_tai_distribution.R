# Plot tAI distribution
# 
# This script creates distribution plots for tAI (tRNA Adaptation Index) values

library(ggplot2)
library(dplyr)
library(readr)

# Function to create tAI distribution plot
plot_tai_distribution <- function(tai_file, output_pdf, output_png) {
  # Read tAI data
  tai_data <- read_tsv(tai_file)
  
  # Create PDF plot
  p_pdf <- ggplot(tai_data, aes(x = tAI)) +
    geom_histogram(bins = 50, fill = "skyblue", color = "black", alpha = 0.7) +
    labs(
      title = "Distribution of tRNA Adaptation Index (tAI) Values",
      x = "tAI Value",
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
  p_png <- ggplot(tai_data, aes(x = tAI)) +
    geom_histogram(bins = 50, fill = "lightgreen", color = "black", alpha = 0.7) +
    labs(
      title = "Distribution of tRNA Adaptation Index (tAI) Values",
      x = "tAI Value",
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
  
  cat("tAI distribution plots saved to:\n")
  cat("- PDF: ", output_pdf, "\n")
  cat("- PNG: ", output_png, "\n")
}

# Main execution
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  # For testing purposes, use snakemake object if available
  if (exists("snakemake")) {
    tai_file <- snakemake@input[['tai_table']]
    output_pdf <- snakemake@output[['tai_plot']]
    output_png <- snakemake@output[['tai_histogram']]
  } else {
    stop("No input file provided and not running in Snakemake context")
  }
} else {
  tai_file <- args[1]
  output_pdf <- args[2]
  output_png <- args[3]
}

plot_tai_distribution(tai_file, output_pdf, output_png)