#!/usr/bin/env Rscript
# Plot CAI Distribution
# Generates publication-quality plots for Codon Adaptation Index (CAI) distribution

# Load required libraries
library(ggplot2)
library(dplyr)
library(readr)
library(gridExtra)
library(scales)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Input and output files
cai_file <- args[1]
output_pdf <- args[2]
output_png <- args[3]

# Validate inputs
if(length(args) < 2) {
  stop("Usage: Rscript plot_cai_distribution.R <cai_values.tsv> <output_pdf> [output_png]")
}

# Read CAI data
if(!file.exists(cai_file)) {
  stop(paste("CAI file does not exist:", cai_file))
}

cai_data <- read_tsv(cai_file, col_types = cols())

# Ensure we have the required columns
if(!"gene_id" %in% names(cai_data) || !"cai_value" %in% names(cai_data)) {
  stop("CAI file must contain 'gene_id' and 'cai_value' columns")
}

# Clean data - remove NA values
cai_data <- cai_data %>%
  filter(!is.na(cai_value))

# Summary statistics
n_genes <- nrow(cai_data)
mean_cai <- mean(cai_data$cai_value, na.rm = TRUE)
median_cai <- median(cai_data$cai_value, na.rm = TRUE)
sd_cai <- sd(cai_data$cai_value, na.rm = TRUE)

# Create distribution plot
p1 <- ggplot(cai_data, aes(x = cai_value)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white", alpha = 0.7) +
  labs(
    title = paste("Distribution of Codon Adaptation Index (CAI)\n(n =", n_genes, "genes)"),
    subtitle = paste("Mean:", round(mean_cai, 3), "| Median:", round(median_cai, 3), "| SD:", round(sd_cai, 3)),
    x = "Codon Adaptation Index (CAI)",
    y = "Gene Count"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13, face = "bold")
  )

# Create density plot
p2 <- ggplot(cai_data, aes(x = cai_value)) +
  geom_density(fill = "lightblue", color = "darkblue", alpha = 0.6) +
  labs(
    title = "Density Plot of Codon Adaptation Index (CAI)",
    x = "Codon Adaptation Index (CAI)",
    y = "Density"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13, face = "bold")
  )

# Create violin plot if there are grouping variables
# For now, just create a basic version
if("gene_class" %in% names(cai_data)) {
  p3 <- ggplot(cai_data, aes(x = factor(1), y = cai_value)) +
    geom_violin(fill = "orange", alpha = 0.6) +
    geom_boxplot(width = 0.1, outlier.shape = NA) +
    labs(
      title = "Boxplot of CAI Values",
      x = "",
      y = "Codon Adaptation Index (CAI)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 13, face = "bold"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
} else {
  # Simple boxplot without grouping
  p3 <- ggplot(cai_data, aes(y = cai_value)) +
    geom_boxplot(fill = "orange", alpha = 0.6, width = 0.3) +
    coord_flip() +
    labs(
      title = "Boxplot of CAI Values",
      x = "Codon Adaptation Index (CAI)",
      y = ""
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 13, face = "bold")
    )
}

# Combine plots
combined_plot <- grid.arrange(p1, p2, p3, ncol = 2, nrow = 2,
                              top = textGrob("Codon Adaptation Index (CAI) Distribution Analysis", 
                                             gp = gpar(fontsize = 16, fontface = "bold")))

# Save PDF
pdf(output_pdf, width = 12, height = 10)
print(combined_plot)
dev.off()

# Save PNG if specified
if(length(args) >= 3 && !is.na(output_png) && output_png != "") {
  png(output_png, width = 12, height = 10, units = "in", res = 300)
  print(combined_plot)
  dev.off()
}

# Print summary statistics
cat("CAI Distribution Summary:\n")
cat("Number of genes:", n_genes, "\n")
cat("Mean CAI:", round(mean_cai, 4), "\n")
cat("Median CAI:", round(median_cai, 4), "\n")
cat("Standard deviation:", round(sd_cai, 4), "\n")
cat("Range:", round(min(cai_data$cai_value, na.rm = TRUE), 4), "-", 
    round(max(cai_data$cai_value, na.rm = TRUE), 4), "\n")

cat("\nPlots saved successfully!\n")