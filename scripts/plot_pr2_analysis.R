#!/usr/bin/env Rscript
# PR2 (Parity Rule 2) Analysis Plot
# Analyzes the equality of nucleotides at 1st and 2nd codon positions

# Load required libraries
library(ggplot2)
library(dplyr)
library(readr)
library(gridExtra)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Input and output files
codon_file <- args[1]
output_pdf <- args[2]

# Validate inputs
if(length(args) < 2) {
  stop("Usage: Rscript plot_pr2_analysis.R <codon_usage_table.tsv> <output_pdf>")
}

# Read codon usage data
if(!file.exists(codon_file)) {
  stop(paste("Codon usage file does not exist:", codon_file))
}

codon_data <- read_tsv(codon_file, col_types = cols())

# Ensure we have the required columns
if(!"gene_id" %in% names(codon_data) || !"codon" %in% names(codon_data) || !"frequency" %in% names(codon_data)) {
  stop("Codon usage file must contain 'gene_id', 'codon', and 'frequency' columns")
}

# Function to calculate nucleotide frequencies at 1st and 2nd positions
calculate_pr2_frequencies <- function(codon_data) {
  # For each gene, calculate A, T, G, C frequencies at 1st and 2nd positions
  pr2_data <- codon_data %>%
    group_by(gene_id) %>%
    summarise(
      # Frequencies at 1st position
      a1_freq = sum(ifelse(str_sub(codon, 1, 1) == "A", frequency, 0)) / sum(frequency),
      t1_freq = sum(ifelse(str_sub(codon, 1, 1) == "T", frequency, 0)) / sum(frequency),
      g1_freq = sum(ifelse(str_sub(codon, 1, 1) == "G", frequency, 0)) / sum(frequency),
      c1_freq = sum(ifelse(str_sub(codon, 1, 1) == "C", frequency, 0)) / sum(frequency),
      
      # Frequencies at 2nd position
      a2_freq = sum(ifelse(str_sub(codon, 2, 2) == "A", frequency, 0)) / sum(frequency),
      t2_freq = sum(ifelse(str_sub(codon, 2, 2) == "T", frequency, 0)) / sum(frequency),
      g2_freq = sum(ifelse(str_sub(codon, 2, 2) == "G", frequency, 0)) / sum(frequency),
      c2_freq = sum(ifelse(str_sub(codon, 2, 2) == "C", frequency, 0)) / sum(frequency),
      .groups = 'drop'
    ) %>%
    filter(!is.na(a1_freq) & !is.na(t1_freq) & !is.na(g1_freq) & !is.na(c1_freq) &
           !is.na(a2_freq) & !is.na(t2_freq) & !is.na(g2_freq) & !is.na(c2_freq))
  
  return(pr2_data)
}

# Calculate PR2 frequencies
pr2_data <- calculate_pr2_frequencies(codon_data)

# Add AT and GC ratios for PR2 plot
pr2_data <- pr2_data %>%
  mutate(
    at1_ratio = (a1_freq + t1_freq) / (a1_freq + t1_freq + g1_freq + c1_freq),
    at2_ratio = (a2_freq + t2_freq) / (a2_freq + t2_freq + g2_freq + c2_freq),
    gc1_ratio = (g1_freq + c1_freq) / (a1_freq + t1_freq + g1_freq + c1_freq),
    gc2_ratio = (g2_freq + c2_freq) / (a2_freq + t2_freq + g2_freq + c2_freq),
    purine1_ratio = (a1_freq + g1_freq) / (a1_freq + t1_freq + g1_freq + c1_freq),
    purine2_ratio = (a2_freq + g2_freq) / (a2_freq + t2_freq + g2_freq + c2_freq),
    pyrimidine1_ratio = (t1_freq + c1_freq) / (a1_freq + t1_freq + g1_freq + c1_freq),
    pyrimidine2_ratio = (t2_freq + c2_freq) / (a2_freq + t2_freq + g2_freq + c2_freq)
  )

# Calculate correlations for PR2 analysis
at_correlation <- cor(pr2_data$a1_freq, pr2_data$t1_freq, use = "complete.obs")
gc_correlation <- cor(pr2_data$g1_freq, pr2_data$c1_freq, use = "complete.obs")
at2_correlation <- cor(pr2_data$a2_freq, pr2_data$t2_freq, use = "complete.obs")
gc2_correlation <- cor(pr2_data$g2_freq, pr2_data$c2_freq, use = "complete.obs")

# Create PR2 plots
p1 <- ggplot(pr2_data, aes(x = a1_freq, y = t1_freq)) +
  geom_point(alpha = 0.6, color = "steelblue", size = 1.5) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    title = paste("PR2 Plot: A vs T at 1st codon position\n(r =", round(at_correlation, 3), ")"),
    x = "Frequency of A at 1st position",
    y = "Frequency of T at 1st position"
  ) +
  xlim(0, 0.5) +
  ylim(0, 0.5) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12, face = "bold")
  )

p2 <- ggplot(pr2_data, aes(x = g1_freq, y = c1_freq)) +
  geom_point(alpha = 0.6, color = "forestgreen", size = 1.5) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    title = paste("PR2 Plot: G vs C at 1st codon position\n(r =", round(gc_correlation, 3), ")"),
    x = "Frequency of G at 1st position",
    y = "Frequency of C at 1st position"
  ) +
  xlim(0, 0.5) +
  ylim(0, 0.5) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12, face = "bold")
  )

p3 <- ggplot(pr2_data, aes(x = a2_freq, y = t2_freq)) +
  geom_point(alpha = 0.6, color = "orange", size = 1.5) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    title = paste("PR2 Plot: A vs T at 2nd codon position\n(r =", round(at2_correlation, 3), ")"),
    x = "Frequency of A at 2nd position",
    y = "Frequency of T at 2nd position"
  ) +
  xlim(0, 0.5) +
  ylim(0, 0.5) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12, face = "bold")
  )

p4 <- ggplot(pr2_data, aes(x = g2_freq, y = c2_freq)) +
  geom_point(alpha = 0.6, color = "purple", size = 1.5) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    title = paste("PR2 Plot: G vs C at 2nd codon position\n(r =", round(gc2_correlation, 3), ")"),
    x = "Frequency of G at 2nd position",
    y = "Frequency of C at 2nd position"
  ) +
  xlim(0, 0.5) +
  ylim(0, 0.5) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12, face = "bold")
  )

# Combine all PR2 plots
combined_pr2 <- grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2,
                             top = textGrob("PR2 Analysis: Parity Rule 2 Assessment", 
                                            gp = gpar(fontsize = 16, fontface = "bold")))

# Save PDF
pdf(output_pdf, width = 14, height = 12)
print(combined_pr2)
dev.off()

# Print summary statistics
n_genes <- nrow(pr2_data)

cat("PR2 Analysis Summary:\n")
cat("Number of genes analyzed:", n_genes, "\n")
cat("A vs T correlation at 1st position:", round(at_correlation, 4), "\n")
cat("G vs C correlation at 1st position:", round(gc_correlation, 4), "\n")
cat("A vs T correlation at 2nd position:", round(at2_correlation, 4), "\n")
cat("G vs C correlation at 2nd position:", round(gc2_correlation, 4), "\n")
cat("\nPR2 interpretation:\n")
cat("  - Points near diagonal line indicate compliance with Parity Rule 2\n")
cat("  - Deviations from diagonal indicate mutational or selective biases\n")
cat("  - PR2 plot helps identify directional mutational pressure\n")

cat("\nPR2 analysis plot saved successfully!\n")