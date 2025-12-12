#!/usr/bin/env Rscript
# Codon Usage Effectiveness Analysis
# Analyzes codon usage patterns and calculates various effectiveness measures

# Load required libraries
library(ggplot2)
library(dplyr)
library(readr)
library(gridExtra)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Input and output files
codon_file <- args[1]
cai_file <- args[2]
output_pdf <- args[3]

# Validate inputs
if(length(args) < 3) {
  stop("Usage: Rscript plot_codon_usage_effectiveness.R <codon_usage_table.tsv> <cai_values.tsv> <output_pdf>")
}

# Read codon usage data
if(!file.exists(codon_file)) {
  stop(paste("Codon usage file does not exist:", codon_file))
}

# Read CAI data
if(!file.exists(cai_file)) {
  stop(paste("CAI file does not exist:", cai_file))
}

codon_data <- read_tsv(codon_file, col_types = cols())
cai_data <- read_tsv(cai_file, col_types = cols())

# Ensure we have the required columns
if(!all(c("gene_id", "codon", "frequency") %in% names(codon_data))) {
  stop("Codon usage file must contain 'gene_id', 'codon', and 'frequency' columns")
}

if(!all(c("gene_id", "cai_value") %in% names(cai_data))) {
  stop("CAI file must contain 'gene_id' and 'cai_value' columns")
}

# Merge codon usage and CAI data
merged_data <- merge(codon_data, cai_data, by = "gene_id")

# Calculate additional metrics
gene_summary <- merged_data %>%
  group_by(gene_id, cai_value) %>%
  summarise(
    total_codons = n(),
    codon_diversity = n_distinct(codon),
    mean_frequency = mean(frequency),
    max_frequency = max(frequency),
    gc_content = sum(ifelse(str_sub(codon, 3, 3) %in% c("G", "C"), frequency, 0)) / sum(frequency),
    .groups = 'drop'
  ) %>%
  filter(!is.na(cai_value))

# Calculate correlation between CAI and various metrics
cor_cai_codon_div <- cor(gene_summary$cai_value, gene_summary$codon_diversity, use = "complete.obs")
cor_cai_gc <- cor(gene_summary$cai_value, gene_summary$gc_content, use = "complete.obs")

# Create plots
p1 <- ggplot(gene_summary, aes(x = cai_value, y = codon_diversity)) +
  geom_point(alpha = 0.6, color = "steelblue", size = 1.5) +
  geom_smooth(method = "lm", color = "red") +
  labs(
    title = paste("CAI vs Codon Diversity\n(r =", round(cor_cai_codon_div, 3), ")"),
    x = "Codon Adaptation Index (CAI)",
    y = "Number of Different Codons"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12, face = "bold")
  )

p2 <- ggplot(gene_summary, aes(x = cai_value, y = gc_content)) +
  geom_point(alpha = 0.6, color = "forestgreen", size = 1.5) +
  geom_smooth(method = "loess", color = "red") +
  labs(
    title = paste("CAI vs GC3 Content\n(r =", round(cor_cai_gc, 3), ")"),
    x = "Codon Adaptation Index (CAI)",
    y = "GC3 Content"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12, face = "bold")
  )

# Distribution of CAI values
p3 <- ggplot(gene_summary, aes(x = cai_value)) +
  geom_histogram(bins = 50, fill = "orange", color = "darkorange", alpha = 0.7) +
  labs(
    title = "Distribution of CAI Values",
    x = "Codon Adaptation Index (CAI)",
    y = "Gene Count"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12, face = "bold")
  )

# Create boxplot of CAI by codon diversity categories
gene_summary$codon_div_cat <- cut(gene_summary$codon_diversity, 
                                  breaks = quantile(gene_summary$codon_diversity, probs = seq(0, 1, 0.25)), 
                                  include.lowest = TRUE, labels = c("Low", "Medium-Low", "Medium-High", "High"))

p4 <- ggplot(gene_summary, aes(x = codon_div_cat, y = cai_value)) +
  geom_boxplot(fill = "purple", alpha = 0.7) +
  geom_jitter(alpha = 0.2, width = 0.1) +
  labs(
    title = "CAI Distribution by Codon Diversity Category",
    x = "Codon Diversity Category",
    y = "Codon Adaptation Index (CAI)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    axis.text = element_text(size = 11, angle = 45, hjust = 1),
    axis.title = element_text(size = 12, face = "bold")
  )

# Combine plots
combined_plot <- grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2,
                              top = textGrob("Codon Usage Effectiveness Analysis", 
                                             gp = gpar(fontsize = 16, fontface = "bold")))

# Save PDF
pdf(output_pdf, width = 14, height = 12)
print(combined_plot)
dev.off()

# Print summary statistics
n_genes <- nrow(gene_summary)
mean_cai <- mean(gene_summary$cai_value, na.rm = TRUE)
median_cai <- median(gene_summary$cai_value, na.rm = TRUE)
sd_cai <- sd(gene_summary$cai_value, na.rm = TRUE)

cat("Codon Usage Effectiveness Analysis Summary:\n")
cat("Number of genes analyzed:", n_genes, "\n")
cat("Mean CAI:", round(mean_cai, 4), "\n")
cat("Median CAI:", round(median_cai, 4), "\n")
cat("SD of CAI:", round(sd_cai, 4), "\n")
cat("Correlation between CAI and codon diversity:", round(cor_cai_codon_div, 4), "\n")
cat("Correlation between CAI and GC3 content:", round(cor_cai_gc, 4), "\n")

cat("\nCodon usage effectiveness plot saved successfully!\n")