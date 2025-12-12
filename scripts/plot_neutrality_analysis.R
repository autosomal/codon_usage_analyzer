#!/usr/bin/env Rscript
# Neutrality Analysis Plot
# Generates GC12 vs GC3 plot for neutrality analysis in codon usage

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
  stop("Usage: Rscript plot_neutrality_analysis.R <codon_usage_table.tsv> <output_pdf>")
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

# Function to calculate GC content at different positions
calculate_gc_content <- function(codon_data) {
  # For each gene, calculate GC1, GC2, and GC3 content
  gc_content <- codon_data %>%
    group_by(gene_id) %>%
    summarise(
      gc1 = sum(ifelse(str_sub(codon, 1, 1) %in% c("G", "C"), frequency, 0)) / sum(frequency),
      gc2 = sum(ifelse(str_sub(codon, 2, 2) %in% c("G", "C"), frequency, 0)) / sum(frequency),
      gc3 = sum(ifelse(str_sub(codon, 3, 3) %in% c("G", "C"), frequency, 0)) / sum(frequency),
      gc12 = (sum(ifelse(str_sub(codon, 1, 1) %in% c("G", "C"), frequency, 0)) + 
              sum(ifelse(str_sub(codon, 2, 2) %in% c("G", "C"), frequency, 0))) / 
             (2 * sum(frequency)),
      .groups = 'drop'
    ) %>%
    filter(!is.na(gc1) & !is.na(gc2) & !is.na(gc3))
  
  return(gc_content)
}

# Calculate GC content
gc_data <- calculate_gc_content(codon_data)

# Calculate correlation between GC12 and GC3
correlation <- cor(gc_data$gc12, gc_data$gc3, use = "complete.obs")
n_genes <- nrow(gc_data)

# Create the main neutrality plot (GC12 vs GC3)
p1 <- ggplot(gc_data, aes(x = gc12, y = gc3)) +
  geom_point(alpha = 0.6, color = "steelblue", size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  labs(
    title = paste("Neutrality Plot: GC12 vs GC3 Content\n(n =", n_genes, "genes, r =", round(correlation, 3), ")"),
    subtitle = "Used to assess the role of mutation bias vs. selection in codon usage",
    x = "GC content at 1st & 2nd positions (GC12)",
    y = "GC content at 3rd position (GC3)"
  ) +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13, face = "bold")
  )

# Create histogram of GC12 content
p2 <- ggplot(gc_data, aes(x = gc12)) +
  geom_histogram(bins = 50, fill = "lightgreen", color = "darkgreen", alpha = 0.7) +
  labs(
    title = "Distribution of GC12 Content",
    x = "GC12 Content",
    y = "Gene Count"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13, face = "bold")
  )

# Create histogram of GC3 content
p3 <- ggplot(gc_data, aes(x = gc3)) +
  geom_histogram(bins = 50, fill = "orange", color = "darkorange", alpha = 0.7) +
  labs(
    title = "Distribution of GC3 Content",
    x = "GC3 Content",
    y = "Gene Count"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13, face = "bold")
  )

# Create hexagonal binning plot for better visualization of density
p4 <- ggplot(gc_data, aes(x = gc12, y = gc3)) +
  geom_hex(bins = 30, show.legend = FALSE) +
  scale_fill_gradient(low = "white", high = "navy") +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  labs(
    title = "Neutrality Plot (Hexagonal Binning)",
    x = "GC12 Content",
    y = "GC3 Content"
  ) +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13, face = "bold")
  )

# Combine plots
combined_plot <- grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2,
                              top = textGrob("Neutrality Analysis: Mutation vs Selection in Codon Usage", 
                                             gp = gpar(fontsize = 16, fontface = "bold")))

# Save PDF
pdf(output_pdf, width = 14, height = 12)
print(combined_plot)
dev.off()

# Print summary statistics
slope <- coef(lm(gc3 ~ gc12, data = gc_data))[2]
intercept <- coef(lm(gc3 ~ gc12, data = gc_data))[1]

cat("Neutrality Analysis Summary:\n")
cat("Number of genes analyzed:", n_genes, "\n")
cat("Correlation between GC12 and GC3:", round(correlation, 4), "\n")
cat("Regression slope:", round(slope, 4), "\n")
cat("Regression intercept:", round(intercept, 4), "\n")
cat("Interpretation:\n")
if(correlation > 0.5) {
  cat("  - High correlation suggests mutation bias is the primary force\n")
} else if(correlation < 0.3) {
  cat("  - Low correlation suggests selection is the primary force\n")
} else {
  cat("  - Moderate correlation suggests both mutation bias and selection are at play\n")
}

cat("\nNeutrality analysis plot saved successfully!\n")