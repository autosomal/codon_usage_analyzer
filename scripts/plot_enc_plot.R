#!/usr/bin/env Rscript
# Plot ENC vs GC3 Content
# Generates the classic ENC-plot for codon usage bias analysis

# Load required libraries
library(ggplot2)
library(dplyr)
library(readr)
library(gridExtra)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Input and output files
enc_file <- args[1]
codon_file <- args[2]
output_pdf <- args[3]

# Validate inputs
if(length(args) < 3) {
  stop("Usage: Rscript plot_enc_plot.R <enc_values.tsv> <codon_usage_table.tsv> <output_pdf>")
}

# Read ENC data
if(!file.exists(enc_file)) {
  stop(paste("ENC file does not exist:", enc_file))
}

enc_data <- read_tsv(enc_file, col_types = cols())

# Read codon usage data
if(!file.exists(codon_file)) {
  stop(paste("Codon usage file does not exist:", codon_file))
}

codon_data <- read_tsv(codon_file, col_types = cols())

# Ensure we have the required columns
if(!"gene_id" %in% names(enc_data) || !"enc_value" %in% names(enc_data)) {
  stop("ENC file must contain 'gene_id' and 'enc_value' columns")
}

if(!"gene_id" %in% names(codon_data) || !"codon" %in% names(codon_data) || !"frequency" %in% names(codon_data)) {
  stop("Codon usage file must contain 'gene_id', 'codon', and 'frequency' columns")
}

# Calculate GC3 content for each gene
# First, identify codons that encode amino acids with synonymous codons
# Focus on 4-fold degenerate sites (ending in T, C, A, G) where third position can vary
# For simplicity, we'll extract GC3 content from codon usage data

# Function to calculate GC3 content
calculate_gc3 <- function(gene_codons) {
  # Filter for codons that end in the third position (synonymous codons)
  # We'll focus on codons that code for amino acids with multiple codons
  
  # Define codons that have synonymous variants (4-fold degenerate sites)
  # These are codons where the third position can vary without changing amino acid
  # For example: Gly (GGN), Ala (GCN), Thr (ACN), Pro (CCN), etc.
  
  # For this implementation, we'll calculate GC content at 3rd position for all codons
  third_pos_codons <- gene_codons %>%
    mutate(third_base = str_sub(codon, 3, 3)) %>%
    filter(third_base %in% c("G", "C", "A", "T"))
  
  # Calculate GC3 content
  gc3_content <- sum(third_pos_codons$frequency[third_pos_codons$third_base %in% c("G", "C")]) /
                 sum(third_pos_codons$frequency)
  
  return(gc3_content)
}

# Calculate GC3 content for each gene
# Group codon data by gene and calculate GC3
gc3_by_gene <- codon_data %>%
  group_by(gene_id) %>%
  summarise(
    gc3_content = sum(ifelse(str_sub(codon, 3, 3) %in% c("G", "C"), frequency, 0)) /
                 sum(frequency),
    .groups = 'drop'
  ) %>%
  filter(!is.na(gc3_content))

# Merge ENC and GC3 data
plot_data <- merge(enc_data, gc3_by_gene, by = "gene_id")
plot_data <- plot_data[complete.cases(plot_data), ]

# Calculate expected ENC values based on GC3 content
# From Wright (1990) formula: ENC = 2 + GC3 + 29/(GC3^2 + (1-GC3)^2)
calculate_expected_enc <- function(gc3) {
  ifelse(gc3 == 0 | gc3 == 1, 2 + 29, 2 + gc3 + 29/(gc3^2 + (1-gc3)^2))
}

# Add expected ENC curve
gc3_seq <- seq(0, 1, 0.01)
expected_enc <- calculate_expected_enc(gc3_seq)

# Create the ENC plot
p1 <- ggplot(plot_data, aes(x = gc3_content, y = enc_value)) +
  geom_point(alpha = 0.6, color = "steelblue", size = 1.5) +
  geom_line(data = data.frame(gc3_seq = gc3_seq, expected_enc = expected_enc),
            aes(x = gc3_seq, y = expected_enc), color = "red", linewidth = 1.2) +
  labs(
    title = "ENC Plot: Effective Number of Codons vs GC3 Content",
    subtitle = "Red line represents theoretical null expectation under no selection",
    x = "GC3 Content",
    y = "Effective Number of Codons (ENC)"
  ) +
  xlim(0, 1) +
  ylim(0, 62) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13, face = "bold")
  )

# Create histogram of ENC values
p2 <- ggplot(plot_data, aes(x = enc_value)) +
  geom_histogram(bins = 50, fill = "lightgreen", color = "darkgreen", alpha = 0.7) +
  labs(
    title = "Distribution of ENC Values",
    x = "Effective Number of Codons (ENC)",
    y = "Gene Count"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13, face = "bold")
  )

# Create histogram of GC3 content
p3 <- ggplot(plot_data, aes(x = gc3_content)) +
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

# Create correlation plot between ENC and GC3
p4 <- ggplot(plot_data, aes(x = gc3_content, y = enc_value)) +
  geom_hex(bins = 30, show.legend = FALSE) +
  scale_fill_gradient(low = "white", high = "navy") +
  geom_smooth(method = "loess", color = "red", se = FALSE) +
  labs(
    title = "ENC vs GC3 Content (Hexagonal Binning)",
    x = "GC3 Content",
    y = "Effective Number of Codons (ENC)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13, face = "bold")
  )

# Combine plots
combined_plot <- grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2,
                              top = textGrob("ENC Analysis: Codon Usage Bias Assessment", 
                                             gp = gpar(fontsize = 16, fontface = "bold")))

# Save PDF
pdf(output_pdf, width = 14, height = 12)
print(combined_plot)
dev.off()

# Calculate summary statistics
n_genes <- nrow(plot_data)
mean_enc <- mean(plot_data$enc_value, na.rm = TRUE)
mean_gc3 <- mean(plot_data$gc3_content, na.rm = TRUE)
correlation <- cor(plot_data$enc_value, plot_data$gc3_content, use = "complete.obs")

# Print summary statistics
cat("ENC-Plot Analysis Summary:\n")
cat("Number of genes analyzed:", n_genes, "\n")
cat("Mean ENC value:", round(mean_enc, 4), "\n")
cat("Mean GC3 content:", round(mean_gc3, 4), "\n")
cat("Correlation between ENC and GC3:", round(correlation, 4), "\n")
cat("Expected ENC range (theoretical): 2-62\n")
cat("Lower observed ENC values indicate stronger codon usage bias\n")

cat("\nENC plot saved successfully!\n")