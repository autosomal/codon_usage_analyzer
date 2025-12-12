#!/usr/bin/env Rscript
# Plot Codon Usage PCA
# Generates PCA plots for codon usage patterns

# Load required libraries
library(ggplot2)
library(dplyr)
library(readr)
library(FactoMineR)
library(factoextra)
library(RColorBrewer)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Input and output files
codon_file <- args[1]
output_pdf <- args[2]
group_var <- ifelse(length(args) > 2, args[3], "none")

# Validate inputs
if(length(args) < 2) {
  stop("Usage: Rscript plot_codon_pca.R <codon_usage_table.tsv> <output_pdf> [group_variable]")
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

# Pivot the data to wide format for PCA analysis
codon_wide <- codon_data %>%
  select(gene_id, codon, frequency) %>%
  filter(!is.na(frequency)) %>%
  tidyr::pivot_wider(names_from = codon, values_from = frequency, values_fill = 0)

# Extract gene IDs and remove the gene_id column for PCA
gene_ids <- codon_wide$gene_id
codon_wide <- codon_wide[, -1]  # Remove gene_id column

# Handle missing values by replacing with 0
codon_wide[is.na(codon_wide)] <- 0

# Check if we have enough codons for meaningful analysis
if(ncol(codon_wide) < 2) {
  stop("Not enough codons with valid frequencies for PCA analysis")
}

if(nrow(codon_wide) < 3) {
  stop("Not enough genes for PCA analysis")
}

# Perform PCA
pca_result <- PCA(as.matrix(codon_wide), scale.unit = TRUE, ncp = 5, graph = FALSE)

# Extract coordinates for plotting
pca_coords <- pca_result$ind$coord
pca_df <- data.frame(
  gene_id = gene_ids,
  PC1 = pca_coords[, 1],
  PC2 = pca_coords[, 2],
  PC3 = pca_coords[, 3]
)

# If grouping variable is specified and exists in original data, add it
if(group_var != "none" && group_var %in% names(codon_data)) {
  # Get unique group assignments per gene
  group_assignments <- codon_data %>%
    select(gene_id, all_of(group_var)) %>%
    distinct() %>%
    filter(!is.na(.data[[group_var]]))
  
  # Merge with PCA data
  pca_df <- merge(pca_df, group_assignments, by = "gene_id", all.x = TRUE)
}

# Get variance explained by each component
variance_explained <- pca_result$eig[, 2]  # Percentage of variance
pc1_var <- round(variance_explained[1], 2)
pc2_var <- round(variance_explained[2], 2)

# Create PCA plot
if(group_var != "none" && group_var %in% names(pca_df)) {
  p1 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = !!sym(group_var))) +
    geom_point(alpha = 0.7, size = 2) +
    labs(
      title = paste("PCA of Codon Usage Patterns\nColored by", group_var),
      subtitle = paste("PC1 explains", pc1_var, "% variance, PC2 explains", pc2_var, "% variance"),
      x = paste("Principal Component 1 (", pc1_var, "%)"),
      y = paste("Principal Component 2 (", pc2_var, "%)")
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 13, face = "bold"),
      legend.title = element_text(size = 12, face = "bold")
    ) +
    guides(color = guide_legend(override.aes = list(size = 3)))
} else {
  p1 <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
    geom_point(alpha = 0.7, color = "steelblue", size = 2) +
    labs(
      title = "PCA of Codon Usage Patterns",
      subtitle = paste("PC1 explains", pc1_var, "% variance, PC2 explains", pc2_var, "% variance"),
      x = paste("Principal Component 1 (", pc1_var, "%)"),
      y = paste("Principal Component 2 (", pc2_var, "%)")
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 13, face = "bold")
    )
}

# Create scree plot showing variance explained by each component
variance_df <- data.frame(
  PC = 1:min(10, ncol(codon_wide)),
  Variance = variance_explained[1:min(10, ncol(codon_wide))]
)

p2 <- ggplot(variance_df, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "orange", alpha = 0.7) +
  labs(
    title = "Scree Plot: Variance Explained by Each PC",
    x = "Principal Component",
    y = "Percentage of Variance Explained (%)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13, face = "bold")
  )

# Create cumulative variance plot
cumvar_df <- data.frame(
  PC = 1:length(variance_explained),
  CumulativeVariance = cumsum(variance_explained)
)

p3 <- ggplot(cumvar_df[1:min(10, length(variance_explained)), ], 
             aes(x = PC, y = CumulativeVariance)) +
  geom_line(color = "purple", linewidth = 1.2) +
  geom_point(color = "purple", size = 2) +
  labs(
    title = "Cumulative Variance Explained",
    x = "Number of Principal Components",
    y = "Cumulative Variance Explained (%)"
  ) +
  ylim(0, 100) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13, face = "bold")
  )

# Create PC3 vs PC1 plot if we have enough components
if(length(variance_explained) >= 3) {
  pc3_var <- round(variance_explained[3], 2)
  p4 <- ggplot(pca_df, aes(x = PC1, y = PC3)) +
    geom_point(alpha = 0.7, color = "darkgreen", size = 2) +
    labs(
      title = "PCA: PC1 vs PC3",
      subtitle = paste("PC1 explains", pc1_var, "% variance, PC3 explains", pc3_var, "% variance"),
      x = paste("Principal Component 1 (", pc1_var, "%)"),
      y = paste("Principal Component 3 (", pc3_var, "%)")
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 13, face = "bold")
    )
} else {
  # If we don't have enough PCs, create a simple histogram of PC1 values
  p4 <- ggplot(pca_df, aes(x = PC1)) +
    geom_histogram(bins = 30, fill = "lightblue", color = "darkblue", alpha = 0.7) +
    labs(
      title = "Distribution of PC1 Scores",
      x = "Principal Component 1 Score",
      y = "Count"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 13, face = "bold")
    )
}

# Combine plots
combined_plot <- gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2,
                              top = textGrob("PCA Analysis of Codon Usage Patterns", 
                                             gp = gpar(fontsize = 16, fontface = "bold")))

# Save PDF
pdf(output_pdf, width = 14, height = 12)
print(combined_plot)
dev.off()

# Print summary statistics
cat("PCA Analysis Summary:\n")
cat("Number of genes analyzed:", nrow(codon_wide), "\n")
cat("Number of codons analyzed:", ncol(codon_wide), "\n")
cat("Total variance explained by first 2 PCs:", round(pc1_var + pc2_var, 2), "%\n")
cat("Proportion of variance explained by PC1:", round(pc1_var, 2), "%\n")
cat("Proportion of variance explained by PC2:", round(pc2_var, 2), "%\n")

cat("\nPCA plot saved successfully!\n")