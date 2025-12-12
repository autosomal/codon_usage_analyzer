#!/usr/bin/env Rscript
# Plot RSCU Correlation
# Generates correlation plots for Relative Synonymous Codon Usage values

# Load required libraries
library(ggplot2)
library(dplyr)
library(readr)
library(corrplot)
library(gridExtra)
library(reshape2)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Input and output files
rscu_file <- args[1]
output_pdf <- args[2]
output_matrix_pdf <- args[3]

# Validate inputs
if(length(args) < 3) {
  stop("Usage: Rscript plot_rscu_correlation.R <rscu_values.tsv> <output_pdf> <output_matrix_pdf>")
}

# Read RSCU data
if(!file.exists(rscu_file)) {
  stop(paste("RSCU file does not exist:", rscu_file))
}

rscu_data <- read_tsv(rscu_file, col_types = cols())

# Ensure we have the required columns
if(!"gene_id" %in% names(rscu_data) || !"codon" %in% names(rscu_data) || !"rscu_value" %in% names(rscu_data)) {
  stop("RSCU file must contain 'gene_id', 'codon', and 'rscu_value' columns")
}

# Pivot the data to wide format for correlation analysis
rscu_wide <- rscu_data %>%
  select(gene_id, codon, rscu_value) %>%
  filter(!is.na(rscu_value)) %>%
  dcast(gene_id ~ codon, value.var = "rscu_value")

# Remove gene_id column for correlation calculation
gene_ids <- rscu_wide$gene_id
rscu_wide <- rscu_wide[, -1]  # Remove gene_id column

# Calculate correlation matrix
cor_matrix <- cor(rscu_wide, use = "pairwise.complete.obs")

# Filter out codons with too many NAs for meaningful correlations
complete_cols <- colSums(!is.na(rscu_wide)) > nrow(rscu_wide) * 0.1  # At least 10% of genes should have values
cor_matrix <- cor_matrix[complete_cols, complete_cols]

# Create the main correlation plot
pdf(output_pdf, width = 14, height = 12)
corrplot(cor_matrix, method = "color", type = "upper", order = "hclust",
         tl.cex = 0.8, tl.col = "black", tl.srt = 45,
         title = "RSCU Correlation Matrix", mar = c(0,0,3,0),
         addCoef.col = "black", number.cex = 0.6)
title(main = "Correlation Between RSCU Values Across Codons", line = -2, cex.main = 1.2)
dev.off()

# Create a more detailed correlation matrix plot with additional information
pdf(output_matrix_pdf, width = 16, height = 14)
# Enhanced correlation plot
par(mar=c(5, 4, 4, 8) + 0.1)  # Add space for right margin labels
corrplot(cor_matrix, method = "color", type = "full", order = "hclust",
         tl.cex = 0.7, tl.col = "black", tl.srt = 45,
         col = colorRampPalette(c("blue", "white", "red"))(200),
         addCoef.col = "black", number.cex = 0.5, diag = TRUE,
         title = "RSCU Correlation Matrix")
title(main = "Correlation Between RSCU Values Across Codons", cex.main = 1.5)
dev.off()

# Create scatter plots for highly correlated codon pairs
high_corr_pairs <- which(abs(cor_matrix) > 0.7 & abs(cor_matrix) < 1.0, arr.ind = TRUE)
high_corr_pairs <- high_corr_pairs[high_corr_pairs[,1] < high_corr_pairs[,2], ]  # Remove duplicates

if(nrow(high_corr_pairs) > 0) {
  # Get top 6 highly correlated pairs
  corr_values <- cor_matrix[high_corr_pairs]
  top_pairs_idx <- order(abs(corr_values), decreasing = TRUE)[1:min(6, length(corr_values))]
  top_pairs <- high_corr_pairs[top_pairs_idx, ]
  
  # Create scatter plots for top correlated pairs
  plots_list <- list()
  for(i in 1:nrow(top_pairs)) {
    codon1 <- colnames(cor_matrix)[top_pairs[i, 1]]
    codon2 <- colnames(cor_matrix)[top_pairs[i, 2]]
    corr_val <- round(cor_matrix[top_pairs[i, 1], top_pairs[i, 2]], 3)
    
    plot_data <- data.frame(
      codon1_values = rscu_wide[[codon1]],
      codon2_values = rscu_wide[[codon2]]
    )
    
    p <- ggplot(plot_data, aes_string(x = "codon1_values", y = "codon2_values")) +
      geom_point(alpha = 0.6, color = "steelblue") +
      geom_smooth(method = "lm", se = TRUE, color = "red") +
      labs(
        title = paste("RSCU Correlation:", codon1, "vs", codon2),
        subtitle = paste("r =", corr_val),
        x = paste("RSCU", codon1),
        y = paste("RSCU", codon2)
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 10),
        plot.subtitle = element_text(hjust = 0.5, size = 9),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9)
      )
    
    plots_list[[i]] <- p
  }
  
  # Save scatter plots if we have any
  if(length(plots_list) > 0) {
    # Calculate how many rows/columns we need for the layout
    n_plots <- length(plots_list)
    n_col <- ifelse(n_plots <= 2, 1, ifelse(n_plots <= 4, 2, 3))
    n_row <- ceiling(n_plots / n_col)
    
    # Create combined plot
    combined_scatter <- do.call(grid.arrange, c(plots_list, ncol = n_col))
    
    # Save to a separate file
    scatter_output <- gsub("\\.[^.]*$", "_scatter.pdf", output_pdf)
    pdf(scatter_output, width = 8 * n_col, height = 6 * n_row)
    print(combined_scatter)
    dev.off()
  }
}

# Print summary statistics
cat("RSCU Correlation Analysis Summary:\n")
cat("Number of codons analyzed:", ncol(rscu_wide), "\n")
cat("Number of genes analyzed:", nrow(rscu_wide), "\n")
cat("Number of codon pairs with correlation > 0.7:", nrow(high_corr_pairs), "\n")
cat("Average absolute correlation:", round(mean(abs(cor_matrix[lower.tri(cor_matrix)]), na.rm = TRUE), 4), "\n")

cat("\nRSCU correlation plots saved successfully!\n")