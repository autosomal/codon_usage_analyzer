#!/usr/bin/env Rscript
# Generate heatmap of codon usage patterns

library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)

# Get input/output from Snakemake
codon_table_file <- snakemake@input[["codon_table"]]
rscu_table_file <- snakemake@input[["rscu_table"]]
output_pdf <- snakemake@output[["heatmap"]]
output_png <- snakemake@output[["heatmap_png"]]

# Get parameters
genetic_code <- snakemake@params[["genetic_code"]]
plot_formats <- snakemake@params[["plot_format"]]

# Read data
codon_df <- read.delim(codon_table_file, stringsAsFactors = FALSE)
rscu_df <- read.delim(rscu_table_file, stringsAsFactors = FALSE)

# Merge data
combined_df <- merge(codon_df, rscu_df[, c("codon", "rscu")], by = "codon", all.x = TRUE)

# Create codon order based on standard genetic code
standard_order <- c(
  # First base T
  "TTT", "TTC", "TTA", "TTG",
  "TCT", "TCC", "TCA", "TCG",
  "TAT", "TAC", "TAA", "TAG",
  "TGT", "TGC", "TGA", "TGG",
  
  # First base C
  "CTT", "CTC", "CTA", "CTG",
  "CCT", "CCC", "CCA", "CCG",
  "CAT", "CAC", "CAA", "CAG",
  "CGT", "CGC", "CGA", "CGG",
  
  # First base A
  "ATT", "ATC", "ATA", "ATG",
  "ACT", "ACC", "ACA", "ACG",
  "AAT", "AAC", "AAA", "AAG",
  "AGT", "AGC", "AGA", "AGG",
  
  # First base G
  "GTT", "GTC", "GTA", "GTG",
  "GCT", "GCC", "GCA", "GCG",
  "GAT", "GAC", "GAA", "GAG",
  "GGT", "GGC", "GGA", "GGG"
)

# Order the data
combined_df$codon <- factor(combined_df$codon, levels = standard_order)
combined_df <- combined_df[order(combined_df$codon), ]

# Prepare data for heatmap
heatmap_data <- combined_df %>%
  select(codon, amino_acid, rscu, frequency, count) %>%
  mutate(
    first_base = substr(codon, 1, 1),
    second_base = substr(codon, 2, 2),
    third_base = substr(codon, 3, 3)
  )

# Create RSCU heatmap
rscu_matrix <- matrix(heatmap_data$rscu, nrow = 4, ncol = 16, byrow = TRUE)
rownames(rscu_matrix) <- c("T", "C", "A", "G")
colnames(rscu_matrix) <- rep(c("T", "C", "A", "G"), 4)

# Create annotation for amino acids
aa_annotation <- matrix(heatmap_data$amino_acid, nrow = 4, ncol = 16, byrow = TRUE)

# Create color palette
my_palette <- colorRampPalette(brewer.pal(9, "YlOrRd"))(100)

# Create heatmap
pheatmap(rscu_matrix,
         color = my_palette,
         border_color = "white",
         cellwidth = 20,
         cellheight = 20,
         display_numbers = aa_annotation,
         number_color = "black",
         fontsize_number = 8,
         main = "Codon Usage Heatmap (RSCU Values)",
         filename = output_pdf,
         width = 8,
         height = 6)

# Also save as PNG if requested
if ("png" %in% plot_formats) {
  pheatmap(rscu_matrix,
           color = my_palette,
           border_color = "white",
           cellwidth = 20,
           cellheight = 20,
           display_numbers = aa_annotation,
           number_color = "black",
           fontsize_number = 8,
           main = "Codon Usage Heatmap (RSCU Values)",
           filename = output_png,
           width = 8,
           height = 6)
}

# Create frequency heatmap
frequency_matrix <- matrix(heatmap_data$frequency * 100, nrow = 4, ncol = 16, byrow = TRUE)

pheatmap(frequency_matrix,
         color = colorRampPalette(brewer.pal(9, "Blues"))(100),
         border_color = "white",
         cellwidth = 20,
         cellheight = 20,
         display_numbers = round(frequency_matrix, 1),
         number_color = "black",
         fontsize_number = 6,
         main = "Codon Frequency Heatmap (%)",
         filename = gsub("heatmap", "frequency_heatmap", output_pdf),
         width = 8,
         height = 6)

# Print summary
cat("Heatmap generation completed successfully\n")
cat(sprintf("RSCU heatmap saved to: %s\n", output_pdf))
if ("png" %in% plot_formats) {
  cat(sprintf("RSCU heatmap (PNG) saved to: %s\n", output_png))
}
cat(sprintf("Analyzed %d codons\n", nrow(combined_df)))