#!/usr/bin/env Rscript
# plot_gene_set_analysis.R
# Description: This script performs gene set analysis for a specified variant type. It calculates the proportion of variants with DAF < Threshold for each gene set, computes standard errors, and generates a plot with error bars.

# Author: Alexandra Chapman
# Date: March 2025

# Load required libraries
library(ggplot2)
library(dplyr)

variant_type <- "stop_gained"  # or whichever variant type you want to analyze

# Set threshold value
Threshold <- 0.01

# Subset the 'variants' data frame for your chosen variant type
variant_subset <- subset(variants, Effect %in% c(variant_type))

# Initialize empty vectors to store gene set names, proportions, and number of variants
gene_set_names <- c()
proportion <- c()
num_vars <- c()

# Loop through each gene set column in 'genesets'
# This assumes that the first column of 'genesets' is 'Gene'
for(i in 2:ncol(genesets)) {
  gene_set_names <- c(gene_set_names, colnames(genesets)[i])
  current_geneset <- genesets$Gene[genesets[, i] == 1]
  x <- subset(variant_subset, Gene %in% current_geneset)
  
  # Calculate the proportion of variants with DAF < Threshold
  # Use nrow(x) as n; if no variants exist, handle division by zero.
  prop_val <- if(nrow(x) > 0) {
    length(x$DAF[x$DAF < Threshold]) / nrow(x)
  } else {
    NA  # Or 0, depending on your preference
  }
  proportion <- c(proportion, prop_val)
  
  # Record the number of variants in this gene set
  num_vars <- c(num_vars, nrow(x))
}

# Create a dataframs for results
results <- data.frame(
  GeneSet = gene_set_names,
  Proportion_DAF_below_Threshold = proportion,
  Num_Variants = num_vars,
  stringsAsFactors = FALSE
)

# Order results by decreasing proportion
results <- results[order(-results$Proportion_DAF_below_Threshold), ]
rownames(results) <- NULL

# Compute standard error for each gene set:
# SE = 2 * sqrt( p * (1 - p) / n )
results$SE <- 2 * sqrt(results$Proportion_DAF_below_Threshold * (1 - results$Proportion_DAF_below_Threshold) / results$Num_Variants)

# Define the desired order for gene sets (optional. This was used in analysis to order all plots in the same order as the missense plot)
desired_order <- c(
  "Transcription_factors", 
  "Ion_channel_activity", 
  "Transporter_activity", 
  "Structural_molecules", 
  "Catalytic_activity", 
  "Binding_activity", 
  "Molecular_transducers"
)
# Set the order of the factor in the results data
results$GeneSet <- factor(results$GeneSet, levels = desired_order)

# Create the plot
p <- ggplot(results, aes(x = GeneSet, y = Proportion_DAF_below_Threshold)) +
  geom_point(aes(color = GeneSet), size = 2) +
  geom_errorbar(aes(
    ymin = Proportion_DAF_below_Threshold - SE,
    ymax = Proportion_DAF_below_Threshold + SE
  ), width = 0.1, color = "black") +
  scale_color_manual(values = setNames(
    c("red3", "yellow3", "green3", "turquoise", "skyblue2", "purple", "orchid"),
    desired_order
  )) +
  labs(
    title = paste("Proportion of", variant_type, "Variants with DAF < 0.01\nbased on Gene Type"),
    x = "Gene Type",
    y = "Proportion"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),  # Hide x-axis text (adjust if needed)
    plot.title = element_text(family = "serif", face = "bold", size = 17),
    legend.title = element_text(family = "serif", face = "bold", size = 14),
    axis.title.x = element_text(family = "serif", face = "bold", size = 12),
    axis.title.y = element_text(family = "serif", face = "bold", size = 12),
    legend.position = "none"
  )

# Print the plot
print(p)


#End of plot_gene_set_analysis.R 
