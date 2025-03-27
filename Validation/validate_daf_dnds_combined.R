#!/usr/bin/env Rscript
# validate_daf_dnds_combined.R
# This script calculates two measures for missense variant genes:
#  (1) The average DAF per gene and the corresponding average dN/dS value.
#  (2) The proportion of DAF values below a threshold (0.01) per gene and the corresponding average dN/dS value.
#   It then creates a combined plot for each species (Mouse and Macaque) that overlays the two relationships.

#   This script assumes that the objects 'missense' (columns "Gene" and "DAF") and 'Missense_genes' (dN/dS values in columns "dn_ds_Mouse" and "dn_ds_Macaque") are already loaded into your R environment.

# Author: Alexandra Chapman
# Date: March 2025

library(ggplot2)
library(dplyr)

# Set the threshold for low DAF
T <- 0.01

# Function to calculate measures per gene and generate a combined plot for a given species.
# 'species' should be either "Mouse" or "Macaque".
analyze_species <- function(species, missense, Missense_genes, threshold = 0.01) {
genes <- unique(missense_data$Gene)
  
  # Initialize data frames for average DAF and proportion of DAF < threshold.
  result_avg <- data.frame(Gene = genes, DAF = rep(NA, length(genes)), dnds = rep(NA, length(genes)), stringsAsFactors = FALSE)
  result_prop <- data.frame(Gene = genes, DAF = rep(NA, length(genes)), dnds = rep(NA, length(genes)), stringsAsFactors = FALSE)
  rownames(result_avg) <- genes
  rownames(result_prop) <- genes
  
  # Loop through each gene to compute the measures.
  for(g in genes) {
    # Average DAF for the gene.
    x <- missense_data$DAF[missense_data$Gene == g]
    result_avg[g, "DAF"] <- mean(x)
    
    # Proportion of DAF values below the threshold for the gene.
    result_prop[g, "DAF"] <- mean(x < threshold)
    
    # Use the appropriate column from gene_data.
    col_name <- paste0("dn_ds_", species)
    y <- gene_data[[col_name]][gene_data$Gene == g]
    result_avg[g, "dnds"] <- mean(y)
    result_prop[g, "dnds"] <- mean(y)
  }
  
# Analyze for Mouse.
mouse_results <- analyze_species("Mouse", missense, Missense_genes, threshold = threshold)
# Analyze for Macaque.
macaque_results <- analyze_species("Macaque", missense, Missense_genes, threshold = threshold)


  # Create the combined plot
  # proportion measure are plotted in one color and average measure in another
  # Both are fitted with a linear model
ggplot() +
    geom_point(data = result_prop, aes(x = DAF, y = dnds),
               color = if (species == "Mouse") "steelblue" else "mediumorchid", size = 2) +
    geom_smooth(data = result_prop, aes(x = DAF, y = dnds),
                method = "lm", se = TRUE, color = "grey20", fullrange = TRUE) +
    geom_point(data = result_avg, aes(x = DAF, y = dnds),
               color = if (species == "Mouse") "seagreen" else "deeppink", size = 2) +
    geom_smooth(data = result_avg, aes(x = DAF, y = dnds),
                method = "lm", se = TRUE, color = "grey20", fullrange = TRUE) +
    theme_minimal() +
    labs(x = "DAF", y = "dN/dS", title = paste("Correlations between DAF and dN/dS in", species, "Genes")) +
    coord_cartesian(ylim = if (species == "Mouse") c(0, 1.25) else c(0, 2.5))
  
  # Return a list with both results and the plot.
  return(list(avg = result_avg, prop = result_prop, plot = p))
}

#End of validate_daf_dnds_combined.R
