#!/usr/bin/env Rscript
# plot_dnds_gene_sets.R
# Description: This script loads dN/dS data for several gene sets, calculates the dn/ds values for Mouse and 
Macaque orthologs and computes the average dn/ds for each gene set, and then creates plots (Mouse and Macaque) showing the average dN/dS values by gene set.

# Author: Alexandra Chapman
# Date: March 2025

# Load required libraries
library(data.table) 
library(dplyr)
library(tidyr)      
library(ggplot2)


# Function: Calculate dn/ds for a Given Gene Set
# Calculate dn/ds for the specified species ("Mouse" or "Macaque")
 # The function checks for missing values or extreme dS values and returns NA in those cases.

calculate_dn_ds <- function(dt, species_prefix) {
  dt[[paste0("dn_ds_", species_prefix)]] <- with(dt, 
    ifelse(is.na(get(paste0("dN with ", species_prefix))) | 
           is.na(get(paste0("dS with ", species_prefix))) | 
           get(paste0("dS with ", species_prefix)) == 0 | 
           get(paste0("dS with ", species_prefix)) > 2, 
           NA, 
           get(paste0("dN with ", species_prefix)) / get(paste0("dS with ", species_prefix))
    ))
  return(dt)
}

# Load and Process Each Gene Set File from BioMart output files
# Same as below for Binding genes, Transcription Factor (TF) genes, Ion Channel genes, Transporter genes, Structural genes and Transducer genes
# This is for catalytic genes
Catalytic_genes <- fread("Catalytic_dnds.txt")
Catalytic_genes <- calculate_dn_ds(Catalytic_genes, "Mouse")
Catalytic_genes <- calculate_dn_ds(Catalytic_genes, "Macaque")

# Calculate Average dN/dS per Gene Set for Mouse and Macaque

Mouse_orthologues <- data.frame(
  Catalytic_Activity      = mean(Catalytic_genes$dn_ds_Mouse, na.rm = TRUE),
  Transcription_Factor    = mean(Tf_genes$dn_ds_Mouse, na.rm = TRUE),
  Binding_Activity        = mean(Binding_genes$dn_ds_Mouse, na.rm = TRUE),
  Ion_channel_Activity    = mean(Ion_genes$dn_ds_Mouse, na.rm = TRUE),
  Transporter_Activity    = mean(Transporter_genes$dn_ds_Mouse, na.rm = TRUE),
  Structural_Molecules    = mean(Structural_genes$dn_ds_Mouse, na.rm = TRUE),
  Molecular_transducers   = mean(Transducer_genes$dn_ds_Mouse, na.rm = TRUE)
)

Macaque_orthologues <- data.frame(
  Catalytic_Activity      = mean(Catalytic_genes$dn_ds_Macaque, na.rm = TRUE),
  Transcription_Factor    = mean(Tf_genes$dn_ds_Macaque, na.rm = TRUE),
  Binding_Activity        = mean(Binding_genes$dn_ds_Macaque, na.rm = TRUE),
  Ion_channel_Activity    = mean(Ion_genes$dn_ds_Macaque, na.rm = TRUE),
  Transporter_Activity    = mean(Transporter_genes$dn_ds_Macaque, na.rm = TRUE),
  Structural_Molecules    = mean(Structural_genes$dn_ds_Macaque, na.rm = TRUE),
  Molecular_transducers   = mean(Transducer_genes$dn_ds_Macaque, na.rm = TRUE)
)

cat("Mouse Orthologues dN/dS Values:\n")
print(Mouse_orthologues)
cat("\nMacaque Orthologues dN/dS Values:\n")
print(Macaque_orthologues)

# Prepare Data for Plotting
# The dataframes were was pivotted Long Format for plotting

# For Mouse:
df_mouse <- Mouse_orthologues %>%
  pivot_longer(cols = everything(), names_to = "GeneSet", values_to = "dNdS_value")
df_mouse$GeneSet <- factor(df_mouse$GeneSet, levels = df_mouse$GeneSet[order(df_mouse$dNdS_value, decreasing = TRUE)])

# For Macaque:
df_macaque <- Macaque_orthologues %>%
  pivot_longer(cols = everything(), names_to = "GeneSet", values_to = "dNdS_value")
df_macaque$GeneSet <- factor(df_macaque$GeneSet, levels = df_macaque$GeneSet[order(df_macaque$dNdS_value, decreasing = TRUE)])

# Create the Plots for Mouse and Macaque
# Define a common color for gene sets
colors <- c(
  "Transcription_Factor"  = "red3", 
  "Ion_channel_Activity"  = "yellow3", 
  "Transporter_Activity"  = "green3", 
  "Structural_Molecules"  = "turquoise", 
  "Catalytic_Activity"    = "skyblue2", 
  "Binding_Activity"      = "purple", 
  "Molecular_transducers" = "orchid"
)

# Mouse plot
Mouse_plot <- ggplot(df_mouse, aes(x = GeneSet, y = dNdS_value, color = GeneSet)) +
  geom_point(size = 3) +
  scale_color_manual(values = colors) +
  theme_minimal() +
  labs(x = "Gene Set", y = "Average dN/dS Value", title = "Average dN/dS Values for Mouse") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"))
print(Mouse_plot)

# Macaque plot
Macaque_plot <- ggplot(df_macaque, aes(x = GeneSet, y = dNdS_value, color = GeneSet)) +
  geom_point(size = 3) +
  scale_color_manual(values = colors) +
  theme_minimal() +
  labs(x = "Gene Set", y = "Average dN/dS Value", title = "Average dN/dS Values for Macaque") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"))
print(Macaque_plot)

#End of plot_dnds_gene_sets.R
