#!/usr/bin/env Rscript
# plot_variant_proportions.R
# This script reads in the variant proportions and standard errors (saved as a .csv file)  and generates a plot with error bars showing the proportion of variants with DAF < 0.01.
# Author: Alexandra Chapman
# Date: March 2025

# Load required libraries
library(ggplot2)
library(dplyr)

# Read in the CSV file produced by the statistical analysis script
variant_propotions <- read.csv("variant_propotions.csv", stringsAsFactors = FALSE)

# Ensure 'Effect' is ordered by decreasing proportion
variant_propotions$Effect <- factor(variant_propotions$Effect, 
levels = variant_propotions$Effect[order(variant_propotions$Proportion_DAF_less_0_01, decreasing = TRUE)])

# Generate the plot with error bars
variant_propotions_plot <- ggplot(variant_propotions, aes(x = Effect, y = Proportion_DAF_less_0_01)) +
  geom_point(aes(color = Effect), size = 2) +
  geom_errorbar(aes(ymin = Proportion_DAF_less_0_01 - SE,
                    ymax = Proportion_DAF_less_0_01 + SE),
                width = 0.1, color = "black") +
  labs(title = "Proportion of Variants with DAF < 0.01",
       x = "Variant Effect Class",
       y = "Proportion") +
  theme_bw() +
  theme(
    plot.title = element_text(family = "serif", face = "bold", size = 18),
    legend.title = element_text(family = "serif", face = "bold", size = 14),
    axis.title.x = element_text(family = "serif", face = "bold", size = 12),
    axis.title.y = element_text(family = "serif", face = "bold", size = 12)
  )
print(variant_propotions_plot)

#End of plot_variant_proportions.R
