#!/usr/bin/env Rscript
# plot_gene_set_thresholds.R
# Description: This script calculates the proportion of variants with DAF < each of several thresholds for multiple variant types (e.g., Synonymous, Stop_gained, 5'UTR, 3'UTR, Intron, Upstream, Downstream). It computes standard errors, combines the results, and then generates a plot comparing the variant types.

# Author: Alexandra Chapman
# Date: March 2025

# Load required libraries
library(ggplot2)
library(dplyr)

# Define a vector of thresholds
thresholds <- c(0.00001, 0.00005, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1)

# Define a function that calculates proportions and standard errors for a given variant subset
calculate_proportions <- function(variant_subset, variant_name, thresholds) {
  df <- data.frame()
  for (thresh in thresholds) {
    num_variants <- sum(variant_subset$DAF < thresh)
    proportion <- if (nrow(variant_subset) > 0) mean(variant_subset$DAF < thresh) else NA
    standard_error <- if(num_variants > 0) 2 * sqrt(proportion * (1 - proportion) / num_variants) else NA
    
    # Append the results
    df <- rbind(df, data.frame(
      Threshold = format(thresh, scientific = FALSE),
      Proportion_DAF_less_Threshold = proportion,
      Num_Variants = num_variants,
      SE = standard_error,
      stringsAsFactors = FALSE
    ))
  }
  # Adds a column for the variant type
  df$Variant_Type <- variant_name
  return(df)
}

# Calculate proportions for each variant type assuming that the following objects exist in your environment): 'synonymous', 'stop_gained', 'fiveprime', 'threeprime', 'upstream', 'downstream', 'intron’.
synonymous_proportions <- calculate_proportions(synonymous, "Synonymous", thresholds)
stop_proportions       <- calculate_proportions(stop_gained, "Stop_gained", thresholds)
fiveprime_proportions  <- calculate_proportions(fiveprime, "5'UTR", thresholds)
threeprime_proportions <- calculate_proportions(threeprime, "3'UTR", thresholds)
upstream_proportions   <- calculate_proportions(upstream, "Upstream", thresholds)
downstream_proportions <- calculate_proportions(downstream, "Downstream", thresholds)
intron_proportions     <- calculate_proportions(intron, "Intron", thresholds)
missense_proportions   <- calculate_proportions(missense, "Missense", thresholds)

# Combine all the proportions into a single data frame
combined_proportions <- rbind(
  synonymous_proportions,
  stop_proportions,
  fiveprime_proportions,
  threeprime_proportions,
  upstream_proportions,
  downstream_proportions,
  intron_proportions,
 missense_proportions  
)

# Convert 'Threshold' to a factor and define desired order
combined_proportions$Threshold <- factor(combined_proportions$Threshold, levels = unique(combined_proportions$Threshold))
combined_proportions$Variant_Type <- factor(combined_proportions$Variant_Type,
                                            levels = c("Stop_gained", “Missense”, "5'UTR", "3'UTR", "Synonymous", "Intron", "Upstream", "Downstream"))

# Create the combined plot
combined_plot <- ggplot(combined_proportions, aes(x = Threshold, y = Proportion_DAF_less_Threshold, color = Variant_Type)) +
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = Proportion_DAF_less_Threshold - SE, ymax = Proportion_DAF_less_Threshold + SE),
                width = 0.005, color = "black") +
  labs(title = "Proportion of Variants with DAF < Varying Thresholds per Variant Type",
       x = "Threshold",
       y = "Proportion") +
  theme_bw() +
  scale_color_manual(values = c(
    "Missense"    = "#DE8C00", 
    "Synonymous"  = "#00C08B",
    "Stop_gained" = "#F8766D",
    "5'UTR"       = "#B79F00",
    "3'UTR"       = "#7CAE00",
    "Upstream"    = "#619CFF",
    "Downstream"  = "#F564E3",
    "Intron"      = "#00B4F0"
  )) +
  theme(
    plot.title = element_text(face = "bold", size = 9),
    axis.text.x = element_text(family = "serif", size = 10),
    axis.title.x = element_text(family = "serif", size = 10),
    axis.title.y = element_text(family = "serif", size = 10)
  )

# Display the plot
print(combined_plot)


#End of plot_gene_set_thresholds.R
