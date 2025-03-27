#!/usr/bin/env Rscript
# chi_square_variant_analysis.R
# This script calculates the proportion of variants with DAF < 0.01 for selected variant types, computes standard errors, and performs a chi-square test using provided rare/common counts. 
# Author: Alexandra Chapman
# Date: March 2025

# Load required libraries and set the threshold value
library(dplyr)
Threshold <- 0.01

# Subset 'variants' for the effects of interest
combined_variants <- subset(variants, Effect %in% c("missense_variant", "synonymous_variant", "stop_gained", 
                                                    "intron_variant", "intron_variant,non_coding_transcript_variant", 
                                                    "upstream_gene_variant", "downstream_gene_variant", 
                                                    "intergenic_variant", "non_coding_transcript_exon_variant", 
                                                    "3_prime_UTR_variant", "5_prime_UTR_variant", 
                                                    "intron_variant,NMD_transcript_variant"))

# Calculate variant proportions for each effect type:
variant_propotions <- combined_variants %>%
  group_by(Effect) %>%
  summarise(Proportion_DAF_less_0_01 = mean(DAF < Threshold)) %>%
  mutate(Effect = factor(Effect, levels = Effect[order(Proportion_DAF_less_0_01, decreasing = TRUE)]))
print("Variant proportions by effect:")
print(variant_propotions)

#For each variant group calculate sqrt((p(1-p))/n) where p is each Proportion_DAF_less_0_01 and n comes from nrow(variant_subset). These are the SE values.

# Manually assign standard errors for each variant type (make sure that the order matches the groups)
variant_propotions$SE <- c(0.00008737557, 0.0000789462, 0.0001033422, 0.0001036285, 
                           0.00004764548, 0.0000903296, 0.00007194132, 0.0001918408, 
                           0.00009319658, 0.0005681913, 0.0004086796, 0.00009367884)

# Multiply SE by 2 (as specified)
variant_propotions$SE <- 2 * variant_propotions$SE


# Create a table of rare and common variants from provided counts
rare_variants <- c(245364, 81262, 1829650, 815011, 8064987, 37222, 3682592, 209510, 258880, 4966, 103396, 2094080)
common_variants <- c(4134.997, 1110.998, 38027.02, 17037.01, 156414.6, 672.0008, 74536.93, 1654.005, 4987.007, 8.000182, 1885.004, 40775.09)
varianttable <- data.frame(Rare = rare_variants, Common = common_variants)
print(varianttable)

# Run chi-square test
chi_square_test <- chisq.test(varianttable)
print(chi_square_test)

#End of chi_square_variant_analysis.R

