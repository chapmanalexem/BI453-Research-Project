#!/usr/bin/env Rscript
# calculate_variant_stats.R
# Description: This script reads in chromosome-specific variant data from TSV files,
# combines them, and calculates various statistics (e.g., total variants, variants per chromosome,
# unique annotation types, rare variant proportions, and gene counts).
# Author: Alexandra Chapman
# Date: March 2025

# 1: Load data

# Load required libraries
library(readr)
library(dplyr)

# Read in chromosome data from TSV files 
# The TSV files had the following columns: Gene, Effect, DAF
chr2  <- readr::read_tsv("Chr2_final.tsv")
chr4  <- readr::read_tsv("Chr4_final.tsv")
chr6  <- readr::read_tsv("Chr6_final.tsv")
chr8  <- readr::read_tsv("Chr8_final.tsv")
chr10 <- readr::read_tsv("Chr10_final.tsv")
chr12 <- readr::read_tsv("Chr12_final.tsv")
chr14 <- readr::read_tsv("Chr14_final.tsv")
chr16 <- readr::read_tsv("Chr16_final.tsv")
chr18 <- readr::read_tsv("Chr18_final.tsv")
chr20 <- readr::read_tsv("Chr20_final.tsv")
chr22 <- readr::read_tsv("Chr22_final.tsv")

# Combine all chromosomes into one data frame
combined_variants <- rbind(chr2, chr4, chr6, chr8, chr10, chr12, chr14, chr16, chr18, chr20, chr22)
print(combined_variants)

# 2: Calculate Variant Statistics 

# Total number of variants across all chromosomes
total_variants <- nrow(combined_variants)

# Number of unique annotation (Effect) types across all variants
num_annotation_types <- combined_variants %>%
  distinct(Effect) %>%
  nrow()

# Overall proportion (%) of variants with DAF < 0.01
prop_rare_variants <- combined_variants %>%
  filter(!is.na(DAF)) %>%
  summarise(prop = (sum(DAF < 0.01) / n()) * 100) %>%
  pull(prop)

# 3: Gene Analysis 

# Remove version numbers from Gene IDs 
variants <- combined_variants
variants$Gene <- sub("\\..*", "", variants$Gene)

# Total number of unique genes
total_genes <- variants %>%
  distinct(Gene) %>%
  nrow()

# 4: Unique Effect Types 

# List all unique effect types
effect_types <- combined_variants %>%
  distinct(Effect) %>%
  pull(Effect)

# Sorted list of effect types
sorted_effect_types <- sort(effect_types)
num_effects <- length(effect_types)

# 5: Print Results 

cat("Total number of variants:", total_variants, "\n")
cat("Number of unique annotation types:", num_annotation_types, "\n")
cat("Proportion of rare variants (DAF < 0.01):", round(prop_rare_variants, 3), "%\n")
cat("Total number of unique genes:", total_genes, "\n")
cat("Total number of unique effect types:", num_effects, "\n")
print(sorted_effect_types)


#End of calculate_variant_stats.R 
