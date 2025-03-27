# BI453-Research-Project
Inference of Evolutionary Selection Pressure acting on Rare Genomic Variants from Derived Allele Frequency Spectra in the UK Biobank

Author: Alexandra Chapman
Course: BI453 
Academic Year: 2024/25

# Overview
This repository contains the code for my final-year project. This project investigates evolutionary selection pressures acting on rare genomic variants by analyzing DAF spectra using data from the UK Biobank. By utilizing Ensembl’s VEP for functional annotation and by analyzing their frequency patterns, this research seeks to identify signals of purifying selection. This study examined selection pressure across different functional variant types including missense, synonymous and stop-gained mutations and compared patterns across different genomic regions. The dN/dS ratios from Mouse and Macaque orthologs were used as a validation metric to confirm that an excess of low-frequency derived alleles provides evidence of purifying selection. Statistical analyses were also conducted to determine significant trends and associations between rare allele frequencies and their functional impacts.

# Repository Folder Overview:

# Data_Analysis
- calculate_variant_stats.R summarizes key variant statistics (e.g., the total number of variants, proportions of rare variants, unique annotation types etc.) by reading in the processed chromosome data and merging it into a combined dataset.
  
# Statistical_Analysis
- chi_square_variant_analysis.R performs chi-squared tests to assess associations between variant categories (e.g., rare vs common variants) and different genomic or functional classes.

# Validation
- validate_daf_dnds_combined.R correlates derived allele frequency (DAF) measures (average or proportion below a threshold) with dN/dS values for genes in Mouse and Macaque orthologs. This script generates combined plots overlaying both of these relationships.
  
# Variant_annotation
- run_vep.sh is a shell script that runs Ensembl's Variant Effect Predictor (VEP) on local variant files to annotate variants with predicted functional consequences (e.g., missense, synonymous).
  
# Visualization
- plot_DAF_spectra_by_effect.R creates histograms (DAF spectra plots) showing the distribution of derived allele frequencies for different variant effects (e.g., missense, synonymous, stop-gained).
  
- plot_gene_set_analysis.R analyzes specific gene sets (extracted from BioMart using Gene Ontology (GO) terms to show how variants in each set differ in proportions of rare alleles.
  
- plot_variant_frequencies.R produces bar plots comparing the counts of various variant categories.
  
- plot_dnds_gene_sets.R shows the average dN/dS values across multiple GO gene sets for Mouse and Macaque orthologs, visualizing potential evolutionary constraints.
  
- plot_gene_set_thresholds.R creates plots illustrating how different thresholds for rare variants (DAF < T) impact the proportions of rare variants within specific gene sets.
  
- plot_variant_proportions.R plots the proportion of variants below a certain DAF threshold (e.g., 0.01) across different variant types, with error bars representing confidence intervals.

# Further Information
This repository contains only the generalised code and scripts used for this project. Raw data are not included due to data-sharing restrictions. Detailed results and the full project write-up are also not publically available in this repository. 

If you have any questions about the methods, results, or wish to view the full project documentation, please contact me directly at A.Chapman3@universityofgalway.ie.
