#!/bin/bash
# run_vep.sh
# This script runs Ensembl's Variant Effect Predictor (VEP) on a VCF file
# Adjust the paths and file names as needed.
# Author = Alexandra Chapman
# Date = March 2025

# Path to VEP executable (change this to actual VEP path)
VEP_PATH="/path/to/vep"

# Input VCF file containing the variants you wish to annotate
INPUT_VCF="input_variants.vcf"

# Reference FASTA file for ancestral sequence 
REFERENCE_FASTA="/path/to/reference.fa.gz”

# Annotation file (GFF) for gene models
GFF_FILE="/path/to/annotations.gff3"

# Output file for VEP results
OUTPUT_FILE="vep_output.txt"

# Run VEP with the desired parameters
$VEP_PATH/vep \
  -i $INPUT_VCF \
  --gff $GFF_FILE \
  --fasta $REFERENCE_FASTA \
  -o $OUTPUT_FILE

# End of script
