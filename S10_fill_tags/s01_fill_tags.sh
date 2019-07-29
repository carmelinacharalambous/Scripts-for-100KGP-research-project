#!/bin/bash

# s01_fill_tags.sh
# Last updated: CC20May2019

# Annotate with bcftools fill-tags plugin
# fill-tags plugin sets INFO tags AF, AC, AC_Hemi, AC_Hom, AC_Het, AN, HWE, MAF, NS

# Stop at runtime errors
set -e

# Start message
echo "Annotate with bcftools fill-tags plugin"
date
echo ""

# Load required module
module load bcftools/1.9
echo ""

# Base folder
base_folder="/re_gecip/inherited_cancer_predisposition/ccharalambous/Final_execution"

# Source VCF
source_folder="${base_folder}/data/S09_add_var_id/results"
source_vcf="${source_folder}/annotated_VEP_ClinVar_VarID.vcf.gz"

# Output VCF
output_folder="${base_folder}/data/S10_fill_tags"
output_vcf="${output_folder}/s01_VEP_clinvar_VarID_tags.vcf.gz"
log="${output_folder}/fill_tags.log"

# Progress report
echo "--- Input and output files ---"
echo ""
#echo "original_vcf: ${original_vcf}"
echo "source_vcf: ${source_vcf}"
echo "output_vcf: ${output_vcf}"
echo ""
echo "Working..."

# Index sourse VCF (if not indexed earlier)
#bcftools index "${source_vcf}"

# Annotate with bcftools fill-tags plugin
bcftools plugin \
  fill-tags \
  --output "${output_vcf}" \
  --output-type z \
  "${source_vcf}" \
  &> "${log}"

# Index sourse VCF
bcftools index "${output_vcf}"

# Completion message
echo ""
echo "Done"
date
echo ""
