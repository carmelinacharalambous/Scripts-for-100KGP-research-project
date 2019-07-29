#!/bin/bash

# Sort VCF file
# Started: CC5May2019
# Last updated: CC5May2019

# Description:
# Requires input of 1 file
# - sorted VCF file
# Outputs filtered VCF file

# Reference(s):
# https://samtools.github.io/bcftools/bcftools.html

# Start message
echo "Filter sorted VCF file"
date
echo ""

# Stop at runtime errors
set -e

# Load bcftools module
module load bcftools/1.9
bcftools --version
echo ""

# Source VCFs
source_folder="/re_gecip/inherited_cancer_predisposition/ccharalambous/Final_execution/data/S04_concat_sort_index"
source_file="${source_folder}/sorted.vcf.gz"

# output_file
output_folder="/re_gecip/inherited_cancer_predisposition/ccharalambous/Final_execution/data/S06_filter_by_qual"
output_file="${output_folder}/filtered_by_qual.vcf.gz"
 
 # Sort VCF
bcftools view  \
  "${source_file}" \
  --include "QUAL>=549" \
  --output-type z \
  --output-file "${output_file}"
   
# Check result 
echo "Num of var in the source file:"
zgrep -v ^# "${source_file}" | wc -l

echo "Num of var in the output file:"
zgrep -v ^# "${output_file}" | wc -l

# Progress report
echo ""
echo "Filtered VCF"
date
   