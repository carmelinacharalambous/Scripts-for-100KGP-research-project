#!/bin/bash

# Sort VCF file
# Started: AL11Dec2018
# Last updated: CC5May2019

# Description:
# Requires input of 1 file
# - concatenated VCF file
# Outputs sorted VCF file

# Reference(s):
# https://samtools.github.io/bcftools/bcftools.html

# Start message
echo "Sort Concatenated VCF using bcftools-sort"
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
source_file="${source_folder}/total_chr.vcf.gz"

# Check number of samples 

x=$(bcftools view -h ${source_file} | tail -n 1 )
echo "Number of samples in the vcf:"
wc -w <<< $x
echo ""

# output_file
output_folder="/re_gecip/inherited_cancer_predisposition/ccharalambous/Final_execution/data/S04_concat_sort_index"
output_file="${output_folder}/sorted.vcf.gz"
 
 
 # Sort VCF
bcftools sort  \
  "${source_file}" \
  --output-type z \
  --output-file "${output_file}"
   
# Progress report
echo ""
echo "Sorted VCF"
date

# Compile source file name
working_folder="/re_gecip/inherited_cancer_predisposition/ccharalambous/Final_execution/scripts/S04_concat_sort_index"
vcf_file="${working_folder}/sorted.vcf.gz"

# Make index VCF
# By defauls bcftools makes "csi" index
# Use -t option create tbi index ("tabix" index)
# Also, remember that you may also use --threads INT option to use multithreading (can make things faster)

echo "Making index ..."

bcftools index -t "${vcf_file}" 
   
# Progress report
echo "Done"
date
