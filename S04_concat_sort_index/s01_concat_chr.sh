#!/bin/bash

# subset_vcd.sh
# Select sample(s) and gene(s) from multisample VCF using bcftools

# Started: AL11Dec2018
# Last updated: CC5May2019

# Use:
# subset_vcd.sh > subset_vcd.log
# !!! This runs the job directly on the login node !!!
# This is only acceptable for very small test jobs (e.g. jobs that take 1-2 min, like the one we do in this file)
# In future we will use compute nodes to run any larger jobs (this is done using bsub submission)

# Description:
# Requires input of 3 files:
# - source VCF file
# - file with genes coordinates
# - file with samples names
# Outputs VCF file

# Reference(s):
# https://samtools.github.io/bcftools/bcftools.html

# Start message
echo "Concatenate VCFs using bcftools"
date
echo ""

# Stop at runtime errors
set -e

# Load bcftools module
module load bcftools/1.9
bcftools --version
echo ""

# Source VCFs
source_folder="/re_gecip/inherited_cancer_predisposition/ccharalambous/Final_execution/data/S03_select_VCFs"
source_list_files="${source_folder}/chrm_list_names.txt"

# output_file
output_folder="/re_gecip/inherited_cancer_predisposition/ccharalambous/Final_execution/data/S04_concat_sort_index"
output_file="${output_folder}/total_chr.vcf.gz"
 
 # Concat VCF
bcftools concat \
    --file-list "${source_list_files}" \
    --output-type z \
    --threads 2 \
     > "${output_file}"
   

# Progress report
echo ""
echo "Concatenated"
date

