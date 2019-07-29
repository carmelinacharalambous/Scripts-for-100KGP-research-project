#!/bin/bash

# Select 8,628 GEL samples reported white british/irish  

# The followig variants are removed after samples selection: 
# - MAF < 0.05
# - call rates < 0.85

# Similar filters had been applied previously to the initial GEL sub-set
# However, MAFs and call rates might slightly change after samples selection

# Started: CC29Apr2019
# Last updated: CC3Jun2019

# Stop at runtime errors
set -e

# Start message
echo "Select 8,628 GEL samples reported white british/irish"
date
echo ""

# Load required module
module load bcftools/1.9

# Files and folders
base_folder="/re_gecip/inherited_cancer_predisposition/ccharalambous/Final_execution"
scripts_folder="${base_folder}/scripts/S14_gel_only_PCA/s01_select_gel_samples"

source_vcf="${base_folder}/data/S13_add_eigenvectors/s01_copy_source_data/all_gel_100kv_b38_genotypes.vcf.gz"
samples_file="${base_folder}/data/S13_add_eigenvectors/s05_make_PCA_plots/reported_white_british_irish.txt"

output_data_folder="${base_folder}/data/S14_gel_only_PCA/s01_select_gel_samples"
interim_vcf="${output_data_folder}/selected_gel_cases_100kv_b38_genotypes.vcf.gz"
output_vcf="${output_data_folder}/selected_gel_cases_filtered_b38_genotypes.vcf.gz"

# Load additional functions to count number of samples and variants in a VCF file
source "${base_folder}/scripts/S14_gel_only_PCA/f01_accessory_functions/f01_functions.sh"
echo "Loaded accessory functions"
echo ""

# Make output data folder
rm -fr "${output_data_folder}" # remove if existed
mkdir -p "${output_data_folder}"

# Progress report
echo "samples_file: ${samples_file}"
echo "source_vcf: ${source_vcf}"
echo "output_vcf: ${output_vcf}"
echo ""

num_of_samples=$(cat ${samples_file} | wc -l) # count samples in the samples file
num_of_samples=$(printf "%'d" "${num_of_samples}") # add comma after thousands
echo "Number of samples in the samples file: ${num_of_samples}"
echo ""

# Select samples
echo "Selecting samples ..."

bcftools view \
  "${source_vcf}" \
  --samples-file "${samples_file}" \
  --output-file "${interim_vcf}" \
  --output-type z

date
echo ""

# Filter variants after the samples selection:
echo "Filtering variants ..."

bcftools view \
  "${interim_vcf}" \
  --min-af 0.05 \
  --max-af 0.95 \
  --exclude "F_MISSING > 0.15" \
  --output-file "${output_vcf}" \
  --output-type z

date
echo ""

# Index
echo "Indexing ..."

bcftools index "${output_vcf}"

date
echo ""

# Progress report

echo "Number of variants and samples in the input and output VCF files:"
echo Source: $(num_of_variants "${source_vcf}") x $(num_of_samples "${source_vcf}")
echo Output: $(num_of_variants "${output_vcf}") x $(num_of_samples "${output_vcf}")
echo ""

# Clean-up
rm "${interim_vcf}"

# Completion message
echo "Done"
date
echo ""
