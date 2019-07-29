#!/bin/bash

# check_vcf_samples.sh
# Check list of samples in VCF using zgrep

# Started: 21Jan2019
# Last updated: 11April2019

# Use:
# check_vcf_samples.sh > check_vcf_samples.log
# !!! This runs the job directly on the login node !!!
# This is only acceptable for very small test jobs (e.g. jobs that take 1-2 min, like the one we do in this file)
# In future we will use compute nodes to run any larger jobs (this is done using bsub submission)

# Description:
# Requires input of 2 files:
# - source VCF file
# - file with samples names (optional)  
# Outputs file with intersects of the samples

# Reference(s):
# Skype VCF file specifications (for the #CHROM line)

# Start message
echo "Check list of samples in VCF"
date
echo ""

# Stop at runtime errors
set -e

# BCFtools
module load bcftools/1.9

# Source VCF
data_folder="/gel_data_resources/main_programme/aggregated_illumina_gvcf/GRCH38/20190228/data"

source_file_name="60k_GRCH38_germline_mergedgVCF_chr14_16089563_16105376.bcf"
source_file="${data_folder}/${source_file_name}"

# User files with genes and samples, file name for output 
working_folder="/re_gecip/inherited_cancer_predisposition/ccharalambous/Final_execution/scripts/S02_select_cases_controls"
intersect_checks_folder="${working_folder}/intersects_check"

mkdir "${intersect_checks_folder}"

samples_file_name="final_cases_control_names.txt"
samples_file="${working_folder}/${samples_file_name}"

output_file_name="intersect_labkey_VCF_samples.txt"
output_file="${intersect_checks_folder}/${output_file_name}"

# Progress report
echo "source_file: ${source_file}"
echo "samples_file: ${samples_file}"
echo "output_file: ${output_file}"
echo ""

# --- Count samples in VCF --- 

# Write the VCF samples to a file ( a line)
bcftools view -h "${source_file}" | tail -n 1 > "${intersect_checks_folder}/vcf_samples_line.txt"

# Convert the line to column (replace blanks with new lines, google sed and regex)
sed 's/\t/\n/g' "${intersect_checks_folder}/vcf_samples_line.txt" > "${intersect_checks_folder}/vcf_samples_column.txt"

# Remove the first 9  lines (not samples)
awk 'NR > 9' "${intersect_checks_folder}/vcf_samples_column.txt" > "${intersect_checks_folder}/vcf_samples.txt"

# Sort vcf samples
sort "${intersect_checks_folder}/vcf_samples.txt" > "${intersect_checks_folder}/sorted_vcf_samples.txt"

echo "Number of samples in VCF file:"
cat "${intersect_checks_folder}/sorted_vcf_samples.txt" | wc -l 
echo ""

# Count samples in the Samples file
echo "Number of samples in the Samples file:"
cat "${samples_file}" | wc -l
echo ""

# Sort samples from the samples file
sort "${samples_file}" > "${intersect_checks_folder}/sorted_selected_samples.txt"

# --- Get the overlap using comm (look at comm --help !) --- #

comm -12 "${intersect_checks_folder}/sorted_vcf_samples.txt" "${intersect_checks_folder}/sorted_selected_samples.txt" > "${output_file}"

echo "Intersect samples:"
# cat "${output_file}"
# echo ""
echo "Number of samples in the intersect:"
cat "${output_file}" | wc -l
echo ""

# --- Samples unique for the selected list --- #

comm -13 "${intersect_checks_folder}/sorted_vcf_samples.txt" "${intersect_checks_folder}/sorted_selected_samples.txt" > "${intersect_checks_folder}/unique_selected.txt"

echo "Samples unique for the Samples file:"
# cat unique_selected.txt
# echo ""
echo "Number of samples unique for the Samples file:"
cat "${intersect_checks_folder}/unique_selected.txt" | wc -l
echo ""

# --- Samples unique for the VCF file --- #
comm -23 "${intersect_checks_folder}/sorted_vcf_samples.txt" "${intersect_checks_folder}/sorted_selected_samples.txt" > "${intersect_checks_folder}/unique_vcf.txt"

#echo "Samples unique for the VCF file:"
#cat unique_vcf.txt
#echo ""
echo "Number of samples unique for the VCF file:"
cat "${intersect_checks_folder}/unique_vcf.txt" | wc -l
echo ""

# Clean-up
# rm -fr "${intersect_checks_folder}" 

# Completion message
echo "Done"
date
