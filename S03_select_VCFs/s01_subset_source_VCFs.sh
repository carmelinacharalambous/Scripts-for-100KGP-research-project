#!/bin/bash

# s02_subset_source_VCFs.sh

# Started: CC18FebJan2019
# Last updated: CC14Mar2019

# Start message
echo "Extract data for required genes and samples from source VCFs"
date
echo ""

# Stop at runtime errors
set -e

# Load bcftools module
module load bcftools/1.9

# Files and folders
base_folder="/re_gecip/inherited_cancer_predisposition/ccharalambous/Final_execution"

source_folder="/gel_data_resources/main_programme/aggregated_illumina_gvcf/GRCH38/20190228/data"
output_folder="${base_folder}/data/S03_select_VCFs"

genes_file="${base_folder}/scripts/S01_select_genes/gene_list_gVCF.bed"
samples_file="${base_folder}/scripts/S02_select_cases_controls/intersects_check/intersect_labkey_VCF_samples.txt"

# Progress report
echo "genes_file: ${genes_file}"
echo "samples_file: ${samples_file}"
echo ""
echo "Started extracting data ..."
echo ""

# Make the output folder
rm -fr "${output_folder}" # remove if existed
mkdir -p "${output_folder}"

# select VCFs
while  read chr start end gene file_suffix
do

  # Compile setings and files
  region="${chr}:${start}-${end}"
  source_file="${source_folder}/60k_GRCH38_germline_mergedgVCF_${file_suffix}.bcf"
  output_file1="${output_folder}/${gene}_raw.bcf"
  log_file1="${output_folder}/${gene}_raw.log"

# Check existance of sourse file
  if [ ! -e "${source_file}" ] 
  then 
    
    # Error message
    echo "----- ${gene} -----"
    echo ""
    echo "The sourse file does not exist!"
    echo ""
    
    # Next gene
    continue
    
  fi

  # Check existance of sourse file index
  if [ ! -e "${source_file}.csi" ] 
  then 
    
    # Error message
    echo "----- ${gene} -----"
    echo ""
    echo "The index for sourse file does not exist!"
    echo ""
    
    # Next gene
    continue
    
  fi

  #Select samples and cases
   bcftools view \
    "${source_file}" \
    --samples-file "${samples_file}" \
    --force-samples \
    --regions "${region}" \
    --output-file "${output_file1}" \
    --output-type b \
    &> "${log_file1}"
    
  # Indexing
  bcftools index "${output_file1}"

  # Compile setings and files
   output_file2="${output_folder}/${gene}_filtered.bcf"
   log_file2="${output_folder}/${gene}_filtered.log"

  # Extract data from source file
  bcftools view \
    "${output_file1}" \
    --exclude 'AC=0' \
    --apply-filters 'PASS' \
    --output-file "${output_file2}" \
    --output-type b \
    &> "${log_file2}"
    
    # Indexing
    bcftools index "${output_file2}"

  # Progress report
  echo "----- ${gene} -----"
  echo ""
  echo "region: ${region}"
  echo "source_file: ${source_file}"
    echo "output_file1: ${output_file1}"
  echo "output_file2: ${output_file2}"
  echo "log_file1: ${log_file1}"
  echo "log_file2: ${log_file2}"
  echo ""
  bcftools view "${output_file2}" | tail | awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}'
  echo ""
  
  # Remove files -
   #rm "${output_file1}"  "${output_file1}.csi"
  
done < "${genes_file}" # Next gene

# Completion message
echo "Done all genes"
date
