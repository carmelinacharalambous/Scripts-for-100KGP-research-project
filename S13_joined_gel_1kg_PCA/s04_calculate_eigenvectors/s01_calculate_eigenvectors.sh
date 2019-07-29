#!/bin/bash

# Started: CC30Apr2019
# Last updated: CC21May2019

# Stop at runtime errors
set -e

# Start message
echo "Calculate eigenvectors for joined dataset of gel and 1000-genomes"
date
echo ""

# Load PLINK module
module load PLINK/1.9b_4.6-x86_64

# Folders
base_folder="//re_gecip/inherited_cancer_predisposition/ccharalambous/Final_execution"
scripts_folder="${base_folder}/scripts/S13_joined_gel_1kg_PCA/s04_calculate_eigenvectors"

###################################################
#                                          Import VCF to PLINK                                               #
###################################################

# Source VCF
source_vcf_folder="${base_folder}/data/s13_joined_gel_1kg_PCA/s03_merge_gel_and_1kg"
source_vcf_file="${source_vcf_folder}/gel_1kg_b38_genotypes.vcf.gz"

# Outputted initial PLINK file set
plink_dataset_folder="${base_folder}/data/S13_joined_gel_1kg_PCA/s04_calculate_eigenvectors/p01_vcf_to_plink"
rm -fr "${plink_dataset_folder}"
mkdir -p "${plink_dataset_folder}"
initial_plink_dataset="${plink_dataset_folder}/gel_1kg_99kv_15ks_b38"

# --- Import VCF to PLINK bed-fam-bim --- #

# --vcf-half-call describes what to do with genotypes like 0/.
# --allow-no-sex suppresses warning about missed sex data
# --double-id puts sample name to both Family-ID and Participant-ID
# --silent suppresses very verbouse output to the "out" file (log file is still available in the data folder)

plink \
  --vcf "${source_vcf_file}" \
  --vcf-half-call "missing" \
  --double-id \
  --allow-no-sex \
  --make-bed \
  --silent \
  --out "${initial_plink_dataset}"

# Progress report
echo "Imported joined gel-1kg VCF to PLINK"
echo "source_vcf_file: ${source_vcf_file}"
echo "initial_plink_dataset: ${initial_plink_dataset}"
date
echo ""

###################################################
#                                           Exclude variants in LD                                            #
###################################################

# Output folder
output_data_folder="${base_folder}/data/S13_joined_gel_1kg_PCA/s04_calculate_eigenvectors/p02_exclude_variants_in_LD"
rm -fr "${output_data_folder}"
mkdir -p "${output_data_folder}"

# Output files
pairphase_LD="${output_data_folder}/pairphase_LD"
LD_pruned_dataset="${output_data_folder}/gel_1kg_70kv_not_in_LD"

# --- Determine variants in LD --- #

# Command indep-pairphase makes two files:
# - list of variants in LD (file with extension .prune.out) 
# - list of variants not in LD ( extension .prune.in)
# The specific parameters 50 5 0.5 are taken from an example 
# discussed in PLINK 1.07 manual for LD pruning

plink \
  --bfile "${initial_plink_dataset}" \
  --indep-pairphase 50 5 0.5 \
  --allow-no-sex \
  --silent \
  --out "${pairphase_LD}"

# --- Make a new bed-fam-bim PLINK dataset excluding variants in LD --- #

# Use list of variants to exclude from the file .prune.out made in the previous step

plink \
  --bfile "${initial_plink_dataset}" \
  --exclude "${pairphase_LD}.prune.out" \
  --allow-no-sex \
  --make-bed \
  --silent \
  --out "${LD_pruned_dataset}"

# Progress report
echo "Excluded variants in pair-wise LD"
echo "LD_pruned_dataset: ${LD_pruned_dataset}"
date
echo ""

###################################################
#                  Thinn dataset to 10k random variants and calculate PCs                  #
###################################################

# Output files and folder
output_data_folder="${base_folder}/data/S13_joined_gel_1kg_PCA/s04_calculate_eigenvectors/p03_thinn_and_pca"
rm -fr "${output_data_folder}"
mkdir -p "${output_data_folder}"
thinned_dataset="${output_data_folder}/gel_1kg_10kv"
output_PCs="${output_data_folder}/gel_1kg_10kv_100PCs"

# --- Select random 10k variants --- #

plink \
  --bfile "${LD_pruned_dataset}" \
  --thin-count 10000 \
  --make-bed \
  --silent \
  --out "${thinned_dataset}"

# Progress report
echo "Selected 10,000 random variants"
echo "thinned_dataset: ${thinned_dataset}"
date
echo ""

# --- Calculate 100 top PC-s --- #

# "header" and "tabs" options are used just to format output

plink \
  --bfile "${thinned_dataset}" \
  --pca 100 header tabs \
  --silent \
  --out "${output_PCs}"

# Progress report
echo "Calculated 100 top PCs based on 10k random common variants not in LD"
echo "output_PCs: ${output_PCs}"
echo ""

# Completion message
echo "Done all tasks"
date
echo ""
