#!/bin/bash
# s01_rsync.sh

# Copy pre-selected overlapping sub-sets of 1000-genomes and GEL projects:
# 100,000 randomly selected biallelic SNPs that have AF>0.05 and call rate >85% in each of the datasets
# Also, copy list of 12,449 samples selected for analysis (add phenotypes later!)

# Started: AL27Apr2019
# Last updated: CC20May2019

# Intended use: RUN ON THE DESKTOP TERMINAL
# s01_rsync.sh > s01_rsync.log

# Notes: 
# Its better to use rsync than cp to copy large files because rsync checks the files after copying 
# Its also a good practice to accompany large files with md5 sums - to have a way to verify the files integrity
# Correct the base folder if run on HPC

# Stop at runtime errors
set -e

# Start message
echo "Copy pre-selected sub-sets of 1000-genomes and GEL data"
date
echo ""

# Files and folders
base_folder="/home/ccharalambous/re_gecip/inherited_cancer_predisposition/ccharalambous"
target_folder="${base_folder}/Final_execution/data/S13_joined_gel_1kg_PCA/s01_copy_source_data"

# Subset of 1000-genomes data
files="/home/alarionov/re_gecip/inherited_cancer_predisposition/tischkowitz/users/alexey/template_scripts/w040_prepare_1kg_genotypes_for_PCA/data/all_1kg_100kv_b38_genotypes_filtered."
rsync -avh "${files}"* "${target_folder}/"
echo ""

# Subset of GEL data
files="/home/alarionov/re_gecip/inherited_cancer_predisposition/tischkowitz/users/alexey/template_scripts/w050_prepare_gel_genotypes_for_PCA/data/all_gel_100kv_b38_genotypes."
rsync -avh "${files}"* "${target_folder}/"
echo ""

# List of selected GEL samples (of April 2019)
samples_file="${base_folder}/Final_execution/scripts/S02_select_cases_controls/intersects_check/intersect_labkey_VCF_samples.txt"
rsync -avh "${samples_file}" "${target_folder}/"
echo ""

# Phenotypes for most of the selected samples (of Mayt 2019)
phenotype_file="${base_folder}/Final_execution/data/S02_select_cases_controls/phenotypes.txt"
rsync -avh "${phenotype_file}" "${target_folder}/"
echo ""

# Completion message
echo "Done"
date
echo ""
