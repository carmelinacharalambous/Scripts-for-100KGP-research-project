#!/bin/bash

# s01_merge_gel_and_1kg

# The script performs 3 consequitive steps:
# 1) Select 1kg variants present in gel vcf
# 2) Merge gel and 1kg vcfs
# 3) Remove multiallelic variants

# It might be possible to intersect, merge and filter in one step, just using additional options in "merge". 
# However, because the bcftools documentation is quite limited, it might be safer to do it step-by-step for now.

# Started: CC30Apr2019
# Last updated: CC30Apr2019

# ----------------------------------------------------------------------------------------------------------- #
#                                         Set environment                                           #
# ----------------------------------------------------------------------------------------------------------- #

# Stop at runtime errors
set -e

# Start message
echo "Merge gel_and 1kg vcfs for PCA"
date
echo ""

# Load required module
module load bcftools/1.9
echo ""

# Files and folders
base_folder="//re_gecip/inherited_cancer_predisposition/ccharalambous/Final_execution"
scripts_folder="${base_folder}/scripts/S13_add_eigenvectors/s03_merge_gel_and_1kg"

source_1kg_folder="${base_folder}/data/S13_add_eigenvectors/s01_copy_source_data"
source_1kg_vcf="${source_1kg_folder}/all_1kg_100kv_b38_genotypes_filtered.vcf.gz"

source_gel_folder="${base_folder}/data/S13_add_eigenvectors/s02_select_gel_samples"
source_gel_vcf="${source_gel_folder}/selected_gel_cases_filtered_b38_genotypes.vcf.gz"

output_data_folder="${base_folder}/data/S13_add_eigenvectors/s03_merge_gel_and_1kg"
rm -fr "${output_data_folder}"
mkdir "${output_data_folder}"

intersect_1kg_vcf="${output_data_folder}/intersect_1kg_b38_genotypes.vcf.gz"
merged_vcf_raw="${output_data_folder}/gel_1kg_b38_genotypes_raw.vcf.gz"
output_vcf="${output_data_folder}/gel_1kg_b38_genotypes.vcf.gz"

# Progress report
echo "scripts_folder: ${scripts_folder}"
echo "source_1kg_vcf: ${source_1kg_vcf}"
echo "source_gel_vcf: ${source_gel_vcf}"
echo "output_vcf: ${output_vcf}"
echo ""

# ----------------------------------------------------------------------------------------------------------- #
#                                   Load accessory functions                                    #
# ----------------------------------------------------------------------------------------------------------- #

# Loaded functions:
# num_of_variants   vcf_file_name
# num_of_samples   vcf_file_name

source "${base_folder}/scripts/S13_add_eigenvectors/f01_accessory_functions/f01_functions.sh"

echo "Loaded accessory functions"
echo ""

# ----------------------------------------------------------------------------------------------------------- #
#                                Get 1kg variants present in gel                             #
# ----------------------------------------------------------------------------------------------------------- #

# A simplified isec version is used because we know in advance that 
# each gel variant is present in 1kg vcf

# Select variants
echo "Get 1kg variants present in gel ..."

bcftools isec \
  --output "${intersect_1kg_vcf}" \
  --output-type z \
  --nfiles=2 \
  --write 1 \
  "${source_1kg_vcf}" \
  "${source_gel_vcf}"

# Index
bcftools index "${intersect_1kg_vcf}"

# Progress report
date
echo ""
echo "Number of variants and samples in the input and output VCF files:"
echo Source 1kg: $(num_of_variants "${source_1kg_vcf}") x $(num_of_samples "${source_1kg_vcf}")
echo Intersect 1kg: $(num_of_variants "${intersect_1kg_vcf}") x $(num_of_samples "${intersect_1kg_vcf}")
echo ""

# ----------------------------------------------------------------------------------------------------------- #
#                                         Merge 1kg and gel                                       #
# ----------------------------------------------------------------------------------------------------------- #

# Merge
echo "Merging 1kg and gel VCFs ..."

bcftools merge \
  "${intersect_1kg_vcf}" \
  "${source_gel_vcf}" \
  --output "${merged_vcf_raw}" \
  --output-type z

# Index
bcftools index "${merged_vcf_raw}"

# Progress report
date
echo ""
echo "Number of variants and samples in the input and output files:"
echo Source 1kg: $(num_of_variants "${intersect_1kg_vcf}") x $(num_of_samples "${intersect_1kg_vcf}")
echo Source GEL: $(num_of_variants "${source_gel_vcf}") x $(num_of_samples "${source_gel_vcf}")
echo Merged: $(num_of_variants "${merged_vcf_raw}") x $(num_of_samples "${merged_vcf_raw}")
echo ""

# ----------------------------------------------------------------------------------------------------------- #
#                                          Filter merged file                                        #
# ----------------------------------------------------------------------------------------------------------- #

# This final filtering assures that vcf contains only commom (MAF>0.05) biallelic SNPs with call rate > 0.85
# Because of the previous filtering steps, it is likely that this step will not to remove any variants
# Hhowever, the merge may create multiallelic sites and slightly change MAF and call rates

# Filter
echo "Filtering merged file ..."

bcftools view \
  "${merged_vcf_raw}" \
  --min-af 0.05 \
  --max-af 0.95 \
  --min-alleles 2 \
  --max-alleles 2 \
  --types snps \
  --exclude "F_MISSING > 0.15" \
  --output-file "${output_vcf}" \
  --output-type z

# Index
bcftools index "${output_vcf}"

# Progress report
date
echo ""
echo "Number of variants and samples in the input and output files:"
echo Raw merged file: $(num_of_variants "${merged_vcf_raw}") x $(num_of_samples "${merged_vcf_raw}")
echo Filtered merged file: $(num_of_variants "${output_vcf}") x $(num_of_samples "${output_vcf}")
echo ""

# ----------------------------------------------------------------------------------------------------------- #
#                                               Clean up                                                #
# ----------------------------------------------------------------------------------------------------------- #

rm "${intersect_1kg_vcf}" "${merged_vcf_raw}" "${intersect_1kg_vcf}.csi" "${merged_vcf_raw}.csi"

# Completion message
echo "Done"
date
echo ""
