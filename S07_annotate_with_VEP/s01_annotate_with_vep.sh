#!/bin/bash
#CarmelinaCharalambous_5_May_2019
# Annotate with VEP

# Note:
# cache version 93 was recommended by the GEL example script, 
# despite the fact that VEP module version is 92. 
# Just in case, we may ask later: why is it so?

# Stop at runtime errors
set -e

# Start message
echo "Annotate with VEP"
date
echo ""

# Load required module
module load vep/92
module show vep/92
echo ""

# Base folder
base_folder="/re_gecip/inherited_cancer_predisposition/ccharalambous/Final_execution"

# Source VCF
data_folder="${base_folder}/data/S06_filter_by_qual"
source_vcf="${data_folder}/filtered_by_qual.vcf.gz"

# Files and folders
working_folder="${base_folder}/data/S07_annotate_with_VEP"
results_folder="${working_folder}/results"

output_vcf="${results_folder}/s01_vep.vcf.gz"
vep_log="${results_folder}/vep_annotation.log"
vep_report="${results_folder}/vep_report.html"

# Make output folder, if it does not exist
if [ !  -d "${results_folder}" ]
then
  mkdir "${results_folder}"
fi

# Cache folder
cache_folder="/tools/apps/vep/92/ensembl-vep/.vep"

# Plugins folder
plugins_folder="${cache_folder}/Plugins"

# CADD plugin data
cadd_data_folder="/public_data_resources/CADD/v1.4/GRCh38"
cadd_snv_data="${cadd_data_folder}/whole_genome_SNVs.tsv.gz"
cadd_indels_data="${cadd_data_folder}/InDels.tsv.gz"

# Progress report
echo "--- Input and output files ---"
echo ""
echo "source_vcf: ${source_vcf}"
echo "output_vcf: ${output_vcf}"
echo "vep_log: ${vep_log}"
echo "vep_report: ${vep_report}"
echo ""
echo "--- Data sources ---"
echo ""
echo "See used VEP cache description in the following file:"
echo "${vep_cache_info}"
echo ""
echo "Used CADD annotation files:"
echo "${cadd_snv_data}"
echo "${cadd_indels_data}"
echo ""
echo "Working..."

# Annotate VCF
vep \
  --input_file "${source_vcf}" \
  --output_file "${output_vcf}" \
  --vcf \
  --force_overwrite \
  --compress_output bgzip \
  --stats_file "${vep_report}" \
  --offline \
  --cache \
  --dir_cache "${cache_folder}" \
  --cache_version 93 \
  --assembly GRCh38 \
  --everything \
  --total_length \
  --check_existing \
  --exclude_null_alleles \
  --pick \
  --gencode_basic \
  --dir_plugins "${plugins_folder}" \
  --plugin CADD,"${cadd_snv_data}","{$cadd_indels_data}" \
  &> "${vep_log}"

#  --nearest symbol \
#  --check_ref \

# Completion message
echo "Done"
date
