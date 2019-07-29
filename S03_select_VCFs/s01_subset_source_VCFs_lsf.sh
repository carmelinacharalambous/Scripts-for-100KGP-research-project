#!/bin/bash

script_folder="/re_gecip/inherited_cancer_predisposition/ccharalambous/Final_execution/scripts/S03_select_VCFs"

script="${script_folder}/s01_subset_source_VCFs.sh"

output_file="${script_folder}/s01_subset_vcf.out"
error_file="${script_folder}/s01_subset_vcd.err"

bsub -q gecip -P re_gecip_inherited_cancer_predisposition -o "${output_file}" -e "${error_file}" "${script}" 