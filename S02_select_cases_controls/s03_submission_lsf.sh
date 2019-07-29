#!/bin/bash

script_folder="/re_gecip/inherited_cancer_predisposition/ccharalambous/Final_execution/scripts/S02_select_cases_controls"
script="${script_folder}/s03_check_vcf_samples.sh"

output_file="${script_folder}/s03_check_vcf_samples.out"
error_file="${script_folder}/s03_check_vcf_samples.err"

bsub -q gecip -P re_gecip_inherited_cancer_predisposition -o "${output_file}" -e "${error_file}" "${script}" 