#!/bin/bash

script_folder="/re_gecip/inherited_cancer_predisposition/ccharalambous/Final_execution"
script="${script_folder}/scripts/S04_concat_sort_index/s02_sort_filter_VCF.sh"

output_file="${script_folder}/data/S04_concat_sort_index/s02_sorted.out"
error_file="${script_folder}/scripts/S04_concat_sort_index/s02_sorted.err"

bsub -q gecip -P re_gecip_inherited_cancer_predisposition -o "${output_file}" -e "${error_file}" "${script}" 