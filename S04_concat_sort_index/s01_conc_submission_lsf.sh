#!/bin/bash

script_folder="/re_gecip/inherited_cancer_predisposition/ccharalambous/Final_execution"
script="${script_folder}/scripts/S04_concat_sort_index/s01_concat_chr.sh"

output_file="${script_folder}/data/S04_concat_sort_index/s01_concat_chr.out"
error_file="${script_folder}/data/S04_concat_sort_index/s01_concat_chr.err"

bsub -q gecip -P re_gecip_inherited_cancer_predisposition -o "${output_file}" -e "${error_file}" "${script}" 