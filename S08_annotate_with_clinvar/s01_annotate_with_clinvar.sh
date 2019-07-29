#!/bin/bash

# Annotate with ClinVar

# Started: CC25Mar2019
# Last updated: CC5May2019

# Stop at runtime errors
set -e

# Start message
echo "Annotate with ClinVar"
date
echo ""

# Load required module
module load bcftools/1.9
echo ""

# Base folder
base_folder="/re_gecip/inherited_cancer_predisposition/ccharalambous/Final_execution"

# Source VCF
data_folder="${base_folder}/data/S07_annotate_with_VEP/results"
source_vcf="${data_folder}/s01_vep.vcf.gz"

# Output VC
results_folder="${base_folder}/data/S08_annotate_with_clinvar"
output_vcf="${results_folder}/s01_VEP_clinvar.vcf.gz"
log="${results_folder}/clinvar_annotation.log"


# ClinVar data
clinvar_data_folder="${base_folder}/scripts/S08_annotate_with_clinvar/clinvar_update"
clinvar_vcf="${clinvar_data_folder}/clinvar_20190219_updated.vcf.gz"

# Progress report
echo "--- Input and output files ---"
echo ""
echo "source_vcf: ${source_vcf}"
echo "output_vcf: ${output_vcf}"
echo "clinvar_vcf: ${clinvar_vcf}"
echo ""
echo "Working..."

# Index the sourse file (because this was not done by VEP, in the previous step)
bcftools index "${source_vcf}"

# Annotate with bcftools
# http://www.htslib.org/doc/bcftools.html#annotate
# Edits ID field (which was empty in the initial file) and adds INFO 
# fileds (because there was no fields with such name in the initial file) 
bcftools annotate \
  --annotations "${clinvar_vcf}" \
  --columns INFO,ID \
  --output "${output_vcf}" \
  --output-type z \
  "${source_vcf}" \
  &> "${log}"

# Completion message
echo ""
echo "Done"
date
echo ""

# --- Description of the added ClinVar annotations --- #

##ID=<Description="ClinVar Variation ID">

##INFO=<ID=ALLELEID,Number=1,Type=Integer,Description="the ClinVar Allele ID">
##INFO=<ID=CLNREVSTAT,Number=.,Type=String,Description="ClinVar review status for the Variation ID">

##INFO=<ID=CLNSIG,Number=.,Type=String,Description="Clinical significance for this single variant">
##INFO=<ID=CLNDN,Number=.,Type=String,Description="ClinVar's preferred disease name for the concept specified by disease identifiers in CLNDISDB">
##INFO=<ID=CLNDISDB,Number=.,Type=String,Description="Tag-value pairs of disease database name and identifier, e.g. OMIM:NNNNNN">

##INFO=<ID=CLNSIGCONF,Number=.,Type=String,Description="Conflicting clinical significance for this single variant">
##INFO=<ID=SSR,Number=1,Type=Integer,Description="Variant Suspect Reason Codes. One or more of the following values may be added: 0 - unspecified, 1 - Paralog, 2 - byEST, 4 - oldAlign, 8 - Para_EST, 16 - 1kg_failed, 1024 - other">

##INFO=<ID=CLNVI,Number=.,Type=String,Description="the variant's clinical sources reported as tag-value pairs of database and variant identifier">
##INFO=<ID=ORIGIN,Number=.,Type=String,Description="Allele origin. One or more of the following values may be added: 0 - unknown; 1 - germline; 2 - somatic; 4 - inherited; 8 - paternal; 16 - maternal; 32 - de-novo; 64 - biparental; 128 - uniparental; 256 - not-tested; 512 - tested-inconclusive; 1073741824 - other">

##INFO=<ID=GENEINFO,Number=1,Type=String,Description="Gene(s) for the variant reported as gene symbol:gene id. The gene symbol and id are delimited by a colon (:) and each pair is delimited by a vertical bar (|)">
##INFO=<ID=CLNVCSO,Number=1,Type=String,Description="Sequence Ontology id for variant type">
##INFO=<ID=CLNVC,Number=1,Type=String,Description="Variant type">
##INFO=<ID=MC,Number=.,Type=String,Description="comma separated list of molecular consequence in the form of Sequence Ontology ID|molecular_consequence">

##INFO=<ID=RS,Number=.,Type=String,Description="dbSNP ID (i.e. rs number)">
##INFO=<ID=CLNHGVS,Number=.,Type=String,Description="Top-level (primary assembly, alt, or patch) HGVS expression.">
##INFO=<ID=DBVARID,Number=.,Type=String,Description="nsv accessions from dbVar for the variant">

##INFO=<ID=AF_ESP,Number=1,Type=Float,Description="allele frequencies from GO-ESP">
##INFO=<ID=AF_EXAC,Number=1,Type=Float,Description="allele frequencies from ExAC">
##INFO=<ID=AF_TGP,Number=1,Type=Float,Description="allele frequencies from TGP">

##INFO=<ID=CLNSIGINCL,Number=.,Type=String,Description="Clinical significance for a haplotype or genotype that includes this variant. Reported as pairs of VariationID:clinical significance.">
##INFO=<ID=CLNDNINCL,Number=.,Type=String,Description="For included Variant : ClinVar's preferred disease name for the concept specified by disease identifiers in CLNDISDB">
##INFO=<ID=CLNDISDBINCL,Number=.,Type=String,Description="For included Variant: Tag-value pairs of disease database name and identifier, e.g. OMIM:NNNNNN">
