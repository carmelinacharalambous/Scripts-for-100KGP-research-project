# s04_filter_by_functional_consequences.R
# last updated: CC24May2019  

# Run on HPC R-3.5.1 requesting 1 core and 30GB RAM for up to 1hrs

# Keep only variants that  
# - belong to either ClinVar, or LoF or FIM (as selected earlier)  
# - are of low-frequency (MAF < 0.05)  

# Input data: 73,907 variants x 12,449 cases
# Output data: 477 variants x 12,449 cases

#####################################################################
#                           Start section                           #
#####################################################################

# Start messages
cat("Filter by functional annotations \n\n")
Sys.time()
version$version.string
gc()
cat("\n")

# Clean-up
rm(list=ls())
graphics.off()

# Base folder
base_folder <- "/re_gecip/inherited_cancer_predisposition/ccharalambous/Final_execution"

# Additional options
options(stringsAsFactors = F,
        warnPartialMatchArgs = T, 
        warnPartialMatchAttr = T, 
        warnPartialMatchDollar = T)

# Progress report
cat("Completed start section \n")
Sys.time()
cat("\n")

#####################################################################
#                          Load source data                         #
#####################################################################

load(paste(base_folder, "data/S12_explore_and_filter/s02_fixed_df.RData", sep="/"))
load(paste(base_folder, "data/S11_VCF_to_R/s01_gt_add_mx.RData", sep="/"))
load(paste(base_folder, "data/S11_VCF_to_R/s01_gq_mx.RData", sep="/"))
load(paste(base_folder, "data/S11_VCF_to_R/s01_dp_mx.RData", sep="/"))
load(paste(base_folder, "data/S11_VCF_to_R/s01_ref_mx.RData", sep="/"))
load(paste(base_folder, "data/S11_VCF_to_R/s01_alt_mx.RData", sep="/"))

cat("Loaded source data \n")
Sys.time()
cat("Size of the loaded objects: \n")
cat("fixed.df:",dim(fixed.df),"\n")
cat("gt_add.mx:",dim(gt_add.mx),"\n")
cat("gq.mx:",dim(gq.mx),"\n")
cat("dp.mx:",dim(dp.mx),"\n")
cat("ref.mx:",dim(ref.mx),"\n")
cat("alt.mx:",dim(alt.mx),"\n\n")

cat("Check that data is in sync: \n")
sum(rownames(fixed.df) != rownames(gt_add.mx))
sum(rownames(fixed.df) != rownames(gq.mx))
sum(rownames(fixed.df) != rownames(dp.mx))
sum(rownames(fixed.df) != rownames(ref.mx))
sum(rownames(fixed.df) != rownames(alt.mx))

sum(colnames(gt_add.mx) != colnames(gq.mx))
sum(colnames(gt_add.mx) != colnames(dp.mx))
sum(colnames(gt_add.mx) != colnames(ref.mx))
sum(colnames(gt_add.mx) != colnames(alt.mx))
cat("\n")

#####################################################################
#                           Select variants                         #
#####################################################################


# Functional annotation
selected_variants <- fixed.df$ClinVar | fixed.df$LoF | fixed.df$FIM

# Allelic frequency
selected_variants <- selected_variants & fixed.df$MAF <= 0.05

cat("Number of selected variants:\n")
sum(selected_variants)
cat("\n")

# Update matrices
cat("Update the data \n")
fixed.df <- fixed.df[selected_variants,]
gt_add.mx <- gt_add.mx[selected_variants,]
gq.mx <- gq.mx[selected_variants,]
dp.mx <- dp.mx[selected_variants,]
ref.mx <- ref.mx[selected_variants,]
alt.mx <- alt.mx[selected_variants,]

cat("Size of the updated objects: \n")
cat("fixed.df:",dim(fixed.df),"\n")
cat("gt_add.mx:",dim(gt_add.mx),"\n")
cat("gq.mx:",dim(gq.mx),"\n")
cat("dp.mx:",dim(dp.mx),"\n")
cat("ref.mx:",dim(ref.mx),"\n")
cat("alt.mx:",dim(alt.mx),"\n\n")

cat("Check that the updated matrices are still in sync: \n")
sum(rownames(fixed.df) != rownames(gt_add.mx))
sum(rownames(fixed.df) != rownames(gq.mx))
sum(rownames(fixed.df) != rownames(dp.mx))
sum(rownames(fixed.df) != rownames(ref.mx))
sum(rownames(fixed.df) != rownames(alt.mx))

sum(colnames(gt_add.mx) != colnames(gq.mx))
sum(colnames(gt_add.mx) != colnames(dp.mx))
sum(colnames(gt_add.mx) != colnames(ref.mx))
sum(colnames(gt_add.mx) != colnames(alt.mx))
cat("\n")

rm(selected_variants)

#####################################################################
#           Remove unnecessary annotations from fixed.df            #
#####################################################################

cat("Remove unnecessary annotations from fixed.df \n\n")

selected_annotations <- c("VarID","CHROM","POS","REF","ALT","STRAND","Existing_variation","VARIANT_CLASS","BIOTYPE",
                          "QUAL","DP","meanDP","meanGQ","Initial_AN","Initial_AC","Initial_AF","MAF","HWE","possibly_multiallelic",
                          "SYMBOL","Allele","IMPACT","Consequence","EXON","INTRON",
                          "cDNA_position","CDS_position","Protein_position","Amino_acids","Codons",
                          "kg_AF","kg_AFR_AF","kg_AMR_AF","kg_EAS_AF","kg_EUR_AF","kg_SAS_AF",
                          "gnomAD_AF","gnomAD_AFR_AF","gnomAD_AMR_AF","gnomAD_ASJ_AF",
                          "gnomAD_EAS_AF","gnomAD_FIN_AF","gnomAD_NFE_AF","gnomAD_OTH_AF","gnomAD_SAS_AF","MAX_AF","MAX_AF_POPS",
                          "CADD_PHRED","SIFT_call","SIFT_score","PolyPhen_call","PolyPhen_score",
                          "CLNSIG","CLNDN","ALLELEID","CLNHGVS","CLNREVSTAT",
                          "ClinVar","LoF","FIM")

fixed.df <- fixed.df[,selected_annotations]

cat("Size of the updated object: \n")
cat("fixed.df:",dim(fixed.df),"\n")

rm(selected_annotations)

#####################################################################
#                                Save                               #
#####################################################################

cat("Save data \n\n")
save(fixed.df, gt_add.mx, dp.mx, gq.mx, ref.mx, alt.mx, 
     file=paste(base_folder,"data/S12_explore_and_filter/s04_filter_by_functional_consequence.RData",sep="/"))

#####################################################################
#                            Final section                          #
#####################################################################

cat("List of objects:\n")
ls()
cat("\n")

cat("sessionInfo:\n")
sessionInfo()
cat("\n")

cat("Memory use:\n")
gc()
cat("\n")

cat("Completed all tasks\n")
Sys.time()
cat("\n")
