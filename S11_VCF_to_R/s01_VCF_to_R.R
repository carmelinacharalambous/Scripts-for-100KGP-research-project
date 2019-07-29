# s01_VCF_to_R.R
# last updated: CC20May2019  

# Run on HPC R-3.5.1 requesting 1 core and 60GB RAM for up to 24hrs
# Assumes that multi-allelic variants are split into separate lines !  

# Expected source data:
# Multi-sample VCF with split multi-allelic variants annotated with VEP
# containing a unique variant ID in INFO field "VarID"  

# Read vcf to R:  
# - Read VCF header to meta.df  
# - Read standard VCF fields to fixed.df  
# - Parse VEP  
# - Read genotypes to numeric (1/0) and character (A/T) matrices  
# - Read gq, dp, ad  
# - Recode genotypes to additive format  
# - Split matrix with allelic depth  
# - Update rownames to VarID  
# - Save as set of RData objects: one for each matrix/dataframe  
# - Save as a full RData object, which contains all the objects  

#####################################################################
#                           Start section                           #
#####################################################################

# Start message
cat("Read vcf to R \n")
Sys.time()
version$version.string
gc()
cat("\n")

# Clean-up
rm(list=ls())
graphics.off()

# Lib folders
additional_personal_lib_folder <- "/re_gecip/inherited_cancer_predisposition/ccharalambous/Final_execution/data/k02_inuvika_r_libs/R_3.5.1"
additional_system_lib_folder <- "/public_data_resources/Rpackages/3.5.1"

.libPaths(c(additional_personal_lib_folder, 
            additional_system_lib_folder, 
            .libPaths()))

cat("Attached lib folders: \n")
.libPaths()
cat("\n")

# Load ibraries
library(vcfR)
library(stringr) # for word (for VEP parsing)
library(tidyr) # for separate (for VEP parsing)
library(dplyr) # for renaming AC, AN and AF

# Working folder
base_folder <- "/re_gecip/inherited_cancer_predisposition/ccharalambous/Final_execution"
setwd(base_folder)

# Additional options
options(stringsAsFactors = F,
        warnPartialMatchArgs = T, 
        warnPartialMatchAttr = T, 
        warnPartialMatchDollar = T)

# Progress report
cat("Completed start section \n")
Sys.time()
cat("\n")

# Clean-up
rm(additional_personal_lib_folder, additional_system_lib_folder)

#####################################################################
#                      Read VCF to vcfR object                      #
#####################################################################

# Source VCF
vcf_file=paste(base_folder, "data/S10_fill_tags/s01_VEP_clinvar_VarID_tags.vcf.gz", sep="/")

# Read vcf to vcfR object
vcfr <- read.vcfR(vcf_file)

# Progress report
cat("Completed reading vcfR \n")
Sys.time()
cat("Size of vcfR object: \n")
dim(vcfr)
cat("\n")

# Clean-up
rm(vcf_file, source_data_folder)

#####################################################################
#    extract meta (vcf header) and fixed (8 standard VCF columns)   #
#####################################################################

# Get data from VCF header and fixed/standard VCF columns
# Warning messages 
# In lapply(ret[names(ns)], as.numeric) : NAs introduced by coercion
# may be generated and ignored at this step
meta_fix <- vcfR2tidy(vcfr, info_only=T)
    
# Get data frame with meta-information from vcf header
meta.df <- meta_fix$meta

# Get data frame with fixed columns (including parsed INFO, convert tibble to data-frame)
fixed.df <- as.data.frame(meta_fix$fix)

# Progress report
cat("Completed extracting meta.df and fixed.df \n")
Sys.time()
cat("Sizes of new dataframes: \n")
cat("meta.df:",dim(meta.df),"\n")
cat("fixed.df:",dim(fixed.df),"\n\n")

# Clean-up
rm(meta_fix)

#####################################################################
#           Rename initial AC and AN, calculate initial AF          #
#####################################################################

# The renaimimg is needed to avoid confusion between different ACs, ANs and AFs later
# for instance, VEP may annotate 1000genomes AF as just "AF" etc

fixed.df <- fixed.df %>% 
  rename(Initial_AC = AC,
         Initial_AN = AN,
         Initial_AF = AF)

# Progress report
cat("Renamed initial AC and AN, calculated initial AF \n")
Sys.time()
cat("\n")

#####################################################################
#                        Parse VEP (CSQ) column                     #
#####################################################################

# Make sure that the VEP annotation field name is CSQ (it may be ANN or something else in some cases)
# Note fixed=T in strsplit: otherwise it would interpreted split as regex

vep_fields <- as.character(meta.df[meta.df$ID=="CSQ","Description"])
vep_fields <- word(vep_fields,-1) # requires stringr
vep_fields <- strsplit(vep_fields, "|", fixed=T)
vep_fields <- unlist(vep_fields)

# Split VEP CSQ column
# Note \\ in sep: this is because the separator is interpreted as regex
fixed.df <- separate(fixed.df, "CSQ", vep_fields, sep="\\|") 

# Progress report
cat("Completed parsing VEP annotation \n")
Sys.time()
cat("Detected", length(vep_fields), "VEP fields \n\n")

# Clean-up
rm(vep_fields)

#####################################################################
#                         Rename VEP AF fileds                      #
#####################################################################

# This renaiming is to avoid confusion between different "AF"-s 

# --- AFs from 1000-genomes project --- #

# It is especially important to rename "AF" field from 1000-genomes !
# This field is emitted by --af VEP option (included in --everything). 
# Unfortunately, VEP decided that AF field should contain the global allele frequency (AF) from 1000 Genomes Phase 3.  
# However, AF may also be used in many other contexts, so it need to be specifyed unambiguously

fixed.df <- fixed.df %>% 
  rename(kg_AF = AF,
         kg_AFR_AF = AFR_AF, 
         kg_AMR_AF = AMR_AF,
         kg_EAS_AF = EAS_AF,
         kg_EUR_AF = EUR_AF,
         kg_SAS_AF = SAS_AF)

# --- AFs from ESP project --- #

fixed.df <- fixed.df %>% 
  rename(esp_AA_AF = AA_AF,
         esp_EA_AF = EA_AF)

# Progress report
cat("Renamed VEP AF fields\n")
Sys.time()
cat("\n")

#####################################################################
#            Recode missed values in VEP annotations to NAs         #
#####################################################################

# Recode blanks as NAs (dot "." might also be checked)
NA -> fixed.df[ fixed.df=="" ] 

# Progress report
cat("Recoded missed VEP values to NA \n")
Sys.time()
cat("\n")

#####################################################################
#                 Extract matrices: gt, gq, dp and ad               #
#####################################################################

# Assuming that multi-allelic variants were split into separate lines !  

gt_num.mx <- extract.gt(vcfr) # numeric codes: 0/1, 1/1 etc
gt_chr.mx <- extract.gt(vcfr, return.alleles = TRUE) # character codes: A/A, T/G etc
NA -> gt_chr.mx [ gt_chr.mx=="." ]

dp.mx <- extract.gt(vcfr, element = "DP", as.numeric = TRUE)
gq.mx <- extract.gt(vcfr, element = "GQ", as.numeric = TRUE)
ad.mx <- extract.gt(vcfr, element = "AD")

# Progress report
cat("Completed extracting matrices: gt_num, gt_chr, gq, dp and ad\n")
Sys.time()
cat("Sizes of the extracted matrices:\n")
cat("gt_num.mx:",dim(gt_num.mx),"\n")
cat("gt_chr.mx:",dim(gt_chr.mx),"\n")
cat("dp.mx:",dim(dp.mx),"\n")
cat("gq.mx:",dim(gq.mx),"\n")
cat("ad.mx:",dim(ad.mx),"\n\n")

# Clean-up
rm(vcfr)

#####################################################################
#   Make matrices with additive, dominant and recessive genotypes   #
#####################################################################

#Assuming that multi-allelic variants are split into separate lines !  

# Recode to additive
gt_add.mx <- gt_num.mx
0 -> gt_add.mx[ gt_num.mx == "0/0" ]
1 -> gt_add.mx[ gt_num.mx == "1/0" ]
1 -> gt_add.mx[ gt_num.mx == "0/1" ]
2 -> gt_add.mx[ gt_num.mx == "1/1" ]

# Convert to numeric
gt_add.mx <- matrix(as.numeric(gt_add.mx),nrow=nrow(gt_add.mx))
colnames(gt_num.mx) -> colnames(gt_add.mx)
rownames(gt_num.mx) -> rownames(gt_add.mx)

# Recode to dominant
gt_dom.mx <- gt_num.mx
0 -> gt_dom.mx[ gt_num.mx == "0/0" ]
1 -> gt_dom.mx[ gt_num.mx == "1/0" ]
1 -> gt_dom.mx[ gt_num.mx == "0/1" ]
1 -> gt_dom.mx[ gt_num.mx == "1/1" ]

# Convert to numeric
gt_dom.mx <- matrix(as.numeric(gt_dom.mx),nrow=nrow(gt_dom.mx))
colnames(gt_num.mx) -> colnames(gt_dom.mx)
rownames(gt_num.mx) -> rownames(gt_dom.mx)

# Recode to recessive 
gt_rec.mx <- gt_num.mx
0 -> gt_rec.mx[ gt_num.mx == "0/0" ]
0 -> gt_rec.mx[ gt_num.mx == "1/0" ]
0 -> gt_rec.mx[ gt_num.mx == "0/1" ]
1 -> gt_rec.mx[ gt_num.mx == "1/1" ]

# Convert to numeric
gt_rec.mx <- matrix(as.numeric(gt_rec.mx),nrow=nrow(gt_rec.mx))
colnames(gt_num.mx) -> colnames(gt_rec.mx)
rownames(gt_num.mx) -> rownames(gt_rec.mx)

# Progress report
cat("Calculated additive, dominant and recessive genotype matrices \n")
Sys.time()
cat("Sizes of the new matrices:\n")
cat("gt_add.mx:",dim(gt_add.mx),"\n")
cat("gt_dom.mx:",dim(gt_dom.mx),"\n")
cat("gt_rec.mx:",dim(gt_rec.mx),"\n\n")

#####################################################################
#                   Split ad (allelic depth) matrix                 #
#####################################################################

# Assuming that multi-allelic variants are split into separate lines !  

# Convert NA (if any) to ".,."
".,." -> ad.mx[is.na(ad.mx)]

# Convert ad matrix to vector
ad <- unlist(strsplit(ad.mx,",")) # split by comma

# Extract vectors for ref and alt elements
ref <- ad[seq(1,length(ad),2)] # uneven elements
alt <- ad[seq(2,length(ad),2)] # even elements

# Convert ref and alt vectors to matrices
ref.mx <- matrix(as.integer(ref), nrow=nrow(ad.mx))
alt.mx <- matrix(as.integer(alt), nrow=nrow(ad.mx))

# Preserve row- and col- names
rownames(ref.mx) <- rownames(ad.mx)
colnames(ref.mx) <- colnames(ad.mx)

rownames(alt.mx) <- rownames(ad.mx)
colnames(alt.mx) <- colnames(ad.mx)

# Progress report
cat("Completed splitting ad matrix \n")
Sys.time()
cat("Sizes of the new matrices:\n")
cat("ref.mx:",dim(ref.mx),"\n")
cat("alt.mx:",dim(alt.mx),"\n\n")

# Clean-up
rm(ad, ref, alt)

#####################################################################
#                           Update rownames                         #
#####################################################################

# Update this section, if there is no VarID field in INFO !

# Generate row names
row_names <- fixed.df$VarID

# Assign the rownames
row_names -> rownames(fixed.df)
row_names -> rownames(gt_num.mx)
row_names -> rownames(gt_chr.mx)
row_names -> rownames(gt_add.mx)
row_names -> rownames(gt_dom.mx)
row_names -> rownames(gt_rec.mx)
row_names -> rownames(dp.mx)
row_names -> rownames(gq.mx)
row_names -> rownames(ad.mx)
row_names -> rownames(ref.mx)
row_names -> rownames(alt.mx)

# Progress report
cat("Updated row names \n")
Sys.time()
cat("\n")

# Clean-up
rm(row_names)

#####################################################################
#                             Save results                          #
#####################################################################

# Save separate files
target_data_folder <- paste(base_folder,"data/S11_VCF_to_R",sep="/")

save(meta.df,file=paste(target_data_folder,"s01_meta_df.RData",sep="/"))
save(fixed.df,file=paste(target_data_folder,"s01_fixed_df.RData",sep="/"))
save(gt_num.mx,file=paste(target_data_folder,"s01_gt_num_mx.RData",sep="/"))
save(gt_chr.mx,file=paste(target_data_folder,"s01_gt_chr_mx.RData",sep="/"))
save(gt_add.mx,file=paste(target_data_folder,"s01_gt_add_mx.RData",sep="/"))
save(gt_dom.mx,file=paste(target_data_folder,"s01_gt_dom_mx.RData",sep="/"))
save(gt_rec.mx,file=paste(target_data_folder,"s01_gt_rec_mx.RData",sep="/"))
save(dp.mx,file=paste(target_data_folder,"s01_dp_mx.RData",sep="/"))
save(gq.mx,file=paste(target_data_folder,"s01_gq_mx.RData",sep="/"))
save(ad.mx,file=paste(target_data_folder,"s01_ad_mx.RData",sep="/"))
save(ref.mx,file=paste(target_data_folder,"s01_ref_mx.RData",sep="/"))
save(alt.mx,file=paste(target_data_folder,"s01_alt_mx.RData",sep="/"))

# Save as one image containing all objects
setwd(target_data_folder)
rm(target_data_folder,base_folder)
save.image("s01_VCF_to_R.RData")

# Progress report
cat("Saved results\n")
Sys.time()
cat("\n")

#####################################################################
#                            Final section                          #
#####################################################################

# Log information about environment

cat("List of",length(ls()),"saved objects:\n")
ls()
cat("\n")

cat("sessionInfo:\n")
sessionInfo()
cat("\n")

cat("Memory info:\n")
gc()
cat("\n")

# Completion message

cat("Completed all tasks\n")
Sys.time()
cat("\n")
