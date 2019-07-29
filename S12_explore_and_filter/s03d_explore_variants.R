# s03d_explore_variants.R
# last updated: CC23May2019  

# Run on HPC R-3.5.1 requesting 1 core and 30GB RAM for up to 1hrs

# - Evaluate FILTER, QUAL, DP and Variants-Call-Rate  

# Suggested additional filter(s):  
# Remove 875 variants with mean DP per sample < 10x  

# No additional filtering by FILTER column is needed:  
# only PASS variants had been retained in the file.

# No additional filtering by QUAL column is needed:  
# filtering by QUAL had already been applied previously (for Ts/Tv >2).

# No need in additional filtering by variant call rates:  
# the call rates per variant are already high (because of the VCF file pre-filtering by GEL)  

# Input and output data: 73,907 variants x 12,449 cases

#####################################################################
#                           Start section                           #
#####################################################################

# Start messages
cat("Explore variants QC metrics \n\n")
Sys.time()
version$version.string
gc()
cat("\n")

# Clean-up
rm(list=ls())
graphics.off()

# Base folder
base_folder <- "/re_gecip/inherited_cancer_predisposition/ccharalambous/Final_execution"

# The plots will be saved to the working folder (the folder should exist):
setwd(paste(base_folder,"scripts/s12_explore_and_filter/plots",sep="/"))

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

cat("Loaded source data \n")
Sys.time()
cat("Size of the loaded objects: \n")
cat("fixed.df:",dim(fixed.df),"\n")
cat("gt_add.mx:",dim(gt_add.mx),"\n\n")

#####################################################################
#                   Explore variants' QC metrics                    #
#####################################################################

# ------ FILTER ------ #

cat("Content of FILTER field: \n")
table(fixed.df$FILTER)
cat("\n")

cat("NAs count in FILTER field: \n")
sum(is.na(fixed.df$FILTER))
cat("\n")

# ------ QUAL ------ #

cat("Quartiles of QUAL: \n")
quantile(fixed.df$QUAL)
cat("\n")

cat("Made QUAL histogram \n\n")
png(filename = "s03d_qual_hist_hr.png", width = 3000, height = 3000, res=300)
  qual.hist <- hist(fixed.df$QUAL)
dev.off()

png(filename = "s03d_qual_hist_5000_hr.png", width = 3000, height = 3000, res=300)
  qual_5000.hist <- hist(fixed.df$QUAL[fixed.df$QUAL<5000])
dev.off()

# ------ DP ------ #

cat("Quartiles of DP: \n")
quantile(fixed.df$DP)
cat("\n")

cat("Number of variants with DP<10: \n")
sum(fixed.df$DP < nrow(fixed.df)*10)
cat("\n")

cat("Made DP histogram \n\n")
png(filename = "s03d_qual_hist_5000_hr.png", width = 3000, height = 3000, res=300)
  dp.hist <- hist(fixed.df$DP)
  abline(v=12449*10, lty=2, col="red")
dev.off()

# Clean-up
rm(fixed.df)

# ------ Variants call rate ------ #

# Function to calculate call rate in a vector: returns proportion of non-NA values
call_rate.udf <- function(x){ sum(!is.na(x)) / length(x)}

# Calculate call rates per variant
var_call_rates <- apply(gt_add.mx, 1, call_rate.udf)

cat("Quartiles of variants call rate: \n")
quantile(var_call_rates)
cat("\n")

cat("Made variants call rate histogram \n\n")
png(filename = "s03d_var_call_rates_hist_hr.png", width = 3000, height = 3000, res=300)
  var_call_rates.hist <- hist(var_call_rates)
dev.off()

# Clean-up
rm(var_call_rates, gt_add.mx, call_rate.udf)

#####################################################################
#                            Save plots                             #
#####################################################################

cat("Saved plots data \n\n")
save(qual.hist, qual_5000.hist, dp.hist, var_call_rates.hist, file="s03d_plots.RData")

#####################################################################
#                            Final section                          #
#####################################################################

cat("sessionInfo:\n")
sessionInfo()
cat("\n")

cat("Memory use:\n")
gc()
cat("\n")

cat("Completed all tasks\n")
Sys.time()
cat("\n")
