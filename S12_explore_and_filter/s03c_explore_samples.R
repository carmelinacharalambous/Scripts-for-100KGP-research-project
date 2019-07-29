# s03c_explore_samples.R
# last updated: CC23May2019  

# Run on HPC R-3.5.1 requesting 1 core and 30GB RAM for up to 1hrs

# - Evaluate call rate and Ts/Tv ratios per sample  

# No need in additional samples filtering (VCF files were pre-filtering by GEL)

# Input and output data: 73,907 variants x 12,449 cases  

#####################################################################
#                           Start section                           #
#####################################################################

# Start messages
cat("Explore samples \n\n")
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
setwd(paste(base_folder,"scripts/S12_explore_and_filter/plots",sep="/"))

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

load(paste(base_folder, "data/s12_explore_and_filter/s02_fixed_df.RData", sep="/"))
load(paste(base_folder, "data/s11_VCF_to_R/s01_gt_add_mx.RData", sep="/"))

cat("Loaded source data \n")
Sys.time()
cat("Size of the loaded objects: \n")
cat("fixed.df:",dim(fixed.df),"\n")
cat("gt_add.mx:",dim(gt_add.mx),"\n\n")

#####################################################################
#                      Evaluate samples call rates                  #
#####################################################################

# Function to calculate call rate in a vector: returns proportion of non-NA values
call_rate.udf <- function(x){ sum(!is.na(x)) / length(x)}

# Calculate call rates per variant
smp_call_rates <- apply(gt_add.mx, 2, call_rate.udf)

cat("Quartiles of samples call rate: \n")
quantile(smp_call_rates)
cat("\n")

cat("Made samples call rate histogram \n\n")
png(filename = "s03c_smp_call_rates_hist_hr.png", width = 3000, height = 3000, res=300)
  smp_call_rates.hist <- hist(smp_call_rates)
dev.off()

# Clean-up
rm(call_rate.udf, smp_call_rates)

#####################################################################
#                  Calculate TsTv ratio per sample                  #
#####################################################################

# ------ Load TsTv function ------ #

source(paste(base_folder,"scripts/s12_explore_and_filter/f01_TsTv_ratio.R",sep="/"))

# ------ Calculate Ts_Tv ------ #

# Each SNP observed in the sample is considered as one transition or one transversion,  
# independently of whether this SNP is a homozygous or heterozygous.  
# Non-SNPs are ignored by the Ts/TV function (loaded in the start section).  

# Prepare empty matrix for writing results
samples_ts_tv.mx <- matrix(ncol=4,nrow=0)
colnames(samples_ts_tv.mx) <- c("case","TsTv","Ts","Tv")

# For each sample
for(sample in colnames(gt_add.mx)){
  
  # Get variants
  vars <- gt_add.mx[,sample] != 0

  # Calculate TsTv
  Ref <- fixed.df[vars,"REF"]
  Alt <- fixed.df[vars,"ALT"]
  ts_tv <- TsTv_ratio(Ref,Alt)
  
  # Write result
  samples_ts_tv.mx <- rbind(samples_ts_tv.mx,c(sample, ts_tv))

}

# Check result
cat("Calculated matrix with Ts/Tv ratios \n")
cat("Size of the matrix: \n")
dim(samples_ts_tv.mx)
cat("\n")

# Clean-up
rm(fixed.df, gt_add.mx, sample, ts_tv, Ref, Alt, TsTv_ratio, vars)

# ------ Evaluate Ts_Tv ------ #

# Convert to data frame (and to numeric where appropriate)
samples_ts_tv.df <- as.data.frame(samples_ts_tv.mx)
samples_ts_tv.df$TsTv <- as.numeric(samples_ts_tv.mx[,2])
samples_ts_tv.df$Ts <- as.numeric(samples_ts_tv.mx[,3])
samples_ts_tv.df$Tv <- as.numeric(samples_ts_tv.mx[,4])

# Summary of TsTv
cat("Quartiles of Ts/Tv ratios: \n")
quantile(samples_ts_tv.df$TsTv)
cat("\n")

# Make TsTv plots
cat("Made histogram and plot of Ts/Tv ratios \n\n")

png(filename = "s03c_smp_TsTv_hist_hr.png", width = 3000, height = 3000, res=300)
  smp_TsTv.hist <- hist(samples_ts_tv.df$TsTv)
dev.off()

png(filename = "s03c_smp_TsTv_plot_hr.png", width = 3000, height = 3000, res=300)
  plot(samples_ts_tv.df$TsTv, main="TsTv ratios per sample", xlab="samples", xaxt='n')
  abline(h=2, col="red", lty=2, lwd=3)
dev.off()

# Clean-up
rm(samples_ts_tv.df, samples_ts_tv.mx, base_folder)
  
#####################################################################
#                            Save plots                             #
#####################################################################

cat("Saved plots data \n\n")
save(smp_call_rates.hist, smp_TsTv.hist, file="s03c_plots.RData")

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
