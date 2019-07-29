# s03b_explore_genotypes.R
# last updated: CC23May2019  

# Run on HPC R-3.5.1 requesting 1 core and 60GB RAM for up to 2hrs

# - evaluate missingness of genotypes  
# - evaluate dp and gq distributions and select filtering thresholds  

# Suggested filters:  
#
# - min gq 20 (removes ~3% of genotypes)  
# - min dp 10 (removes ~2% of genotypes)  

# Input and output data: 73,907 variants x 12,449 cases  

#####################################################################
#                           Start section                           #
#####################################################################

# Start messages
cat("Explore genotypes \n\n")
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

load(paste(base_folder, "data/S11_VCF_to_R/s01_gt_add_mx.RData", sep="/"))
load(paste(base_folder, "data/S11_VCF_to_R/s01_gq_mx.RData", sep="/"))
load(paste(base_folder, "data/S11_VCF_to_R/s01_dp_mx.RData", sep="/"))

cat("Loaded source data \n")
Sys.time()
cat("Size of the loaded objects: \n")
cat("gt_add.mx:",dim(gt_add.mx),"\n")
cat("gq.mx:",dim(gq.mx),"\n")
cat("dp.mx:",dim(dp.mx),"\n\n")

#####################################################################
#                   Explore variants' QC metrics                    #
#####################################################################

# ------ Genotypes NA rate ------ #

cat("Genotypes NA rate:", sum(is.na(gt_add.mx)) / (nrow(gt_add.mx) * ncol(gt_add.mx)),"\n\n")

# ------ Make matrix of indices for non-NA genotypes ------ #

gt_ok <- ! is.na(gt_add.mx)
num_gt_ok <- sum(gt_ok)

# ------ Proportions of HomRef, Het and HomAlt genotypes ------ #

cat("Proportion of HomRef genotypes:", sum(gt_add.mx==0, na.rm=T)/num_gt_ok,"\n")
cat("Proportion of Het genotypes:", sum(gt_add.mx==1, na.rm=T)/num_gt_ok,"\n")
cat("Proportion of HomAlt genotypes:", sum(gt_add.mx==2, na.rm=T)/num_gt_ok,"\n\n")

# ------ gq ------ #

cat("Quartiles of gq in non-NA genotypes: \n")
quantile(gq.mx[gt_ok], na.rm = T)
cat("\n")

# Low quality genotypes
cat("Proportion of genotypes with gq < 20:",sum(gq.mx[gt_ok] < 20, na.rm=T) / num_gt_ok,"\n\n")

cat("Made histograms of gq in non-NA genotypes \n\n")

png(filename = "s03b_gq_hist_hr.png", width = 3000, height = 3000, res=300)
  gq.hist <- hist(gq.mx[gt_ok])
dev.off()

png(filename = "s03b_gq_500_hist_hr.png", width = 3000, height = 3000, res=300)
  gq_500.hist <- hist(gq.mx[gt_ok & gq.mx < 500])
  abline(v=20, lty=2, col="red")
dev.off()

# ------ dp ------ #

cat("Quartiles of dp in non-NA genotypes: \n")
quantile(dp.mx, na.rm=T)
cat("\n")

# Low covered genotypes
cat("Proportion of genotypes with dp < 10:",sum(dp.mx[gt_ok] < 10, na.rm=T) / num_gt_ok,"\n\n")

cat("Made histogram of dp in non-NA genotypes \n\n")
png(filename = "s03b_gt_dp_hist_hr.png", width = 3000, height = 3000, res=300)
  gt_dp.hist <- hist(dp.mx[gt_ok])
  abline(v=10, lty=2, col="red")
dev.off()

# ------ Occasional Inconsistencies between gt and gq/dp ------ #

cat("Some additional checks: \n\n")
cat("Total number of genotypes:", nrow(gt_add.mx) * ncol(gt_add.mx), "\n\n")

cat("Num of genotypes missed despite the good quality (gq > 20):",sum(!gt_ok & gq.mx > 20),"\n")
cat("Num of genotypes missed despite the good coverage (dp > 10):",sum(!gt_ok & dp.mx > 10),"\n\n")

cat("Num of genotypes assigned despite gq==0:",sum(gq.mx[gt_ok] == 0, na.rm=T),"\n")
cat("Num of genotypes assigned despite dp==0:",sum(dp.mx[gt_ok] == 0, na.rm=T),"\n\n")

cat("Num of genotypes assigned despite the missed gq:",sum(is.na(gq.mx[gt_ok])),"\n")
cat("Num of genotypes assigned despite the missed dp:",sum(is.na(dp.mx[gt_ok])),"\n\n")

#####################################################################
#                            Save plots                             #
#####################################################################

cat("Saved plots data \n\n")
save(gq.hist, gq_500.hist, gt_dp.hist, file="s03b_plots.RData")

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
