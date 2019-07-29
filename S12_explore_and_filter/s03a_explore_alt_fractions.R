# s03a_explore_alt_fractions.R
# last updated: CC23May2019  

# Run on HPC R-3.5.1 requesting 1 core and 120GB RAM for up to 2hrs

# - evaluate concordance of the sum of alt fractions with dp
# - evaluate distributions and select filtering thresholds for allelic fractions

# Suggested filters  
# - dp != alt + ref (<<1%% of genotypes)
# - for HomRefs: Alt reads fraction >5% (<<1%% of HomRefs)
# - for Hets: Alt reads fraction 25-75% (~5%% of Hets)
# - for HomAlts: Alt reads fraction <95% (~2% of HomAlts)

# For simplicity, we dont explore additional thresholds for AF counts here

# Input and output data: 73,907 variants x 12,449 cases  

#####################################################################
#                           Start section                           #
#####################################################################

# Start messages
cat("Explore fractions of Alt reads \n\n")
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
load(paste(base_folder, "data/S11_VCF_to_R/s01_dp_mx.RData", sep="/"))
load(paste(base_folder, "data/S11_VCF_to_R/s01_alt_mx.RData", sep="/"))
load(paste(base_folder, "data/S11_VCF_to_R/s01_ref_mx.RData", sep="/"))

cat("Loaded source data \n")
Sys.time()
cat("Size of the loaded objects: \n")
cat("gt_add.mx:",dim(gt_add.mx),"\n")
cat("dp.mx:",dim(dp.mx),"\n")
cat("alt.mx:",dim(alt.mx),"\n")
cat("ref.mx:",dim(ref.mx),"\n\n")

#####################################################################
#                Explore excess of dp over alt+ref                  #
#####################################################################

# Make matrix of indices for non-NA genotypes
gt_ok <- ! is.na(gt_add.mx)
num_gt_ok <- sum(gt_ok)

# Calculate alt+ref and excess of dp over alt+ref
alt_plus_ref.mx <- alt.mx + ref.mx
dp_excess_abs.mx <- dp.mx - alt_plus_ref.mx
dp_excess_rel.mx <- dp.mx / alt_plus_ref.mx

# dp is equal to alt + ref in absolute majotiry of cases
cat("Proportion of genotypes whith dp == alt+ref:", sum(dp_excess_abs.mx == 0, na.rm=T) / num_gt_ok, "\n\n")

# In a minorIty cases dp is larger than alt+ref

cat("Quartiles for absolute dp excess over alt+ref \n")
quantile(dp_excess_abs.mx, na.rm = T)
cat("\n")

cat("Quartiles for relative dp excess over alt+ref \n")
quantile(dp_excess_rel.mx, na.rm = T)
cat("\n")

cat("Made histograms for dp excess over alt+ref \n\n")

png(filename = "s03a_dp_excess_abs_hist_hr.png", width = 3000, height = 3000, res=300)
  dp_excess_abs.hist <- hist(dp_excess_abs.mx, main="Excess of dp over ref+alt (absolute)")
dev.off()

png(filename = "s03a_dp_excess_rel_hist_hr.png", width = 3000, height = 3000, res=300)
  dp_excess_rel.hist <- hist(dp_excess_rel.mx, main="Excess of dp over ref+alt (relative)")
dev.off()

# Clean-up
rm(dp.mx, dp_excess_abs.mx, dp_excess_rel.mx)

#####################################################################
#                     Alt fraction by genotype                      #
#####################################################################

# Prepare additional matrices
alt_fraction.mx <- alt.mx / alt_plus_ref.mx
hom_refs <- gt_add.mx == 0 & gt_ok
hets <- gt_add.mx == 1 & gt_ok
hom_alts <- gt_add.mx == 2 & gt_ok

cat("Prepared the following additional matrices:\n")
cat("alt_fraction.mx:",dim(alt_fraction.mx),"\n")
cat("hom_refs:",dim(hom_refs),"\n")
cat("hets:",dim(hets),"\n")
cat("hom_alts:",dim(hom_alts),"\n\n")

cat("Quartiles of alt fraction in HomRef genotypes:\n")
quantile(alt_fraction.mx[hom_refs], na.rm=T)
cat("\n")

cat("Quartiles of alt fraction in Het genotypes:\n")
quantile(alt_fraction.mx[hets], na.rm=T)
cat("\n")

cat("Quartiles of alt fraction in HomAlt genotypes:\n")
quantile(alt_fraction.mx[hom_alts], na.rm=T)
cat("\n")

cat("Made histograms for Alt fractions in different genotypes\n\n") 

png(filename = "s03a_alt_fr_HomRefs_hist_hr.png", width = 3000, height = 3000, res=300)
  alt_fr_HomRefs.hist <- hist(alt_fraction.mx[hom_refs], xlim=c(0,1))
  abline(v=0.05, lty=2, col="red")
dev.off()

png(filename = "s03a_alt_fr_Hets_hist_hr.png", width = 3000, height = 3000, res=300)
  alt_fr_Hets.hist <- hist(alt_fraction.mx[hets], xlim=c(0,1))
  abline(v=0.15, lty=2, col="red")
  abline(v=0.85, lty=2, col="red")
dev.off()

png(filename = "s03a_alt_fr_HomAlts_hist_hr.png", width = 3000, height = 3000, res=300)
  alt_fr_HomAlts.hist <- hist(alt_fraction.mx[hom_alts], xlim=c(0,1))
  abline(v=0.95, lty=2, col="red")
dev.off()

cat("Fractions of genotypes that would be removed at suggested thresholds:\n")
cat("Fraction of HomRefs with > 5% of Alt reads:", sum(alt_fraction.mx[hom_refs] > 0.05, na.rm = T) / sum(hom_refs),"\n")
cat("Fraction of Hets with < 15% or >85% of Alt reads:", sum(alt_fraction.mx[hets] < 0.15 | alt_fraction.mx[hets] > 0.85, na.rm = T) / sum(hets), "\n")
cat("Fraction of HomAlts with < 95% of Alt reads:", sum(alt_fraction.mx[hom_alts] < 0.95, na.rm = T) / sum(hom_alts),"\n\n")

#####################################################################
#                            Save plots                             #
#####################################################################

cat("Saved plots data \n\n")
save(dp_excess_abs.hist, dp_excess_rel.hist, 
     alt_fr_HomRefs.hist, alt_fr_Hets.hist, alt_fr_HomAlts.hist, 
     file="s03a_plots.RData")

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
