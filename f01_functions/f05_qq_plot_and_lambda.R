# Make QQ-plot and estimate genomic inflation factor lambda 

# Started: Alexey Larionov, 24Jun2019  
# Last updated: Alexey Larionov, 24Jun2019  

# Stops if length p_obs != p_exp
# Stops if there are p-values <0 or >1
# Removes pairs where one of p values is NA
# Allows to remove zero p-values, or substitute them by min_p / 10
# Allows to remove or to keep p==1

qq_plot_and_lambda <- function(p_obs, p_exp, p_exp_05, p_exp_95, exclude_zeroes=F, exclude_ones=F, main="QQ plot"){

  # --- Check length of input data --- #
  
  if(length(p_obs) != length(p_exp)) stop("p_obs and p_exp should  the same length")
  if(length(p_obs) <= 10) stop("Length of p-value vectors should be > 10")
  
  # --- Stop if there are p-values < 0 --- #
  
  neg_p_obs <- p_obs < 0
  neg_p_exp <- p_exp < 0
  neg_p_num <- sum(neg_p_obs) + sum(neg_p_exp)
  
  if(neg_p_num > 0) stop(paste("Detected",neg_p_num,"negative p-values"))
  
  # --- Stop if there are p-values > 1 --- #
  
  # Check if p is significantly > 1 (more than a computing approximation error)
  large_p_obs <- p_obs > 1.00001
  large_p_exp <- p_exp > 1.00001
  large_p_num <- sum(large_p_obs) + sum(large_p_exp)
  
  if(large_p_num > 0) stop(paste("Detected",large_p_num,"p-values > 1"))
  
  # Update p if its just slightly > 1
  1 -> p_obs[p_obs>1]
  1 -> p_exp[p_exp>1]
  
  # --- Exclude the pair, if one of the pair is NA --- #
  
  missed <- is.na(p_obs) | is.na(p_exp)
  if(sum(missed)>0){
    warning(paste("Excluded",sum(missed),"pairs of p-values because of missed data"))
    p_obs <- p_obs[!missed]
    p_exp <- p_exp[!missed]
    if(length(p_obs) <= 10) stop("Number of valid pairs of p-values should be > 10")
  }
  
  # --- Deal with zero p-values --- #
  
  zero_obs <- p_obs == 0
  zero_exp <- p_exp == 0
  
  if(exclude_zeroes){
    zero_obs_or_exp <- zero_obs | zero_exp
    p_obs <- p_obs[! zero_obs_or_exp]
    p_exp <- p_exp[! zero_obs_or_exp]   
    if(length(p_obs) <= 1) stop("Number of non-emty non-zero pairs of p-values should be > 1")
    if(sum(zero_obs_or_exp) > 0) warning(paste("Detected and removed", sum(zero_obs_or_exp), "pairs of p-values with p=0"))
  }else{
    zero_substitute <- min(p_exp, p_obs) / 10
    p_obs[zero_obs] <- zero_substitute
    p_exp[zero_exp] <- zero_substitute
    if(sum(zero_obs) > 0) warning(paste("Detected ", sum(zero_obs), "zeroes in p_obs, substituted to", zero_substitute))
    if(sum(zero_exp) > 0) warning(paste("Detected ", sum(zero_exp), "zeroes in p_exp, substituted to", zero_substitute))
  }

# --- Deal with p-values = 1 --- #

  one_obs <- p_obs == 1
  one_exp <- p_exp == 1
  
  if(exclude_ones){
    one_obs_or_exp <- one_obs | one_exp
    p_obs <- p_obs[! one_obs_or_exp]
    p_exp <- p_exp[! one_obs_or_exp]   
    if(length(p_obs) <= 1) stop("Number of non-emty non-one pairs of p-values should be > 1")
    if(sum(one_obs_or_exp) > 0) warning(paste("Detected and removed", sum(one_obs_or_exp), "pairs of p-values with p=1"))
  }else{
    if(sum(one_obs) > 0) warning(paste("Detected and kept", sum(one_obs), "p=1 in p_obs"))
    if(sum(one_exp) > 0) warning(paste("Detected and kept", sum(one_exp), "p=1 in p_exp"))
  }
  
  # --- Check the number of remained pairs of p-values
  
  if(length(p_obs) <= 10) stop("Number of valid pairs of p-values should be > 10")
  
  # --- Prepare data for plot --- #
  
  neg_log_p_obs <- -log10(p_obs)
  neg_log_p_exp <- -log10(p_exp)
  neg_log_p_exp_05 <- -log10(p_exp_05)
  neg_log_p_exp_95 <- -log10(p_exp_95)
  
  #ax_lim <- c(0, max(neg_log_p_obs, neg_log_p_exp))
  
  fit <- lm(neg_log_p_obs ~ 0 + neg_log_p_exp)
  
  lambda_est <- summary(fit)$coefficients[1,"Estimate"]
  lambda_se <- summary(fit)$coefficients[1,"Std. Error"]
  
  #alt_lambda_est <- median(p_obs) / median(p_exp) ???
  
  legend_title = vector("expression", 1)
  legend_title[1] = substitute(expression(lambda == L_EST %+-% L_SE), 
                               list(L_EST = round(lambda_est,2), 
                                    L_SE = round(lambda_se,2)))[2]
  
  # --- Draw QQ-plot --- #
  
  # Make plot
  plot(neg_log_p_obs ~ neg_log_p_exp, #xlim=ax_lim, ylim=ax_lim,
       xlab="-log10(p_exp)", ylab="-log10(p_obs)", 
       col="blue",pch=20,
       main=main)
  
  # Add diagonal and fit
  abline(0,1,lty=2)
  abline(fit,col="blue")

  # Add conf intervals
  sm_05 <- smooth.spline(neg_log_p_exp_05 ~ neg_log_p_exp, spar = 2/3)
  sm_95 <- smooth.spline(neg_log_p_exp_95 ~ neg_log_p_exp, spar = 2/3)
  lines(sm_05, lty=3, lwd=0.5)
  lines(sm_95, lty=3, lwd=0.5)

  legend("topleft",
         legend=c("expected under the null","observed"),col=c("black","blue"),lty=c(2,1),
         title=legend_title,
         bty="n")
    
  # Compile and return lambda_est, lambda_se and an alternative estimate for lambda
  result <- c(lambda_est=round(lambda_est,3), lambda_se=round(lambda_se,3))
  
  return((result))
  
}