# Function for calculatintg Fisher exact p-values and the corresponding  
# p-values expected under the null (obtained by permutation of case/control lables) 

# Started: Alexey Larionov, 22Jun2019
# Last updated: Alexey Larionov, 22Jun2019

# genotypes.mx is a matrix of additive genotypes (NA,0,1,2) 
# with genes in rows, samples in columns

# is_case is vector of case/control lables (0,1)

# n_perm is the number of required permutations

# Takes ~15 min to permute 1000 times a matrix of 2,000 variants x 500 samples

# Consider additional checks with source data and intolerance to warnings ...

fisher_variant_qq_perm <- function(genotypes.mx, cases, n_perm=1000, show_progress=F){
  
  # To show progress (not used by default)
  if(show_progress){library(svMisc)}
  
  # -------------------------------------------------------------- #
  # Function to get sorted vector of p-values for a given additive #
  #   genotypes matrix and a given vector of case/control lables   #
  # -------------------------------------------------------------- #
  
  get_or_and_p_values <- function(genotypes, is_case){
    
    # Convert cases to logical
    is_case <- as.logical(is_case)
    
    # Make an empty vectors for odd-ratios and p-values for each variant
    p_values <- rep(NA, nrow(genotypes))
    or_values <- rep(NA, nrow(genotypes))
    
    # For each variant (row) in the matrix
    for(var in 1:nrow(genotypes)){
      
      # Get vector of genotypes for the variant
      gtv <- genotypes[var,]
      
      # Count Alt and Ref alleles in cases
      cases_alt <- sum(gtv[is_case],na.rm=T)
      cases_an <- 2*sum(!is.na(gtv[is_case]))
      cases_ref <- cases_an - cases_alt
      
      # Count Alt and Ref alleles in controls
      controls_alt <- sum(gtv[!is_case],na.rm=T)
      controls_an <- 2*sum(!is.na(gtv[!is_case]))
      controls_ref <- controls_an - controls_alt
      
      # Calculate and add Fisher p-value to the vector of p-values
      ft <- fisher.test(matrix(c(cases_alt,cases_ref,controls_alt,controls_ref),nrow=2))
      
      p_values[var] <- ft$p.value
      or_values[var] <- ft$estimate
    }
    
    # Add names, sort p-values, sync or-values
    names(p_values) <- rownames(genotypes)
    names(or_values) <- rownames(genotypes)
    p_values <- sort(p_values)
    or_values <- or_values[names(p_values)]
    
    return(list(p_values=p_values,or_values=or_values))
  }
  
  # -------------------------------------------------------------- #
  #                      Sorted observed p-values                  #
  # -------------------------------------------------------------- #
  
  observed <- get_or_and_p_values(genotypes.mx, cases)
  
  # -------------------------------------------------------------- #
  #  Sorted p-values expected under the null (by permuting lables) #
  # -------------------------------------------------------------- #
  
  # Make empty matrix to accumulate p-values during permutations
  p.mx <- matrix(NA, ncol=nrow(genotypes.mx), nrow=0)
  
  # Accumulate sorted permuted p-values for the required number of permutations
  for(p in 1:n_perm){
    
    # Permute lables using sample() function, 
    # calculate the vector of Fisher P-values using local get_p_values() function
    # add the vector to permutations matrix
    permuted <- get_or_and_p_values(genotypes.mx,sample(cases))
    p.mx <- rbind(p.mx, permuted$p_values)
    
    # Show progress (requires library svMisc, not used by default)
    if(show_progress){progress(value = p, max.value = n_perm)}
  }
  
  # Calculate expected p-values as mean (and 95CI) of the permuted p-values
  p_exp <- apply(p.mx,2,mean,na.rm=T)
  p_exp_95pp <- apply(p.mx,2,quantile,probs=0.95,na.rm=T)
  p_exp_05pp <- apply(p.mx,2,quantile,probs=0.05,na.rm=T)
  names(p_exp) <- NULL
  names(p_exp_95pp) <- NULL
  names(p_exp_05pp) <- NULL
  
  # Compile "expected" as a list with the central expectation of p-values and 95CI
  expected <- list(p_values=p_exp,p_exp_05pp=p_exp_05pp,p_exp_95pp=p_exp_95pp)
  
  # Return observed (or and p) and expected (p)
  return(list(observed=observed, expected=expected))
  
}
