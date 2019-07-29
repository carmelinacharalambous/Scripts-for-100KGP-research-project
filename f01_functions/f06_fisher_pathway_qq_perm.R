# Function for calculatintg Fisher exact p-values and the corresponding  
# p-values expected under the null (obtained by permutation of case/control lables) 
# for variants aggregated per pathways 

# THE DESCRIPTION YET TO BE UPDATED

# Started: Alexey Larionov, 24Jun2019
# Last updated: Alexey Larionov, 24Jun2019

# genotypes.mx is a matrix of additive genotypes (NA,0,1,2) 
# with genes in rows, samples in columns

# is_case is vector of case/control lables (0,1)

# n_perm is the number of required permutations

# Takes ~15 min to permute 1000 times a matrix of 2,000 variants x 500 samples

# Consider additional checks with source data and intolerance to warnings ...

# test data:
#load("/home/alarionov/re_gecip/inherited_cancer_predisposition/tischkowitz/users/alexey/fa_analysis/data/s19_association_tests/s02_variant_tests.RData")

fisher_pathway_qq_perm <- function(genotypes.mx, variants.df, cases, n_perm=1000, show_progress=F){
  
  # To show progress (not used by default)
  if(show_progress){library(svMisc)}
  
  # -------------------------------------------------------------- #
  #  Function to get sorted vector of p-values for a given set of  #
  #            genotypes, variants and case/control lables         #
  # -------------------------------------------------------------- #
  
  get_or_and_p_values <- function(genotypes, variants, is_case){
    
    # Convert cases to logical
    is_case <- as.logical(is_case)
    
    # Get list of pathways
    pathways <- unique(variants$pathways)
    
    # Make an empty vectors for odd-ratios and p-values for each pathway
    p_values <- rep(NA,  length(pathways))
    pathways -> names(p_values) 
    or_values <- rep(NA, length(pathways))
    pathways -> names(or_values)
    
    # For each  pathway
    for(pathway in pathways){
      
      # Get variansts for the  pathway
      pathway_vars <- variants$pathways == pathway
      
      # Get matrix of genotypes for the  pathway
      pathway_data <- genotypes[pathway_vars,,drop=F]
      
      # Count Alt and Ref alleles in cases
      cases_alt <- sum( pathway_data[,is_case],na.rm=T)
      cases_an <- 2*sum(!is.na( pathway_data[,is_case]))
      cases_ref <- cases_an - cases_alt
      
      # Count Alt and Ref alleles in controls
      controls_alt <- sum( pathway_data[,!is_case],na.rm=T)
      controls_an <- 2*sum(!is.na( pathway_data[,!is_case]))
      controls_ref <- controls_an - controls_alt
      
      # Calculate and add Fisher p-value to the vector of p-values
      ft <- fisher.test(matrix(c(cases_alt,cases_ref,controls_alt,controls_ref),nrow=2))
      
      p_values[ pathway] <- ft$p.value
      or_values[ pathway] <- ft$estimate
    }
    
    # Sort p-values, sync or-values
    p_values <- sort(p_values)
    or_values <- or_values[names(p_values)]
    
    return(list(p_values=p_values,or_values=or_values))
  }
  
  # -------------------------------------------------------------- #
  #                      Sorted observed p-values                  #
  # -------------------------------------------------------------- #
  
  observed <- get_or_and_p_values(genotypes.mx, variants.df, cases)
  
  # -------------------------------------------------------------- #
  #  Sorted p-values expected under the null (by permuting lables) #
  # -------------------------------------------------------------- #
  
  # Make empty matrix to accumulate p-values during permutations
  pathways <- unique(variants.df$pathways)
  p.mx <- matrix(NA, ncol=length(pathways), nrow=0)
  
  # Accumulate sorted permuted p-values for the required number of permutations
  for(p in 1:n_perm){
    
    # Permute lables using sample() function, 
    # calculate the vector of Fisher P-values using local get_p_values() function
    # add the vector to permutations matrix
    permuted <- get_or_and_p_values(genotypes.mx,variants.df,sample(cases))
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
