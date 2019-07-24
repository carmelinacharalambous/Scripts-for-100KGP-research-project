# Function for handling errorsa and warnings  
# Motivated by demo(error.catching)  
# Started: Alexey Larionov, 09Mar2017  
# Last updated: Alexey Larionov, 06Jun2019  

# Description:
# Tries to execute an expression and returns a list with 3 elements:
# - value: if expression succeeded without error (NA if error)  
# - status: (one of the 3 words: "succeeded", "warning" or "error"  
# - message: error/warning message, if generated (NA if no error/warning)  

# Notes:
# 1) The function does NOT handle multiple errors or multiple warnings 
#    or simulteneous error(s)+warning(s)*

# 2) Sometime errors are not generated when expected,  
#    this cannot be handled with tryCatch, e.g:  
#    tryCatchAdv(1/0)  

# (*) See more about double-error, which generates error(s) + warning here:  
# http://stackoverflow.com/questions/20596902/r-avoiding-restarting-interrupted-promise-evaluation-warning  

# Tests:  
# tryCatchAdv(1/2)  
# tryCatchAdv(1/"A")  
# tryCatchAdv(chisq.test(matrix(c(1,2,3,4,5,6), nrow=2)))  
# tryCatchAdv(1/0) # succeded - the result is num Inf !  

tryCatchAdv <- function(expr)
{
  
  # Initial settings
  V <- NA
  S <- "succeeded"
  M <- NA
  
  # Warning handler
  w.handler <- function(w){
    
    # Record information about warning
    S <<- "warning"
    M <<- w
    # <<- is used for assignment outside the function scope (i.e. in the external environment)
    # http://stackoverflow.com/questions/2628621/how-do-you-use-scoping-assignment-in-r
    
    # Execute expression again, suppressing warnings
    invokeRestart("muffleWarning")
    
  }
  
  # Error handler
  e.handler <- function(e){
    
    # Record information about error
    S <<- "error"
    M <<- e
    
    # Return NA as result
    return(NA)
    
  }
  
  # Try to execute the expression with the above handlers
  V <- withCallingHandlers(tryCatch(expr, error=e.handler), warning=w.handler)
  
  # Return value
  list(value = V,
       status = S,
       message = M)
}
