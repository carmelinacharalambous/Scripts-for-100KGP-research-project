# Function to calculate Transitions to Transversions ratio
# Started: Alexey Larionov 07Sep2018
# Last updated: Carmelina Charalambous 5May2019

# Takes vectors of Ref and Alt alleles coded as A,T,C,G  

TsTv_ratio <- function(Ref,Alt){
  
  # Transitions
  Ts_AG <- sum(Ref=="A" & Alt=="G", na.rm=T)
  Ts_GA <- sum(Ref=="G" & Alt=="A", na.rm=T)
  Ts_CT <- sum(Ref=="C" & Alt=="T", na.rm=T)
  Ts_TC <- sum(Ref=="T" & Alt=="C", na.rm=T)
  
  Ts <- sum(Ts_AG, Ts_GA, Ts_CT, Ts_TC, na.rm=T)
  
  # Transversions
  Tv_AC <- sum(Ref=="A" & Alt=="C", na.rm=T)
  Tv_AT <- sum(Ref=="A" & Alt=="T", na.rm=T)
  Tv_GC <- sum(Ref=="G" & Alt=="C", na.rm=T)
  Tv_GT <- sum(Ref=="G" & Alt=="T", na.rm=T)
  Tv_CA <- sum(Ref=="C" & Alt=="A", na.rm=T)
  Tv_CG <- sum(Ref=="C" & Alt=="G", na.rm=T)
  Tv_TA <- sum(Ref=="T" & Alt=="A", na.rm=T)
  Tv_TG <- sum(Ref=="T" & Alt=="G", na.rm=T)
  
  Tv <- sum(Tv_AC, Tv_AT, Tv_GC, Tv_GT, Tv_CA, Tv_CG, Tv_TA, Tv_TG, na.rm=T)
  
  # Ratio
  if(Tv == 0){
    TsTv <- NA
  }else{
    TsTv <- Ts / Tv
  }
  
  # Compile result
  result <- c(TsTv, as.integer(Ts), as.integer(Tv))
  names(result) <- c("TsTv", "Ts", "Tv")
  
  # Return result
  return(result)
  
}
