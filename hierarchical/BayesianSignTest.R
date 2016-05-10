BayesianSignTest <- function(diffVector,rope_min,rope_max) {

  library(MCMCpack)
  
  samples <- 3000

  #build the vector 0.5 1 1 ....... 1 
  weights <- c(0.5,rep(1,length(diffVector)))
  
  #add the fake first observation in 0
  diffVector <- c (0, diffVector)  
  
  
  #for the moment we implement the sign test. Signedrank will follows
  probLeft <- mean (diffVector < rope_min)
  probRope <- mean (diffVector > rope_min & diffVector < rope_max)
  probRight <- mean (diffVector > rope_max)
  
  
  
  results = list ("probLeft"=probLeft, "probRope"=probRope,
                  "probRight"=probRight)
  
  return (results)
  
}

