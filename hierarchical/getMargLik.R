getMargLik <- function(hierarchicalResults,xTest,rho) {
  library(mvtnorm)
  #folds to be predicted
  foldsPred<-ncol(xTest)
  dsets<-nrow(xTest)
  samples <- length(hierarchicalResults$stanResults$delta0)
  stanResults<-hierarchicalResults$stanResults
  logLik<-matrix(ncol = dsets, nrow =   samples)
  
  
  for (currentDset in 1:dsets){
    for (currentSample in 1:samples){
      currentMean<-stanResults$delta[currentSample,currentDset]
      currentStd<-stanResults$sigma[currentSample,currentDset]
      
      currentCovar <- matrix (currentStd^2*rho,ncol=foldsPred,nrow = foldsPred)
      #build covariance matrix
      for (i in 1:foldsPred){
            currentCovar[i,i] <- currentStd^2
      }
      logLik[currentSample,currentDset] <- dmvnorm (xTest[currentDset,], rep(currentMean, length = foldsPred), currentCovar, log=TRUE)
    }
  }
  margLik <- colSums(logLik)
  results = list ("margLik"=margLik)
  return (margLik)
}