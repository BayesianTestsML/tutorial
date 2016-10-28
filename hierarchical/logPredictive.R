logPredictive <- function  (fittedModel,testX, rho){
  #computes the log-predictive of the fitted model on the test set (testX)
  library('mvtnorm')

    buildCovarMatrix <- function (){
      #build a matrix full of sigma^2
    covarMatrix <- matrix (data = rep (rho * currentSigma^2,  instancesEachDset * instancesEachDset), nrow=instancesEachDset)
    for (i in 1:instancesEachDset){
        covarMatrix[i,i] <- currentSigma^2 
      }
    return (covarMatrix)
  }
  
  
  dsets <- nrow(testX)
  samples <- dim(fittedModel$stanResults$delta)[1]
  #initialized as a vector of 0s, one for each dset
  postLogPredictive <- rep (0,dsets) 
  #test instances available on each dset, for which to compute the log predictive
  instancesEachDset <- ncol(testX)
  
  for (currentDset in 1:dsets){
   for (currentSample in 1:samples){
     currentMu = fittedModel$stanResults$delta[currentSample,currentDset]
     #make currentMu a vector
     currentMu <- rep(currentMu, instancesEachDset)
     currentSigma = fittedModel$stanResults$sigma[currentSample,currentDset]
     covarMatrix <- buildCovarMatrix()
     postLogPredictive[currentDset] <- postLogPredictive[currentDset] + log (dmvnorm( x=testX[currentDset,], mean=currentMu,sigma = covarMatrix) )
     if (is.infinite (dmvnorm( x=testX[currentDset,], mean=currentMu,sigma = covarMatrix))) {
       browser()
     }
   } 
  }
  postLogPredictive <- postLogPredictive / samples
  return (postLogPredictive)
  
}