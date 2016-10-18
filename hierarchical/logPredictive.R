logPredictive <- function  (fittedModel,testX, rho){
  #computes the log-predictive of the fitted model on the test set (testX)
  library('mvtnorm')

    buildCovarMatrix <- function (){
    covarMatrix <- matrix (nrow=instancesEachDset, ncol = instancesEachDset)
    for (i in 1:instancesEachDset){
      for (j in 1:instancesEachDset){
        covarMatrix[i,j] <- ifelse(i==j, currentSigma^2, currentSigma^2 * rho)
      }
    }
    return (covarMatrix)
  }
  
  
  dsets <- nrow(testX)
  samples <- dim(fittedModel$stanResults$delta)[1]
  postLogPredictive <- vector (length = samples)
  instancesEachDset <- ncol(testX)
  
  for (currentDset in 1:dsets){
   for (currentSample in 1:samples){
     currentMu = fittedModel$stanResults$delta[currentSample][currentDset]
     #make currentMu a vector
     currentMu <- rep(currentMu, instancesEachDset)
     currentSigma = fittedModel$stanResults$sigma[currentSample][currentDset]
     covarMatrix <- buildCovarMatrix()
     postLogPredictive[currentSample] <- log (dmvnorm( x=testX[currentDset,], mean=currentMu,sigma = currentSigma) )
   } 
  }
  return (sum(postLogPredictive))
  
}