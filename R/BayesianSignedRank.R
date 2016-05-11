BayesianSignedRank <- function(diffVector,rope_min,rope_max) {
  
  library(MCMCpack)
  
  samples <- 30000
  
  #build the vector 0.5 1 1 ....... 1
  weights <- c(0.5,rep(1,length(diffVector)))
  
  #add the fake first observation in 0
  diffVector <- c (0, diffVector)
  
  sampledWeights <- rdirichlet(samples,weights)
  
  winLeft <- vector(length = samples)
  winRope <- vector(length = samples)
  winRight <- vector(length = samples)
  
  for (rep in 1:samples){
    currentWeights <- sampledWeights[rep,]
    for (i in 1:length(currentWeights)){
      for (j in 1:length(currentWeights)){
        product= currentWeights[i] * currentWeights[j]
        if (diffVector[i]+diffVector[j] > (2*rope_max) ) {
          winRight[rep] <- winRight[rep] + product
        }
        else if (diffVector[i]+diffVector[j] > (2*rope_min) ) {
          winRope[rep] <- winRope[rep] + product
        }
        else {
          winLeft[rep] <- winLeft[rep] + product
        }
      }
    }
    maxWins=max(winRight[rep],winRope[rep],winLeft[rep])
    winners = (winRight[rep]==maxWins)*1 + (winRope[rep]==maxWins)*1 + (winLeft[rep]==maxWins)*1
    winRight[rep] <- (winRight[rep]==maxWins)*1/winners
    winRope[rep] <- (winRope[rep]==maxWins)*1/winners
    winLeft[rep] <- (winLeft[rep]==maxWins)*1/winners
  }
  
  
  results = list ("winLeft"=mean(winLeft), "winRope"=mean(winRope),
                  "winRight"=mean(winRight) )
  
  return (results)
  
}