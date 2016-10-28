#store multiple useful functions: getActuals, KL

getActuals <- function (){
  
  actualFileName <- 'csvResults/actualAccFriedman1.csv'
  actualAccFriedman1 <- read.csv(actualFileName)
  actualFileName <- 'csvResults/actualAccFriedman2.csv'
  actualAccFriedman2 <- read.csv(actualFileName)
  actualFileName <- 'csvResults/actualAccFriedman3.csv'
  actualAccFriedman3 <- read.csv(actualFileName)
  actualAccFriedman <- rbind (actualAccFriedman1, actualAccFriedman2, actualAccFriedman3)
  ropeMin <- - 0.01
  ropeMax <- 0.01
  
  diff <- actualAccFriedman$ldaAccuracy - actualAccFriedman$cartAccuracy
  actuals <- list (delta0 = mean(diff))
  actuals$pDeltaRight <- length(which(diff>ropeMax)) / length(diff)
  actuals$pDeltaLeft <- length(which(diff<ropeMin)) / length(diff)
  actuals$pDeltaRope <- length(which(diff>ropeMin & diff<ropeMax )) / length(diff)
  actuals$eachSetting <- actualAccFriedman
  
  return (actuals)
  
}

#computes the KL between the actual distribution of the deltaI
#in the population and the distribution estimated by the model, sampling
#from the posterior. We consider rope, left, right
KL <- function (fittedModel,actuals, sampling='student'){
  ropeMin <- - 0.01
  ropeMax <- 0.01
  #here we sample the Deltas from the posterior
  postSamples <- length(fittedModel$stanResults$delta0)
  sampledDeltas <- vector (length = postSamples)
  
  std <- fittedModel$stanResults$std0
  mu  <- fittedModel$stanResults$delta0
  if (sampling=='student'){
    nu  <- fittedModel$stanResults$nu
    for (r in 1:postSamples){
      sampledDeltas[r] <- 1000
      while (abs(sampledDeltas[r]) > 1) {
        sampledDeltas[r]   <- (rt (1, nu[r]) * std[r] + mu[r]) * fittedModel$stdX
      }
    } 
  }
  if (sampling=='gaussian'){
    for (r in 1:postSamples){
      sampledDeltas[r] <- 1000
      while (abs(sampledDeltas[r]) > 1) {
        sampledDeltas[r]   <- rnorm (1, sd= std[r], mean= mu[r]) * fittedModel$stdX
      }
    }
  }
  pSampleRight <- length(which(sampledDeltas>ropeMax)) / postSamples
  pSampleLeft <- length(which(sampledDeltas<ropeMin)) / postSamples
  pSampleRope <- length(which(sampledDeltas>ropeMin & sampledDeltas<ropeMax )) / postSamples
  
  actualP <- c(actuals$pDeltaRight,actuals$pDeltaRope,actuals$pDeltaLeft)
  sampleP <- c(pSampleRight, pSampleRope, pSampleLeft)
  eps <- 0.0001
  sampleP[sampleP<eps] <- eps
  actualP[actualP<eps] <- eps
  klDiv <- sum( actualP * log (actualP/sampleP) )
  return (klDiv)
  
}

getActualTrainDeltaI <- function (actuals, trainSettings){
  for (i in 1:nrow(trainSettings)){
    idx = which ( 
      actuals$eachSetting$redundantFeats == trainSettings$redundantFeats[i] &
        actuals$eachSetting$sampleSize == trainSettings$sampleSize[i] &
        actuals$eachSetting$friedmanSd == trainSettings$friedmanSd[i])
  }
  actualDiff <- actuals$eachSetting$ldaAccuracy[i] - actuals$eachSetting$cartAccuracy[i] 
}


plotPosterior <- function  (hierPosterior, hierPosteriorSens, class1, class2) {
  
  postSamples <- length(hierPosterior$stanResults$delta0)
  sampleDelta <- vector (length = postSamples)
  std <- hierPosterior$stanResults$std0
  mu  <- hierPosterior$stanResults$delta0
  nu  <- hierPosterior$stanResults$nu
  redone <- -1 * postSamples
  redoneSens <- -1 * postSamples
  deltaShr <- hierPosterior$delta_each_dset
  for (r in 1:postSamples){
    sampleDelta[r] <- 1000
    while (abs(sampleDelta[r]) > 1) {
      sampleDelta[r]   <- (rt (1, nu[r]) * std[r] + mu[r]) * hierPosterior$stdX
      redone <- redone + 1
    }
  }
  
  sampleDeltaSens <- vector (length = postSamples)
  std <- hierPosteriorSens$stanResults$std0
  mu  <- hierPosteriorSens$stanResults$delta0
  nu  <- hierPosteriorSens$stanResults$nu
  
  for (r in 1:postSamples){
    sampleDeltaSens[r] <- 1000
    while (abs(sampleDeltaSens[r]) > 1) {
      sampleDeltaSens[r]   <- (rt (1, nu[r]) * std[r] + mu[r]) * hierPosteriorSens$stdX
      redoneSens <- redoneSens + 1
    }
  }
  
  d1 <- density(deltaShr)
  d2 <- density(sampleDelta)
  d3<- density(sampleDeltaSens)
  
  filename <- paste('plotPost',class1,class2,'.pdf',sep = '')
  plot(d1, col=1, xlim = c(-.1,.1))
  lines(d2,col=2)
  lines(d3,col=3)
  legend(0.02,30,legend=c('shr','hier','hierSens'), lty=c(1,1,1),col=c(1,2,3))
  dev.off()
}