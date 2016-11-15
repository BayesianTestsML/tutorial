#store multiple useful functions: getActuals, KL

#computes MSEshr and MSEmle for simulation with hierarchical models, 
#assuming data to be laready loaded
computeMSE <- function(){
  reps <- length(hierModels)
  mseShr <- vector(length = reps)
  mseMle <- vector(length = reps)
  for (i in 1:reps){
    actual <- trainData[[i]]$actualTrainDeltaI
    mle <- trainData[[i]]$currentMleDiffLdaCart
    shrinkage <- hierModels[[i]]$delta_each_dset
    mseMle[i] <- mean ( (actual - mle)^2 )
    mseShr[i] <- mean ( (actual - shrinkage)^2 )
  }
  df <- data.frame('mseMle'=mseMle,'mseShr'=mseShr)
  write.table(df,sep=',',col.names = TRUE, row.names = FALSE, file = 'mse.tex')
}


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
  actualDiff <- vector (length = nrow(trainSettings))
  for (i in 1:nrow(trainSettings)){
    idx = which ( 
      actuals$eachSetting$redundantFeats == trainSettings$redundantFeats[i] &
        actuals$eachSetting$sampleSize == trainSettings$sampleSize[i] &
        actuals$eachSetting$friedmanSd == trainSettings$friedmanSd[i])
    
    actualDiff[i] <- actuals$eachSetting$ldaAccuracy[idx] - actuals$eachSetting$cartAccuracy[idx] 
  }
  return(actualDiff)
}


plotPosterior <- function  (hierPosterior, hierPosteriorSens, class1=99, class2=99) {
  
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
  
  d1 <- density(deltaShr, from = 2*min(deltaShr), to=2*max(deltaShr), n = 512 )
  d2 <- density(sampleDelta, from = 2*min(deltaShr), to=2*max(deltaShr), n = 512 )
  d3<- density(sampleDeltaSens, from = 2*min(deltaShr), to=2*max(deltaShr), n = 512  )
  
  filename <- paste('plotPost',class1,class2,'.pdf',sep = '')
  plot(d1, col=1, xlim = c(-.1,.1))
  lines(d2,col=2)
  lines(d3,col=3)
  legend(0.02,10,legend=c('shr','hier','hierSens'), lty=c(1,1,1),col=c(1,2,3))
  # dev.off()
}

#this function parses the list of hierarchical models and return relevant facts about the estimated probabilities.
postEstimatedDelta <- function (hierModels){
  #gets a list of hierarchical models and estimates via sampling the probability of the next delta being left, rope and right.
  reps <- length(hierModels)
  EstimPLeft <- vector (length = reps)
  EstimPRope <- vector (length = reps)
  EstimPRight <- vector (length = reps)
  EstimPLeftBias <- vector (length = reps)
  EstimPRopeBias <- vector (length = reps)
  EstimPRightBias <- vector (length = reps)
  
  postSamples <- length(hierModels[[1]]$stanResults$delta0)
  sampledDeltas <- vector (length = postSamples)
  
  for (i in 1:reps){
    fittedModel <- hierModels[[i]]
    std <- fittedModel$stanResults$std0
    mu  <- fittedModel$stanResults$delta0
    nu  <- fittedModel$stanResults$nu
    for (r in 1:postSamples){
      sampledDeltas[r] <- 1000
      while (abs(sampledDeltas[r]) > 1) {
        sampledDeltas[r]   <- (rt (1, nu[r]) * std[r] + mu[r]) * fittedModel$stdX
      }
    } 
    #probabilities of the sampled deltas to be left rigth rope
    EstimPLeft[i] <- mean (sampledDeltas < -0.01)
    EstimPRight[i] <- mean (sampledDeltas > 0.01)
    EstimPRope[i] <- mean (sampledDeltas < 0.01 & sampledDeltas > -0.01)
    EstimPLeftBias[i] <- hierModels[[i]]$nextDelta$left
    EstimPRopeBias[i] <- hierModels[[i]]$nextDelta$rope
    EstimPRightBias[i] <- hierModels[[i]]$nextDelta$right
    
  }
  actuals <- getActuals()
  #probability of delta belonging to left, rope and right
  estimatedDelta <- data.frame('EstimPLeft'=EstimPLeft, 'EstimPRight'= EstimPRight, 'EstimPRope' = EstimPRope)
  #probability that we return as inference
  estimatedDeltaBias <- data.frame('pLeft'=EstimPLeftBias, 'pRight'= EstimPRightBias, 'pRope' = EstimPRopeBias)
  
  return (list('actuals'=actuals, 'estimatedDelta'=estimatedDelta, 'estimatedDeltaBias'=estimatedDeltaBias))
}




#computes the Kl div between the posterior and the shrinkage estimates
KLPostShrinkage <- function (fittedModel,shrinkageEstimates){
  library('flexmix')
  #from the posterior t or the posterior gaussian 
  postSamples <- length(fittedModel$stanResults$delta0)
  sampledDeltas <- vector (length = postSamples)
  
  std <- fittedModel$stanResults$std0
  mu  <- fittedModel$stanResults$delta0
  
  
  #if the model is Student
  if (!is.null(fittedModel$nu)){
    nu  <- fittedModel$stanResults$nu
    for (r in 1:postSamples){
      sampledDeltas[r] <- 1000
      while (abs(sampledDeltas[r]) > 1) {
        sampledDeltas[r]   <- (rt (1, nu[r]) * std[r] + mu[r]) * fittedModel$stdX
      }
    } 
  }
  #otherwise gaussian sampling
  else{
    for (r in 1:postSamples){
      sampledDeltas[r] <- 1000
      while (abs(sampledDeltas[r]) > 1) {
        sampledDeltas[r]   <- rnorm (1, sd= std[r], mean= mu[r]) * fittedModel$stdX
      }
    }
  }
  d1 <- density(shrinkageEstimates, from = 2*min(shrinkageEstimates), to=2*max(shrinkageEstimates), n = 512 )
  d2 <- density(sampledDeltas, from = 2*min(shrinkageEstimates), to=2*max(shrinkageEstimates), n = 512 )
  estimatedKL <- KLdiv(cbind(shrink=d1$y,fitted=d2$y))
  return(estimatedKL)
}



#produces boxplot and scatter plots for
#comparing MLE and shrunken estimates
plotShrinkageMle <- function (){
  library(tikzDevice)
  load("~/Documents/devel/tutorialML/hierarchical/Rdata/cvalFriedmanPredictive.Rdata")
  tikz("bplotHierShr.tex", width= 6.5, height=4.5)
  boxplot( hierModels[[2]]$delta_each_dset, trainData[[2]]$currentMleDiffLdaCart, 
           xlab=c('hier','mle'), ylab="Estimate of $delta_i$")
  dev.off()
  
  tikz("scatterHierActual.tex", width= 6.5, height=4.5)
  plot(dset$averageTime ~ dset$AthleteName, xlab="", ylab="Execution time")
  dev.off()
  
  tikz("scatterMleActual.tex", width= 6.5, height=4.5)
  plot(dset$averageTime ~ dset$AthleteName, xlab="", ylab="Execution time")
  dev.off()
  
}

plotPosteriorGGplot2 <- function  (hierPosterior, hierPosteriorKru, hierPosteriorJua, suffix) {
  # library(tikzDevice)
  library(ggplot2)
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
  
  sampleDeltaKru <- vector (length = postSamples)
  std <- hierPosteriorKru$stanResults$std0
  mu  <- hierPosteriorKru$stanResults$delta0
  nu  <- hierPosteriorKru$stanResults$nu
  redoneKru <- -1 * postSamples
  for (r in 1:postSamples){
    sampleDeltaKru[r] <- 1000
    while (abs(sampleDeltaKru[r]) > 1) {
      sampleDeltaKru[r]   <- (rt (1, nu[r]) * std[r] + mu[r]) * hierPosteriorKru$stdX
      redoneKru <- redoneKru + 1
    }
  }
  
  sampleDeltaJua <- vector (length = postSamples)
  std <- hierPosteriorJua$stanResults$std0
  mu  <- hierPosteriorJua$stanResults$delta0
  nu  <- hierPosteriorJua$stanResults$nu
  for (r in 1:postSamples){
    sampleDeltaJua[r] <- 1000
    while (abs(sampleDeltaJua[r]) > 1) {
      sampleDeltaJua[r]   <- (rt (1, nu[r]) * std[r] + mu[r]) * hierPosteriorJua$stdX
    }
  }
  
  limit<-max(abs(min(deltaShr)),abs(max(deltaShr)))
  d1 <- density(deltaShr, from = -2*limit, to=2*limit, n = 512 )
  d2 <- density(sampleDelta, from = -2*limit, to=2*limit, n = 512 )
  d3 <- density(sampleDeltaKru, from = -2*limit, to=2*limit, n = 512  )
  d4 <- density(sampleDeltaJua, from = -2*limit, to=2*limit, n = 512  )
  
  #save the density to file
  stopifnot(mean(d1$x==d2$x)==1)
  stopifnot(mean(d1$x==d3$x)==1)
  df<- data.frame(x=d1$x,shr=d1$y,hier=d2$y,kru=d3$y,jua=d4$y) 
  fileName <- paste('csvResults/density',suffix,'.tex',sep='')
  write.table(df,file=fileName,sep=',',row.names = FALSE, quote = FALSE)
  
  
  filename <- paste('plotPostHierKru',suffix,'.pdf',sep = '')
  pdf(file=filename)
  plot(d1, col='black', xlim = c(-.1,.1), main = '' )
  lines(d2,col='blue')
  lines(d3,col='green')
  # lines(d3,col='red')
  legend('topright',c('shr','hier','jua'),lty=c(1,1,1),
         col=c('black','blue', 'green'))
  #    legend('topright',c('shr','hier','kru','jua'),lty=c(1,1,1),
  #           col=c('black','blue', 'green','red'))
  
  
  
  dev.off()
  
  # deltaShrNA<- c(deltaShr,rep(NA,length(sampleDelta) - length(deltaShr)))
  # m <- data.frame(hier=sampleDelta,kru=sampleDeltaKru,shrinkage=deltaShrNA)
  # dfs <- stack(m)
  # outFile <- paste('densityPlot',suffix,'.pdf',sep='')
  # pdf(outFile)
  # tikz(outFile, width= 6.5, height=4.5)
  # p <- ggplot(dfs, aes(x=values)) + geom_density(aes(group=ind, colour=ind, fill=ind), alpha=0.3) + xlim(-0.1, 0.1)
  
  # dev.off()
}
