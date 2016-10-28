genFriedmanSettings <- function (friedmanType=1) {
  
  #generate epxerimental settings for Friedman family:
  #levels of sample size, std dev, redundant features and the discretization threshold.
  #settings to be expanded later
  source('generateFriedmanData.R')
  redundantFeats <- c (0,20)
  sampleSize <- c(30, 100, 1000)
  
  if (friedmanType==1) {
    friedmanSd <- c(0.5, 1, 2) 
  }
  else  if (friedmanType==2) {
    friedmanSd <- c(62.5, 125, 250) 
  }
  else  if (friedmanType==3) {
    friedmanSd <- c(0.05, 0.1, 0.2) 
  }
  
  
  #data frame containing all the experimental setups
  #settings and testSettings differ as for the sampleSize, which in trainSetginds varies while in testSettings is always 1000.
  settings <- expand.grid(redundantFeats=redundantFeats,sampleSize=sampleSize,friedmanSd=friedmanSd)
  # settings$friedmanType <- rep(friedmanType, length(settings$redundantFeats))
  #now we need to accurately estimate the discr threshold.
  #to this end we generate 10000 samples .The threshold is equal for
  #all settings of the same Friedman function, as it is not affected by sample size,
  #redundant Feats or noise, as it is the median of y.
  
  # freezing the seed is necessary to produce the same threshold for cross-validation
  # and actualFriedmanAccuracy experiments 
  
  set.seed(42)
  tmpSettings <- settings[1,]
  tmpSettings$sampleSize <- 10000
  data <- generateFriedmanData(friedmanType,tmpSettings)
  threshold <- median (data$y)
  settings$threshold <- rep(threshold,dim(settings)[1])
  return (settings)
}