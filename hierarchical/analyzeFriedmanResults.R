analyzeFriedmanResults <- function (friedmanType=1) {
  
  actualFileName <- paste('actualAccFriedman',friedmanType,'.csv',sep = '')
  actualAccFriedman <- read.csv(actualFileName)
  
  cvalFileName <- paste('cvalAccFriedman',friedmanType,'.csv',sep = '')
  cvalAccFriedman <- read.csv(cvalFileName)
  
  
  
  #mean estimation error on the whole family, averaged over the repetitions
  #check how many repetitions: how many times the first dset has been run
  idx <- 1 #which setting we use as referneces
  repetitions <- sum(
    cvalAccFriedman$redundantFeats == actualAccFriedman$redundantFeats[1] &
      cvalAccFriedman$sampleSize == actualAccFriedman$sampleSize[1] &
      cvalAccFriedman$friedmanSd == actualAccFriedman$friedmanSd[1])
  #proportion of times hier leads lower error than MLE (out of all the dsets of the family), averaged over repetitions
  
  maeHier <- vector(length = repetitions)
  maeMle <- vector(length = repetitions)
  proportionHierBeatsMLe <- vector(length = repetitions)
  
  counter <- 1
  howManySettings <- dim(actualAccFriedman)[1]
  for (currentRep in 1:repetitions){
    #those two variables will store the Mae (aggregated over all data sets for the current iteration)
    currentMaeHier <- vector (length = howManySettings)
    currentMaeMle <- vector  (length = howManySettings)
    #now loop over the data sets (each data set is also a setting)
    for (currentSetting in 1:howManySettings){
      currentMleEstim <- cvalAccFriedman$mleDiffLdaCart[counter]
      currentHierEstim <- cvalAccFriedman$hierDiffLdaCart[counter]  
      idx = which (actualAccFriedman$redundantFeats == cvalAccFriedman$redundantFeats[counter] &
                     actualAccFriedman$sampleSize == cvalAccFriedman$sampleSize[counter] &
                     actualAccFriedman$friedmanSd == cvalAccFriedman$friedmanSd[counter] 
      )
      actualDifference <- actualAccFriedman$ldaAccuracy[idx] - actualAccFriedman$cartAccuracy[idx] 
      currentMaeHier[currentSetting] <-  abs (currentHierEstim - actualDifference)
      currentMaeMle[currentSetting] <-  abs (currentMleEstim - actualDifference)
      counter <- counter + 1
    }
    maeHier[currentRep] <- mean(currentMaeHier)
    maeMle[currentRep] <- mean(currentMaeMle)
    proportionHierBeatsMLe[currentRep] <- mean (currentMaeHier < currentMaeMle)
  }
  #at this points we save the result to file
  csvFilename <- paste('csvResults/comparisonMleHierFriedman',friedmanType,".csv",sep='')
  results <- data.frame(maeHier, maeMle, proportionHierBeatsMLe)
  write.matrix(results,file=csvFilename, sep=",")
}