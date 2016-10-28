analyzeFriedmanResults <- function (friedmanType=1) {
  
  #the case in which individual families have to be analyzed
  if (friedmanType < 4){
    actualFileName <- paste('csvResults/actualAccFriedman',friedmanType,'.csv',sep = '')
    actualAccFriedman <- read.csv(actualFileName)
    
    cvalFileName <- paste('csvResults/cvalAccFriedman',friedmanType,'.csv',sep = '')
    cvalAccFriedman <- read.csv(cvalFileName)
  }
  
  #the case in which the three families are jointly experimented
  #following code is a bit rough but works
  if (friedmanType== 4){
    actualFileName <- 'csvResults/actualAccFriedman1.csv'
    actualAccFriedman1 <- read.csv(actualFileName)
    actualFileName <- 'csvResults/actualAccFriedman2.csv'
    actualAccFriedman2 <- read.csv(actualFileName)
    actualFileName <- 'csvResults/actualAccFriedman3.csv'
    actualAccFriedman3 <- read.csv(actualFileName)
    actualAccFriedman <- rbind (actualAccFriedman1, actualAccFriedman2, actualAccFriedman3)
    
    cvalFileName <- paste('csvResults/cvalAccFriedman',friedmanType,'.csv',sep = '')
    cvalAccFriedman <- read.csv(cvalFileName)
  }
  
  #mean estimation error on the whole family, averaged over the repetitions
  #check how many repetitions: how many times the first dset has been run
  idx <- 1 #which setting we use as referneces
  repetitions <- sum(
    cvalAccFriedman$redundantFeats == actualAccFriedman$redundantFeats[1] &
      cvalAccFriedman$sampleSize == actualAccFriedman$sampleSize[1] &
      cvalAccFriedman$friedmanSd == actualAccFriedman$friedmanSd[1])
  #proportion of times hier leads lower error than MLE (out of all the dsets of the family), averaged over repetitions
  
  maeHier <- vector(length = repetitions)
  maeHierGauss <- vector(length = repetitions)
  maeHierKru <- vector(length = repetitions)
  maeHierJua <- vector(length = repetitions)
  maeMle <- vector(length = repetitions)
  proportionDsetsHierBeatsMLe <- vector(length = repetitions)
  HierBeatsMLeJointMae <- vector(length = repetitions)
  HierBeatsKruJointMae <- vector(length = repetitions)
  KruBeatsGaussJointMae <- vector(length = repetitions)
  KruBeatsJuaJointMae <- vector(length = repetitions)
  
  counter <- 1
  howManySettings <- dim(actualAccFriedman)[1]
  for (currentRep in 1:repetitions){
    #those two variables will store the Mae (aggregated over all data sets for the current iteration)
    currentMaeHier <- vector (length = howManySettings)
    currentMaeMle <- vector  (length = howManySettings)
    currentMaeGauss <- vector  (length = howManySettings)
    currentMaeKru <- vector  (length = howManySettings)
    currentMaeJua <- vector  (length = howManySettings)
    
    #now loop over the data sets (each data set is also a setting)
    for (currentSetting in 1:howManySettings){
      currentMleEstim <- cvalAccFriedman$mleDiffLdaCart[counter]
      currentHierEstim <- cvalAccFriedman$hierDiffLdaCart[counter]
      currentGaussEstim <- cvalAccFriedman$gaussDiffLdaCart[counter]
      currentKruEstim <- cvalAccFriedman$kruDiffLdaCart[counter]
      currentJuaEstim <- cvalAccFriedman$juaDiffLdaCart[counter]
      
      idx = which (actualAccFriedman$redundantFeats == cvalAccFriedman$redundantFeats[counter] &
                     actualAccFriedman$sampleSize == cvalAccFriedman$sampleSize[counter] &
                     actualAccFriedman$friedmanSd == cvalAccFriedman$friedmanSd[counter] 
      )
      actualDifference <- actualAccFriedman$ldaAccuracy[idx] - actualAccFriedman$cartAccuracy[idx] 
      currentMaeHier[currentSetting] <-  abs (currentHierEstim - actualDifference)
      currentMaeMle[currentSetting] <-  abs (currentMleEstim - actualDifference)
      currentMaeGauss[currentSetting] <-  abs (currentGaussEstim - actualDifference)
      currentMaeKru[currentSetting] <-  abs (currentKruEstim - actualDifference)
      currentMaeJua[currentSetting] <-  abs (currentJuaEstim - actualDifference)
      
      counter <- counter + 1
    }
    maeHier[currentRep] <- mean(currentMaeHier)
    maeMle[currentRep] <- mean(currentMaeMle)
    maeHierKru[currentRep] <- mean(currentMaeKru)
    maeHierJua[currentRep] <- mean(currentMaeJua)
    maeHierGauss[currentRep] <- mean(currentMaeGauss)
    proportionDsetsHierBeatsMLe[currentRep] <- mean (currentMaeHier < currentMaeMle)
    HierBeatsMLeJointMae[currentRep] <- ifelse ( maeHier[currentRep] < maeMle[currentRep], 1, 0)
    HierBeatsKruJointMae[currentRep] <- ifelse ( maeHier[currentRep] < maeHierKru[currentRep], 1, 0)
    KruBeatsGaussJointMae[currentRep] <- ifelse ( maeHierKru[currentRep] < maeHierGauss[currentRep], 1, 0)
    KruBeatsJuaJointMae[currentRep] <- ifelse ( maeHierKru[currentRep] < maeHierJua[currentRep], 1, 0)
  }
  #at this points we save the result to file
  csvFilename <- paste('csvResults/comparisonMleHierFriedman',friedmanType,".csv",sep='')
  dsetResults <- data.frame(maeMle, maeHierGauss, maeHier, maeHierKru, maeHierJua, HierBeatsMLeJointMae, HierBeatsKruJointMae, 
                            KruBeatsGaussJointMae,KruBeatsJuaJointMae)
  write.matrix(dsetResults,file=csvFilename, sep=",")
  
  results <- list ('dsetResults'=dsetResults, 'HierBeatsMLeJointMae' = HierBeatsMLeJointMae)
  return (results)
}