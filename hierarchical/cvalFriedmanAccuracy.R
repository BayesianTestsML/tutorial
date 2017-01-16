cvalFriedmanAccuracy <- function(friedmanType=1, reps=250) {
  #estimates  accuracy of  classifiers on data sets of the Friedman family via cross-validation 
  #(10 runs of 10 folds)
  #friedmanType is the friedman function (1,2,3) while reps is the number of repetitions
  #of the process data generation/cross-validation.
  #the family contains right onw 18 data sets.
  library(mlbench)
  library(caret)
  source('generateFriedmanData.R')
  source('genFriedmanSettings.R')
  source('hierarchical_test.R')
  
  
  if (friedmanType < 4){
    settings <- genFriedmanSettings(friedmanType)
    friedmanTypeVec <- rep (friedmanType,nrow(settings))
  }
  
  if (friedmanType == 4){
    settings <- rbind(genFriedmanSettings(1),genFriedmanSettings(2),genFriedmanSettings(3)) 
    friedmanTypeVec <- 
      cbind ( t(rep(1,nrow(genFriedmanSettings(1)))), t(rep(2,nrow(genFriedmanSettings(2)))),
              t(rep(3,nrow(genFriedmanSettings(3))))  )
    #friedman Type is a vector, containing the friedman Type of each setting
  }
  
  totExperiments <- dim(settings)[1] * reps
  mleDiffLdaCart <-  vector (length = totExperiments)
  hierDiffLdaCart <- vector (length = totExperiments)
  gaussDiffLdaCart <- vector (length = totExperiments)
  kruDiffLdaCart <- vector (length = totExperiments)
  juaDiffLdaCart <- vector (length = totExperiments)
  #knnAccuracy <-  vector (length = totExperiments)
  redundantFeats <- vector (length = totExperiments)
  sampleSize <- vector (length = totExperiments)
  friedmanSd <- vector (length = totExperiments)
  
  #10 runs of 10-folds cross-validation, hard-coded
  nRuns <- 10
  nFolds <- 10
  control <- trainControl(method = "repeatedcv", number=nFolds, repeats=nRuns)
  
  for (currentRep in 1:reps){
    cat('Repetition:', currentRep,'\n')
    crossValResults <-  matrix(nrow=dim(settings)[1], ncol = nRuns * nFolds);
    
    for (currentSetting in 1:dim(settings)[1]){
      
      #the following functions should generate the data according to the given setting
      data <- generateFriedmanData(friedmanTypeVec[currentSetting],settings[currentSetting,])
      
      #for simplicity we focus on the difference between  lda and cart
      #we need to set the seed in order to pair the folds
      #we also need the seed to be different between each execution
      currentSeed <- as.numeric(Sys.time())
      set.seed(currentSeed)
      
      fit.lda <- train(data$x, data$class, method="lda", trControl=control, tuneLength = 1)
      set.seed(currentSeed)
      fit.cart <- train(data$x, data$class, method="rpart1SE", trControl=control, tuneLength = 1)
      crossValResults[currentSetting,] <- fit.lda$resample$Accuracy - fit.cart$resample$Accuracy
      
      #we exploit the fact that results vectors are initially filled with FALSE
      #and we save the settings of the experiments (equal for all data sets)
      firstAvailable <- min (which (sampleSize == FALSE))
      redundantFeats[firstAvailable] <- settings[currentSetting,]$redundantFeats
      sampleSize[firstAvailable] <- settings[currentSetting,]$sampleSize
      friedmanSd[firstAvailable] <- settings[currentSetting,]$friedmanSd
      
      
    }
    currentMleDiffLdaCart <- apply (crossValResults, 1, mean )
    
    #at this point all the simulation for the given setting and repetitions have been done.
    stanFileName <- paste ('Stan',friedmanType, sep='')
    stanResults <- hierarchical.test(x = crossValResults, rho = 1/nFolds, sample_file = stanFileName, chains=4)
    currentHierDiffLdaCart <- stanResults$delta_each_dset
    
#     #Gaussian hierarchical model.
#     stanGaussianResults <-  hierarchical.test(x = crossValResults, rho = 1/nFolds, sample_file = stanFileName, chains=4, samplingType = "gaussian")
#     currentHierGaussDiffLdaCart <- stanGaussianResults$delta_each_dset
#     
#     #Kru prior
#     kruResults <-  hierarchical.test(x = crossValResults, rho = 1/nFolds, sample_file = stanFileName, chains=4, samplingType = "studentKruschke")
#     currentHierKruDiffLdaCart <- kruResults$delta_each_dset
#     
#     #Jua prior
#     juaResults <-  hierarchical.test(x = crossValResults, rho = 1/nFolds, sample_file = stanFileName, chains=4, samplingType = "studentJuanez")
#     currentHierJuaDiffLdaCart <- juaResults$delta_each_dset
    
    #we exploit the fact that both vectors are initially filled with FALSE
    firstAvailable <- min (which (mleDiffLdaCart == FALSE))
    mleDiffLdaCart  [ firstAvailable : (firstAvailable + length(currentMleDiffLdaCart) -1 ) ] <- currentMleDiffLdaCart
    hierDiffLdaCart [ firstAvailable : (firstAvailable + length(currentHierDiffLdaCart) -1) ] <- currentHierDiffLdaCart
#     gaussDiffLdaCart [ firstAvailable : (firstAvailable + length(currentHierGaussDiffLdaCart) -1) ] <- currentHierGaussDiffLdaCart
#     kruDiffLdaCart [ firstAvailable : (firstAvailable + length(currentHierKruDiffLdaCart) -1) ] <- currentHierKruDiffLdaCart
#     juaDiffLdaCart [ firstAvailable : (firstAvailable + length(currentHierJuaDiffLdaCart) -1) ] <- currentHierJuaDiffLdaCart
    
  }
  
  #at this points we save the result to file
  csvFilename <- paste('csvResults/cvalAccFriedman',friedmanType,".csv",sep='')
  # results <- data.frame(redundantFeats, sampleSize, friedmanSd, mleDiffLdaCart, hierDiffLdaCart,gaussDiffLdaCart,kruDiffLdaCart,juaDiffLdaCart)
  results <- data.frame(redundantFeats, sampleSize, friedmanSd, mleDiffLdaCart, hierDiffLdaCart)
  write.matrix(results,file=csvFilename, sep=",")
  
}