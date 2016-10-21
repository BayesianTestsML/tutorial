cvalFriedmanPredictive <- function(reps=50) {
  #estimates  accuracy of  classifiers on data sets of the Friedman family via cross-validation 
  #(10 runs of 10 folds)
  #we jointly run all the three Friedman families. 
  #tracks the MAE of MLE and of different variants of hier models on the delta_i
#tracks MAE on delta0 and log_loss on P(delta) belonging to left, right, rope 
library(mlbench)
  library(caret)
  source('generateFriedmanData.R')
  source('genFriedmanSettings.R')
  source('hierarchical_test.R')
  source('selectTrainSettings.R')
  source('Utils.R')
  
  settings <- rbind(genFriedmanSettings(1),genFriedmanSettings(2),genFriedmanSettings(3)) 
  friedmanTypeVec <- 
    cbind ( t(rep(1,nrow(genFriedmanSettings(1)))), t(rep(2,nrow(genFriedmanSettings(2)))),
            t(rep(3,nrow(genFriedmanSettings(3))))  )
  actuals <- getActuals()
  
  totExperiments <- dim(settings)[1] * reps
  mleMaeTrainDeltaI <-  vector (length = reps)
  hierMaeTrainDeltaI <- vector (length = reps)
  gaussMaeTrainDeltaI <- vector (length = reps)
  kruMaeTrainDeltaI <- vector (length = reps)
  hierMaeDelta0 <- vector (length = reps)
  gaussMaeDelta0 <- vector (length = reps)
  kruMaeDelta0 <- vector (length = reps)
  hierKLDeltaI <- vector (length = reps)
  gaussKLDeltaI <- vector (length = reps)
  kruKLDeltaI <- vector (length = reps)
  gaussModels <- list()
  hierModels <- list()
  kruModels <- list()
  trainData <- list()
  
  
  #10 runs of 10-folds cross-validation, hard-coded
  nRuns <- 10
  nFolds <- 10
  control <- trainControl(method = "repeatedcv", number=nFolds, repeats=nRuns)
  
  for (currentRep in 1:reps){
    set.seed(currentRep)
    cat('Repetition:', currentRep,'/',reps,'\n')
    trainIdx <- selectTrainSettings(friedmanTypeVec)
    crossValResults <-  matrix(nrow=length(trainIdx), ncol = nRuns * nFolds)
    trainSettings <- settings[trainIdx,]
    trainFriedmanTypeVec <- friedmanTypeVec[trainIdx] 
    
    
    for (currentSetting in 1:nrow(trainSettings)){
      
      #the following functions should generate the data according to the given setting
      data <- generateFriedmanData(trainFriedmanTypeVec[currentSetting],trainSettings[currentSetting,])
      
      #for simplicity we focus on the difference between  lda and cart
      #we need to set the seed in order to pair the folds
      #we also need the seed to be different between each execution
      set.seed(currentRep)
      fit.lda <- train(data$x, data$class, method="lda", trControl=control, tuneLength = 1)
      set.seed(currentRep)
      fit.cart <- train(data$x, data$class, method="rpart1SE", trControl=control, tuneLength = 1)
      crossValResults[currentSetting,] <- fit.lda$resample$Accuracy - fit.cart$resample$Accuracy
    }
    currentMleDiffLdaCart <- apply (crossValResults, 1, mean )
    
    #at this point all the simulation for the given setting and repetitions have been done.
    stanFileName <- paste ('StanHier', sep='')
    alphaBeta = list('lowerAlpha' =0.5,'upperAlpha'= 3,'lowerBeta' = 0.005,'upperBeta' = 0.05)
    #infer the different models
    hierModel <- hierarchical.test (x = crossValResults, sample_file = stanFileName,  chains=4, samplingType = "student", alphaBeta = alphaBeta)
    gaussModel <-  hierarchical.test(x = crossValResults, sample_file = stanFileName, chains=4, samplingType = "gaussian")
    kruModel <-  hierarchical.test(x = crossValResults, sample_file = stanFileName, chains=4, samplingType = "studentKruschke")
    
    #track the error on Delta0
    stdX <- hierModel$stdX
    stopifnot(hierModel$stdX == gaussModel$stdX)
    stopifnot(hierModel$stdX == kruModel$stdX)
    hierMaeDelta0[currentRep]<- abs(actuals$delta0 - mean (hierModel$stanResults$delta0 *stdX))
    gaussMaeDelta0[currentRep] <- abs(actuals$delta0 - mean (gaussModel$stanResults$delta0 *stdX)) 
    kruMaeDelta0[currentRep] <- abs(actuals$delta0 - mean (kruModel$stanResults$delta0 *stdX))
    
    #track the  KL on prob of DeltaI
    hierKLDeltaI[currentRep] <- KL(hierModel,actuals)
    gaussKLDeltaI[currentRep] <- KL(gaussModel,actuals,sampling='gaussian')
    kruKLDeltaI[currentRep] <- KL(kruModel,actuals)
    
    #track the mae on DeltaI of the dsets analyzed via CV
    actualTrainDeltaI <- getActualTrainDeltaI (actuals, trainSettings)
    mleMaeTrainDeltaI[currentRep] <- mean (abs (actualTrainDeltaI - currentMleDiffLdaCart)       )
    hierMaeTrainDeltaI[currentRep] <- mean (abs (actualTrainDeltaI - hierModel$delta_each_dset)  ) 
    gaussMaeTrainDeltaI[currentRep]     <-  mean (abs (actualTrainDeltaI - gaussModel$delta_each_dset) )
    kruMaeTrainDeltaI[currentRep]       <-  mean (abs (actualTrainDeltaI - kruModel$delta_each_dset)   )
    
    #store results for later analysis
    currentTrainList <- list(crossValResults,currentMleDiffLdaCart,actualTrainDeltaI)
    trainData <- list (trainData,  currentTrainList)
    gaussModels <- list(gaussModels, gaussModel)
    hierModels <- list( hierModels, hierModel)
    kruModels <- list(kruModels,kruModel)
    
    #save tmp results
    tmpFilename <- paste('Rdata/cvalFriedmanPredictiveTMP',currentRep,'.Rdata',sep = '')
    save (trainData, gaussModels, hierModels, kruModels, 
          mleMaeTrainDeltaI, hierMaeTrainDeltaI, gaussMaeTrainDeltaI, kruMaeTrainDeltaI, 
          hierMaeDelta0, gaussMaeDelta0, kruMaeDelta0, hierKLDeltaI, gaussKLDeltaI, 
          kruKLDeltaI, file = tmpFilename)  
    
  }
  
  maeTrainDeltaI <- list(
    mleMaeTrainDeltaI=mleMaeTrainDeltaI,
    hierMaeTrainDeltaI=hierMaeTrainDeltaI,
    gaussMaeTrainDeltaI=gaussMaeTrainDeltaI,
    kruMaeTrainDeltaI=kruMaeTrainDeltaI)
  
  KLDelta <- list(hierKLDeltaI=hierKLDeltaI,gaussKLDeltaI=gaussKLDeltaI,kruKLDeltaI=kruKLDeltaI)
  maeDelta0 <- list(hierMaeDelta0=hierMaeDelta0, gaussMaeDelta0=gaussMaeDelta0, kruMaeDelta0=kruMaeDelta0)
  #at this points we save the result to file
  filename <- 'Rdata/cvalFriedmanPredictive.Rdata'
  save (trainData, gaussModels, hierModels, kruModels, maeTrainDeltaI,KLDelta,maeDelta0,file = filename)  
}