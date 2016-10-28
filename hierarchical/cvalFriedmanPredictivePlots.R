cvalFriedmanPredictivePlots <- function(reps=1) {
  #run one single runs of cross-validation producing 
  #scatterplots of mle and hier estimates,
  #and comparison of the actual and estimated densities of delta_i
  #estimates  accuracy of  classifiers on data sets of the Friedman family via cross-validation 
  #(10 runs of 10 folds)
  #we jointly run all the three Friedman families. 
  #tracks the MAE of MLE and of different variants of hier models on the delta_i
  #tracks MAE on delta0 and log_loss on P(delta) belonging to left, right, rope 
  library(mlbench)
  library(caret)
  library(tikzDevice)
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
  hierMaeDelta0 <- vector (length = reps)
  hierKLDeltaI <- vector (length = reps)
  hierModels <- list()
  trainData <- list()
  set.seed(42)
  
  
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
    hierModel <- hierarchical.test (x = crossValResults, sample_file = stanFileName,  chains=4, samplingType = "student", alphaBeta = alphaBeta)
    
    stdX <- hierModel$stdX
    hierMaeDelta0[currentRep]<- abs(actuals$delta0 - mean (hierModel$stanResults$delta0 *stdX))
    hierKLDeltaI[currentRep] <- KL(hierModel,actuals)
    
    #track the mae on DeltaI of the dsets analyzed via CV
    actualTrainDeltaI <- getActualTrainDeltaI (actuals, trainSettings)
    mleMaeTrainDeltaI[currentRep] <- mean (abs (actualTrainDeltaI - currentMleDiffLdaCart)       )
    hierMaeTrainDeltaI[currentRep] <- mean (abs (actualTrainDeltaI - hierModel$delta_each_dset)  ) 
    
    #store results for later analysis
    currentTrainList <- list('crossValResults'=crossValResults,'currentMleDiffLdaCart'=currentMleDiffLdaCart,'actualTrainDeltaI'=actualTrainDeltaI)
    trainData <- list (trainData,  currentTrainList, actualTrainDeltaI)
    hierModels <- list( hierModels, hierModel)
    
    df<- data.frame('actual'=actualTrainDeltaI,'MLE'=currentMleDiffLdaCart,
                    'Hier'=hierModel$delta_each_dset)
    fileName <- 'plots/scatterHierMle.csv'
    write.matrix(df,file = fileName)
        
#     yLimit=c(min(currentMleDiffLdaCart),max(currentMleDiffLdaCart))
#     xLimit <- yLimit
#     #boxplot of the session-average results
#      tikz("plots/scatterHierActual.tex", width= 3, height=1.5)
#     plot(actualTrainDeltaI,hierModel$delta_each_dset, 
#          xlab='Actual', ylab= 'Hier estimate',xlim=xLimit, ylim=yLimit)
#     abline(0,1)
#      dev.off()
#     
#     tikz("plots/scatterHierMle.tex", width= 3, height=1.5)
#     plot(actualTrainDeltaI,currentMleDiffLdaCart, 
#          xlab='Actual', ylab= 'MLE estimate', xlim=xLimit, ylim=yLimit)
#     abline(0,1)
#     dev.off()
    
    #save results to file
    tmpFilename <- paste('cvalFriedmanPredictivePlot',currentRep,'.Rdata',sep = ',')
    save (trainData, hierModels,  
          mleMaeTrainDeltaI, hierMaeTrainDeltaI, 
          hierMaeDelta0, hierKLDeltaI, 
          file = tmpFilename) 
    
    #now plot the posterior
    postSamples <- length(hierModel$stanResults$delta0)
    sampleDelta <- vector (length = postSamples)
    std <- hierModel$stanResults$std0
    mu  <- hierModel$stanResults$delta0
    nu  <- hierModel$stanResults$nu
    deltaShr <- hierModel$delta_each_dset
    for (r in 1:postSamples){
      sampleDelta[r] <- 1000
      while (abs(sampleDelta[r]) > 1) {
        sampleDelta[r]   <- (rt (1, nu[r]) * std[r] + mu[r]) * hierModel$stdX
      }
    }
    deltaShr <- hierModel$delta_each_dset
    d1 <- density(deltaShr, from = 2*min(deltaShr), to=2*max(deltaShr), n = 512 )
    d2 <- density(sampleDelta, from = 2*min(deltaShr), to=2*max(deltaShr), n = 512 )
    
    filename <- paste('plots/densityShrHier.pdf')
    pdf(file=filename)
    plot(d1, col=1, xlim = c(-.1,.1))
    lines(d2,col=2)
    legend(0.02,10,legend=c('shr','hier'), lty=c(1,1),col=c(1,2))
    dev.off()
  }
  
  maeTrainDeltaI <- list(
    mleMaeTrainDeltaI=mleMaeTrainDeltaI,
    hierMaeTrainDeltaI=hierMaeTrainDeltaI,
    kruMaeTrainDeltaI=kruMaeTrainDeltaI)
  
  KLDelta <- list(hierKLDeltaI=hierKLDeltaI,kruKLDeltaI=kruKLDeltaI)
  maeDelta0 <- list(hierMaeDelta0=hierMaeDelta0,kruMaeDelta0=kruMaeDelta0)
  #at this points we save the result to file
  filename <- 'Rdata/cvalFriedmanPredictive.Rdata'
  save (trainData, hierModels, kruModels, maeTrainDeltaI,KLDelta,maeDelta0,file = filename)  
}