actualFriedmanAccuracy <- function(friedmanType=1, reps=1000) {
  #estimates actual accuracy of real classifiers on data sets of the Friedman family
  #friedmanType is the friedman function (1,2,3) while reps is the number of repetitions asked for the assessment.
  library(mlbench)
  library(caret)
  source('generateFriedmanData.R')
  
  #settings to be expanded later
  redundantFeats <- c (0,20)
  sampleSize <- c(30, 100, 1000)
  
  
  set.seed(42)
  
  
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
  testSettings <- settings
  testSettings$sampleSize <- rep(5000,length(testSettings$sampleSize)) 
  #the threshold will be dinamically filled
  settings$threshold <- vector (length = length(testSettings$sampleSize)) 
  settings$ldaAccuracy <- vector (length = length(testSettings$sampleSize))
  settings$cartAccuracy <- vector (length = length(testSettings$sampleSize))
  settings$knnAccuracy <- vector (length = length(testSettings$sampleSize))
  
  classifier <- c('lda','cart','knn')
  
  #making sure that no cross-validation is performed while fitting the models.
  control <- trainControl(method="none") 
  
  
  for (currentSetting in 1:dim(settings)[1]){
    
    cat('setting: ',currentSetting,'/',dim(settings)[1],'\n')
    testData <- generateFriedmanData(friedmanType,testSettings[currentSetting,])
    settings$threshold[currentSetting] <- median (testData$y)
    testSettings$threshold[currentSetting] <- median (testData$y)
    ldaAccuracy <- vector(length = reps)
    cartAccuracy  <- vector(length = reps)
    knnAccuracy  <- vector(length = reps)
    
    for (currentRep in 1:reps){
      cat(currentRep,'\n')
      #the following functions should generate the data according to the given settings and discretize according to
      #the supplied threshold
      trainingData <- generateFriedmanData(friedmanType,settings[currentSetting,])
      testData <- generateFriedmanData(friedmanType,testSettings[currentSetting,])
      
      
      #fit lda
      fit.lda <- train(trainingData$x, trainingData$class, method="lda", trControl=control, tuneLength = 1)
      ldaPredictions <- predict(fit.lda,testData$x)
      ldaAccuracy[currentRep] <- mean (ldaPredictions == testData$class)
      
      fit.cart <- train(trainingData$x, trainingData$class, method="rpart1SE", trControl=control, tuneLength = 1)
      cartPredictions <- predict (fit.cart,testData$x)
      cartAccuracy[currentRep] <- mean (cartPredictions == testData$class)
      
      fit.knn <- train(trainingData$x, trainingData$class, method="knn", trControl=control, tuneLength = 1)
      knnPredictions <- predict (fit.knn,testData$x)
      knnAccuracy[currentRep] <- mean (knnPredictions == testData$class)      
    }
    #at this point all the simulation for the given settings have been done.
    #let's store the results
    settings$ldaAccuracy[currentSetting] <- mean (ldaAccuracy)
    settings$cartAccuracy[currentSetting] <- mean (cartAccuracy)
    settings$knnAccuracy[currentSetting] <- mean (knnAccuracy)
  }
  settings$repetitions <- rep(reps,length(testSettings$sampleSize))
  settings$FriedmanType <- rep(friedmanType,length(testSettings$sampleSize))
  #at this points we save the result to file
  csvFilename <- paste('actualAccFriedman',friedmanType,".csv",sep='')
  write.matrix(settings,file=csvFilename, sep=",")
  
  
  
}