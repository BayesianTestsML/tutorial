actualFriedmanAccuracy <- function(friedmanType=1, reps=500) {
  #estimates actual accuracy of real classifiers on data sets of the Friedman family
  #friedmanType is the friedman function (1,2,3) while reps is the number of repetitions asked for the assessment.
  library(mlbench)
  library(caret)
  source('generateFriedmanData.R')
  source('genFriedmanSettings.R')
  set.seed(42)
  
  
  settings <-  genFriedmanSettings(friedmanType)
  testSettings <- settings
  #every test set will have size 5000
  testSettings$sampleSize <- rep(5000,length(testSettings$sampleSize)) 
  settings$ldaAccuracy <- vector (length = length(testSettings$sampleSize))
  settings$cartAccuracy <- vector (length = length(testSettings$sampleSize))
  settings$knnAccuracy <- vector (length = length(testSettings$sampleSize))
  
  
  #making sure that no cross-validation is performed while fitting the models.
  control <- trainControl(method="none") 
  
  
  for (currentSetting in 1:dim(settings)[1]){
    
    cat('setting: ',currentSetting,'/',dim(settings)[1],'\n')
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
  csvFilename <- paste('csvResults/actualAccFriedman',friedmanType,".csv",sep='')
  write.matrix(settings,file=csvFilename, sep=",")
  
  
  
}