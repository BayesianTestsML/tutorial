friedmanFamily <- function() {
  #implements cross-validation and estimation of the real accuracies
  #on the family of Friedman data sets, adopting the classifiers from the caret package
  library(mlbench)
  library(caret)
  
  #settings to be expanded later
  redundantFeats <- c (0,10,100)
  sampleSizes <- c(30, 300)
  noise <- c(1) 
  classifier <- c('lda','cart')
  cvFolds <- 10
  cvRuns <- 10
  
  
  #function called later in the main
  #simulate 100 times with a large test and average the accuracies
  getActualAcc <- function (classifierIdx,featIdx,sampleIdx,noiseIdx,seed){
    #generate the data, fixing the seed so that they are equal for each classifier
  }
  
  
  
  
  #actualAcc is matrix, whose rows are the different classifiers and whose columns are the different settings
  actualAcc <- matrix (ncol = length(classifier), nrow = length(redundantFeats) * length(sampleSizes) * length(noise))
  cvalAcc <- array (dim=c(cvFolds * cvRuns, dim(actualAcc)))
  seed <- 42 #to generate the same data
  
  for (classifierIdx in 1:length(classifier)){
    currentClassifier <- classifier (classifierIdx)
    for (featIdx in 1:length(redundantFeats)) {
      currentFeats <- redundantFeats ( featIdx )
      for (noiseIdx in 1:length(noise)){
        currentNoise <- noise[noiseIdx]
        for (sampleIdx in 1:length(sampleSizes)){
          currentSample <- sampleSizes[sampleIdx]
          actualAcc [classifierIdx,featIdx,sampleIdx,noiseIdx] <- 
            getActualAcc(classifier[classifierIdx],redundantFeats[featIdx],sampleSizes[sampleIdx],noise[noiseIdx],seed)
        }
      }
    }
  } 
  
  
  
  
  control <- trainControl(method="repeatedcv", number=10, repeats=3)
  currentSeed<-7 
  
  data <- mlbench.friedman1(10000, sd=10)
  summary(data$y)
  
  cut.vec <- cut(data$y, breaks=quantile(data$y, 0:3/3), include.lowest=T, right=FALSE) 
  data$class <- cut.vec
  
  
  #we set the seed before running each classifier in order to generate paired folds.cart[1]
  set.seed(currentSeed)
  fit.cart <- train(data$x, data$class, method="rpart1SE", trControl=control,tuneGrid=NULL, tuneLength = 1)
  cartAccuracy <- fit.cart$resample$Accuracy
  
  
  #LDA
  set.seed(currentSeed)
  fit.lda <- train(data$x, data$class, method="lda", trControl=control, tuneLength = 1)
  ldaAccuracy <- fit.lda$resample$Accuracy
  
  #SVM
  set.seed(currentSeed)
  fit.svm <- train(diabetes~., data=PimaIndiansDiabetes, method="svmRadial", trControl=control, tuneLength = 1)
  svmAccuracy <- fit.svm$resample$Accuracy
  
  
  # kNN
  set.seed(currentSeed)
  fit.knn <- train(diabetes~., data=PimaIndiansDiabetes, method="knn", trControl=control, tuneLength = 1)
  knnAccuracy <- fit.knn$resample$Accuracy
  
  # Random Forest
  set.seed(currentSeed)
  fit.rf <- train(diabetes~., data=PimaIndiansDiabetes, method="rf", trControl=control, tuneLength = 1)
  rfAccuracy <- fit.rf$resample$Accuracy
  
}