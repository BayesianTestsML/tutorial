friedman1Simul <- function(){
  #implements cross-validation adopting the classifiers from the caret package
  library(mlbench)
  library(caret)
  data(PimaIndiansDiabetes)
  
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