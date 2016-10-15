#compares the predictive estimates yielded by the MLE and  the hierarchical model. 
#It infers the mean difference on the first 10 cross-validation  runs
#(10 x 10 = 100 observations)
#and predicts the mean of the whole  500 runs, which supposedly is the actual difference of accuracy.

predictive_Cval500 <- function() {
  library(MASS)
  library(rstan)
  rstan_options(auto_write = TRUE)
  
  options(mc.cores = parallel::detectCores())
  source ("hierarchical_test.R")
  std_upper_bound <- 1000
  standardization <- 1 
  
  #utiliy function which gets the average score of the given classifier on each data set
  #this is assumed to be the real accuracy
  getAvgScore <- function (currentClassifier,score) {
    avg_score <- vector(length = how_many_dsets, mode = "double");
    for ( ii in 1:how_many_dsets ){
      avg_score[ii] <- mean ( score [as.numeric(classifierId)==currentClassifier & as.numeric(dsetId)==dsetsList[ii]] )
    }
    return (avg_score)
  }
  

  
  #after having loaded the csv, this code build and save the useful results
#   dsetId<- as.factor(resultsWeka500Runs$Key_Dataset)
#   classifierId <- as.factor(resultsWeka500Runs$Key_Scheme)
#   score <- resultsWeka500Runs$Percent_correct
#   foldID <- resultsWeka500Runs$Key_Fold
#   runID <- resultsWeka500Runs$Key_Run
  #the fold ID discriminates between different folds and different runs

  load("weka500runs.Rdata")
  nFolds <- 10 #hard-coded
  rho=1/nFolds #hard-coded, at the moment foldId is not useful
  trainResults <- 100 #how many results to infer MLE and hierarchical
  
  #variables classifierID, dsetID, foldID, runID and score are available
  #they are all factors
  score <- score/100
  filename <- paste("results500runs")
  Rdata_filename <- paste (filename,"Rdata",sep=".")
  
  
 
  
  rope_min  <- -0.01
  rope_max  <- 0.01
  
  how_many_classifiers <- max(as.numeric(classifierId))
  dsetsList <- unique(as.numeric(dsetId))
  how_many_dsets <- length(dsetsList)
  
  #debug
  #how_many_dsets <- 3
  #how_many_classifiers <- 3
  
  #generate always the same permutations
  set.seed(42)
  
  classifier_names <- levels(classifierId)
  how_many_comparisons <- how_many_classifiers*(how_many_classifiers-1)/2 #number of pairwise comparisons
  
  stanResults <- vector(mode='list', length=how_many_comparisons) 
  
  #fields to be filled during the multiple comparisons
  howManyPredictions <- how_many_comparisons*how_many_dsets
  classifierI <- vector()
  classifierJ <- vector()
  dset <- vector()
  hierDeltaEachDset <- vector()
  mleDeltaEachDset <- vector ()
  actual <- vector()
  
  counter <- 1
  
  for (i in 1: (how_many_classifiers-1) ) {
    for (j in (i+1) : how_many_classifiers){
      print('classifiers')
      show(c(i,j))
      
      #prepare the data for  the hierarchical test
      results <- cbind (as.numeric(classifierId), as.numeric(dsetId), score, foldID)
      resultsI  <- results[as.numeric(classifierId)==i,] 
      resultsJ  <- results[as.numeric(classifierId)==j,]
      
      #check results are consistently paired as for dataset and foldID
      #data are already properly sorted, so this control pases
      #otherwise we should sort the matrixes
      stopifnot( mean (resultsI[,c(2,4)]==resultsJ[,c(2,4)]) == 1)
      #diffIJ is a matrix, the first column stores the dSetID, the second the 100 differences for that dset
      diffIJ <- cbind (resultsI[,2] , resultsI[,3]-resultsJ[,3])
      
      
      #build matrix of total results
      x<-matrix(ncol  = length(classifierId), nrow = how_many_dsets )
      for (dsetIdx in 1:how_many_dsets) {
        tmp <- diffIJ [diffIJ[,1] == dsetIdx,] 
        x[dsetIdx,]  <- t (tmp [,2])
      }
      
      #xTrain  contains the first 100 results on each dset
      #we could eventually build a for cycle on this 
      xTrain<-x[,1:trainResults]
      
      #run the hierarchical test
      #we do not provide a simulation ID as this is run locally
      simulationID <- paste(as.character(i*10 + j),as.character(std_upper_bound),sep = "-")
      
      #debug
      #xTrain<-xTrain[1:3,1:5]
      #debug 2 chains
      currentResults <- hierarchical.test (xTrain,rho,rope_min,rope_max,simulationID,std_upper_bound)
      #store the results
      stanResults[[counter]] <- currentResults
      hierDeltaEachDset <- c(hierDeltaEachDset, currentResults$delta_each_dset)
      currentMleDeltaEachDset <- apply(xTrain,1,mean)
      mleDeltaEachDset <- c(mleDeltaEachDset, currentMleDeltaEachDset)
      actual <- c (actual, getAvgScore(i,score) - getAvgScore(j,score))
      classifierI <- c(classifierI,rep(i,how_many_dsets))
      classifierJ <- c(classifierJ,rep(j,how_many_dsets))
      dset <- c(dset,1:how_many_dsets)
    }
  }
  
  errHier <- hierDeltaEachDset - actual
  errMle <- mleDeltaEachDset -actual
  predictiveHalfCval <- data.frame(classifierI,classifierJ,dset,hierDeltaEachDset,mleDeltaEachDset,actual)
  predictiveHalfCval$errMle <- predictiveHalfCval$mleDeltaEachDset - predictiveHalfCval$actual
  predictiveHalfCval$errHier <- predictiveHalfCval$hierDeltaEachDset - predictiveHalfCval$actual
  write.matrix(predictiveHalfCval,file="predictiveHalfCval.csv",sep=",")
  
  save(stanResults, file = Rdata_filename)
  
  return (results)
  
}

