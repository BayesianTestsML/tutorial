#THIS CODE IS UNDER DEVELOPMENT, NOT YET STABLE
#compares the predictive estimates yielded by the MLE and  the hierarchical model. 
#It infers the mean difference on the first 50 cross-validation  observations 
#and predicts the mean of the remaining 50 results.
#samplingType might be gaussian or student
predictive_halfCvalFriedman <- function(samplingType="student") {
  library(MASS)
  library(rstan)
  rstan_options(auto_write = TRUE)
  
  options(mc.cores = parallel::detectCores())
  source ("hierarchical_test.R")
  std_upper_bound <- 1000
  standardization <- 1 
  
  #utiliy function which gets the average score of the given classifier on each data set
  getAvgScore <- function (currentClassifier,score) {
    avg_score <- vector(length = how_many_dsets, mode = "double");
    for ( ii in 1:how_many_dsets ){
      avg_score[ii] <- mean ( score [classifierId==currentClassifier & dsetId==dsetsList[ii]] )
    }
    return (avg_score)
  }
  
  load("friedmanData.RData")
  
  
  #for each couple of classifiers: p_value_sign_rank,prob_classA_delta0,prob_rope_delta0,prob_classB_delta0
  dsetName <- factor(resultsWeka$Key_Dataset)
  dsetId<- as.numeric(dsetName)
  classifierName <- factor(resultsWeka$Key_Scheme)
  classifierId <- as.numeric(classifierName)
  #we run on accuracy
  score <- resultsWeka$Percent_correct
  
  foldID <- resultsWeka$Key_Fold
  nFolds <- max (foldID)
  rho=1/nFolds
  
  filename <- paste("friedmanResultsHalfCval",samplingType,sep="_")
  Rdata_filename <- paste (filename,"Rdata",sep=".")
  
  
  runID <- resultsWeka$Key_Run
  #the fold ID discriminates between different folds and different runs
  foldID <- runID*10+foldID-10
  
  rope_min  <- -0.01
  rope_max  <- 0.01
  
  how_many_classifiers <- max(unique(classifierId))
  dsetsList <- unique(dsetId)
  how_many_dsets <- length(dsetsList)
  
  #debug
  #how_many_dsets <- 3
  #how_many_classifiers <- 3
  
  
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
      results <- cbind (classifierId, dsetId, score, foldID)
      resultsI  <- results[classifierId==i,] 
      resultsJ  <- results[classifierId==j,]
      
      #check results are consistently paired as for dataset and foldID
      #data are already properly sorted, so this control pases
      #otherwise we should sort the matrixes
      stopifnot( mean (resultsI[,c(2,4)]==resultsJ[,c(2,4)]) == 1)
      #diffIJ is a matrix, the first column stores the dSetID, the second the 100 differences for that dset
      diffIJ <- cbind (resultsI[,2] , resultsI[,3]-resultsJ[,3])
      
      
      #build matrix of total results
      x<-matrix(ncol  = max(foldID), nrow = how_many_dsets )
      for (dsetIdx in 1:how_many_dsets) {
        tmp <- diffIJ [diffIJ[,1] == dsetIdx,] 
        x[dsetIdx,]  <- t (tmp [,2])
      }
      
      xTrain<- x[,1:10]
      #xTest  contains all the results (100 for each dset) for the dset in the second half of the permutation
      xTest<- x[,11:100]
      
      #run the hierarchical test
      #we do not provide a simulation ID as this is run locally
      simulationID <- paste(as.character(i*10 + j),as.character(std_upper_bound),sep = "-")
      
      #debug
      #xTrain<-xTrain[1:3,1:5]
      #debug 2 chains
      currentResults <- hierarchical.test (xTrain,rho,rope_min,rope_max,simulationID,std_upper_bound,samplingType)
      #store the results
      stanResults[[counter]] <- currentResults
      hierDeltaEachDset <- c(hierDeltaEachDset, currentResults$delta_each_dset)
      currentMleDeltaEachDset <- apply(xTrain,1,mean)
      mleDeltaEachDset <- c(mleDeltaEachDset, currentMleDeltaEachDset)
      actual <- c (actual, apply(xTest,1,mean))
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
  csvFilename = paste("predictiveHalfCvalFriedman",samplingType,sep='-')
  csvFilename = paste(csvFilename,"csv",sep='.')
  write.matrix(predictiveHalfCval,csvFilename,sep=",")
  
  save(stanResults, file = Rdata_filename)
  
  return (results)
  
}

