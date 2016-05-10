#compares the predictive estimates yielded by the Gaussian and by the Student. It infers the model on the first 50 cross-validation  observations 
#and computed the marginal likelihood of the remaining 50
predictive_halfCval <- function() {
  
  library(MASS)
  library(matrixStats)
  library(rstan)
  rstan_options(auto_write = TRUE)
  
  options(mc.cores = parallel::detectCores())
  source ("hierarchical_test.R")
  source ("getMargLik.R")
  std_upper_bound <- 1000
  standardization <- 1 #we run on standardized data
  
  #utiliy function which gets the average score of the given classifier on each data set
  getAvgScore <- function (currentClassifier,score) {
    avg_score <- vector(length = how_many_dsets, mode = "double");
    for ( ii in 1:how_many_dsets ){
      avg_score[ii] <- mean ( score [classifierId==currentClassifier & dsetId==dsetsList[ii]] )
    }
    return (avg_score)
  }
  
  load("uci_data.RData")
  
  #output parameters 
  
  #check arguments
  if (std_upper_bound<=1){
    stop("std_std_upper_bound should be larger than 1")
  }
  
  
  #for each couple of classifiers: p_value_sign_rank,prob_classA_delta0,prob_rope_delta0,prob_classB_delta0
  
  dsetId<- uci_classification$DatasetID
  classifierId <- uci_classification$ClassifierID
  #we run on accuracy
  score <- uci_classification$Percent.correct
  
  foldID <- uci_classification$Key.Fold
  nFolds <- max (foldID)
  rho=1/nFolds
  
  filename <- paste("uciResults_predictive")
  Rdata_filename <- paste (filename,"Rdata",sep=".")
  
  
  runID <- uci_classification$Key.Run
  #the fold ID discriminates between different folds and different runs
  foldID <- runID*10+foldID-10
  
  rope_min  <- -0.01
  rope_max  <- 0.01
  
  how_many_classifiers <- max(unique(classifierId))
  dsetsList <- unique(dsetId)
  how_many_dsets <- length(dsetsList)
  classifier_names <- c('naive Bayes','aode','hnb','j48','j48_grafted')  #hard-coded
  how_many_comparisons <- how_many_classifiers*(how_many_classifiers-1)/2 #number of pairwise comparisons
  
  #fields to be filled during the multiple comparisons
  margLikNormal<-matrix(nrow = how_many_comparisons,ncol = how_many_dsets)
  margLikStudent<-matrix(nrow = how_many_comparisons,ncol = how_many_dsets)
  classifierI <- vector(length = how_many_comparisons, mode = "integer")
  classifierJ <- vector(length = how_many_comparisons, mode = "integer")
  studentTime <- vector(length = how_many_comparisons)
  normalTime <- vector(length = how_many_comparisons)
  hierarchicalResultsStudent <- list()
  hierarchicalResultsNormal <- list()
  
  counter <- 1
  for (i in 1: (how_many_classifiers-1) ) {
    for (j in (i+1) : how_many_classifiers){
      #debug
      i<-4
      j<-5
      
      show(c(i,j))
      
      classifierI[counter] <- i
      classifierJ[counter] <- j
      
      
      #prepare the data for  the hierarchical test
      results <- cbind (classifierId, dsetId, score, foldID)
      resultsI  <- results[classifierId==i,] 
      resultsJ  <- results[classifierId==j,]
      
      #sort by dataset and then by fold
      #resultsI<-mat.sort(resultsI, c(2,4))
      #resultsJ<-mat.sort(resultsJ, c(2,4))
      
      #check results are consistently paired as for dataset and foldID
      #data are already properly sorted, so this control pases
      #otherwise we should sort the matrixes
      stopifnot( mean (resultsI[,c(2,4)]==resultsJ[,c(2,4)]) == 1)
      diffIJ <- cbind (resultsI[,2] , resultsI[,3]-resultsJ[,3])
      
      
      #build matrix of total results
      x<-matrix(ncol  = max(foldID), nrow = how_many_dsets )
      for (dsetIdx in 1:how_many_dsets) {
        tmp <- diffIJ [diffIJ[,1] == dsetIdx,] 
        x[dsetIdx,]  <- t (tmp [,2])
      }
      
      xTrain<-x[,1:ncol(x)/2]
      xTest<- x[,(ncol(x)/2+1):ncol(x)]
      
      #run the hierarchical test
      #we do not provide a simulation ID as this is run locally
      samplingType<-"student"
      simulationID <- paste(as.character(i*10 + j),as.character(std_upper_bound),samplingType,sep = "-")
      startTime<-Sys.time()
      hierarchicalResultsStudent[[counter]] <- hierarchical.test (xTrain,rho,rope_min,rope_max,simulationID,std_upper_bound,samplingType)
      stopTime<-Sys.time()
      studentTime[counter]<-(stopTime-startTime)
      #getMargLik returns a row vector 1 x nDsets
      margLikStudent[counter,]<-getMargLik(hierarchicalResultsStudent[[counter]],xTest,rho);
      write.matrix(margLikStudent,file="margLikStudent.csv",sep=",")
      
      samplingType<-"normal"
      simulationID <- paste(as.character(i*10 + j),as.character(std_upper_bound),samplingType,sep = "-")
      startTime<-Sys.time()
      hierarchicalResultsNormal[[counter]] <- hierarchical.test (xTrain,rho,rope_min,rope_max,simulationID,std_upper_bound,samplingType)
      stopTime<-Sys.time()
      normalTime[counter]<-(stopTime-startTime)
      margLikNormal[counter,]<-getMargLik(hierarchicalResultsNormal[[counter]],xTest,rho);
      write.matrix(margLikNormal,file="margLikNormal.csv",sep=",")
      
      counter <- counter + 1
    }
  }
  
  results <- list() 
  results[[1]] <- list('margLikNormal'=margLikNormal,
                       'margLikStudent'=margLikStudent)
  
  save(results, file = Rdata_filename)
  
  return (results)
  
}

