#compares the predictive estimates yielded by the MLE and  the hierarchical model. 
#It infers the mean difference on the first 50 cross-validation  observations 
#and predicts the mean of the remaining 50 results.

#samplingType might be gaussian or student
#groups might be "all" (all dsets), 5, 10, 25, 50 (which selects only dsets with 5/10/25/50 features)
#for samples: 50/250/500/1000 
predictive_halfCvalFriedman <- function(samplingType="student", group='all') {
  library(MASS)
  library(rstan)
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  source ("hierarchical_test.R")
  std_upper_bound <- 1000
  standardization <- 1 
  load("friedmanData.RData")
  set.seed(0)
  
  
  #for each couple of classifiers: p_value_sign_rank,prob_classA_delta0,prob_rope_delta0,prob_classB_delta0
  dsetName <- factor(resultsWeka$Key_Dataset)
  
  #SELECTION OF THE DATA SETS TO BE ANALYZED: ALL OR A SUB-GROUP
  token <- "fri" #if no group is set, we choose all dset
  
  #we build a token like _5-
  if (group != 'all') {
    token <- paste('_',group,sep='')
    #to group on the number of features
    token <- paste(token,'-',sep='')
    #to group on the number of instances
    #    token <- paste(token,'_',sep='')
  }
  
  #arrays of booleans, if the dsetName contains the token
  #we will select the rows from weka results using dsetIsValid
  dsetIsValid <- grepl(token,dsetName)  
  dsetId<- as.numeric(dsetName[dsetIsValid])
  classifierName <- factor(resultsWeka$Key_Scheme)
  classifierName <- classifierName[dsetIsValid]
  
  classifierId <- as.numeric(classifierName)
  #we run on accuracy
  score <- resultsWeka$Percent_correct
  score <- score[dsetIsValid]
  # to manage weka saved accuracies whihc range between 0 and 100
  if (max(score)>1) {
    score <- score/100
  }
  
  
  foldID <- resultsWeka$Key_Fold
  foldID <- foldID[dsetIsValid]
  nFolds <- max (foldID)
  rho=1/nFolds
  
  filename <- paste("friedmanResultsHalfCval",samplingType,sep="_")
  filename <- paste(filename,group,sep="_")
  Rdata_filename <- paste (filename,"Rdata",sep=".")
  
  
  runID <- resultsWeka$Key_Run
  runID <- runID[dsetIsValid]
  #the fold ID discriminates between different folds and different runs
  foldID <- runID*10+foldID-10
  
  rope_min  <- -0.01
  rope_max  <- 0.01
  
  how_many_classifiers <- max(unique(classifierId))
  dsetsList <- unique(dsetId)
  how_many_dsets <- length(dsetsList)
  
  
  
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
      validDsets <- unique (dsetId)
      #check results are consistently paired as for dataset and foldID
      #data are already properly sorted, so this control pases
      #otherwise we should sort the matrixes
      stopifnot( mean (resultsI[,c(2,4)]==resultsJ[,c(2,4)]) == 1)
      #diffIJ is a matrix, the first column stores the dSetID, the second the 100 differences for that dset
      diffIJ <- cbind (resultsI[,2] , resultsI[,3]-resultsJ[,3])
      
      
      #build matrix of total results
      x<-matrix(ncol  = max(foldID), nrow = how_many_dsets )
      for (idx in 1:how_many_dsets) {
        currentDset <- validDsets[idx] 
        tmp <- diffIJ [diffIJ[,1] == currentDset,] 
        x[idx,]  <- t (tmp [,2])
      }
      
      xTrain<- x[,1:5]
      #xTest  contains all the results (100 for each dset) for the dset in the second half of the permutation
      xTest<- x[,6:100]
      
      #debug code on synthetic data
      #xTrain <- matrix(data= rnorm(n=400, mean = 0, sd = 0.1), nrow = 80, ncol = 5)
      #debug code on synthetic data
      #xTest <- matrix(data=rnorm(n=800, mean = 0, sd = 0.1), nrow = 80, ncol = 10)
      
      #run the hierarchical test
      #we do not provide a simulation ID as this is run locally
      simulationID <- paste(as.character(i*10 + j),as.character(std_upper_bound),sep = "-")
      
      #debug
      #we are running only 4 chains for speed.
      currentResults <- hierarchical.test (xTrain,rho,rope_min,rope_max,simulationID,std_upper_bound,samplingType,chains=4)
      #store the results
      stanResults[[counter]] <- currentResults
      hierDeltaEachDset <- c(hierDeltaEachDset, currentResults$delta_each_dset)
      currentMleDeltaEachDset <- apply(xTrain,1,mean)
      mleDeltaEachDset <- c(mleDeltaEachDset, currentMleDeltaEachDset)
      actual <- c (actual, apply(xTest,1,mean))
      classifierI <- c(classifierI,rep(i,how_many_dsets))
      classifierJ <- c(classifierJ,rep(j,how_many_dsets))
      dset <- c(dset,1:how_many_dsets)
      
      #debug code to be removed
      maeMle <- mean(abs (actual - mleDeltaEachDset)) 
      maeHier <- mean(abs (actual - hierDeltaEachDset)) 
      print (paste ('maeMle', maeMle))
      print (paste ('maeHier', maeHier))
    }
  }
  
  errHier <- hierDeltaEachDset - actual
  errMle <- mleDeltaEachDset -actual
  predictiveHalfCval <- data.frame(classifierI,classifierJ,dset,hierDeltaEachDset,mleDeltaEachDset,actual)
  predictiveHalfCval$errMle <- predictiveHalfCval$mleDeltaEachDset - predictiveHalfCval$actual
  predictiveHalfCval$errHier <- predictiveHalfCval$hierDeltaEachDset - predictiveHalfCval$actual
  csvFilename = paste(filename,"csv",sep='.')
  write.matrix(predictiveHalfCval,csvFilename,sep=",")
  save(stanResults, file = Rdata_filename)
  
  return (results)
  
}

