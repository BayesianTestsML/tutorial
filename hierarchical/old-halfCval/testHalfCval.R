#compares the predictive estimates yielded by the MLE and  the hierarchical model. 
#on a single dset
testHalfCval <- function(samplingType="student", group='all') {
  library(MASS)
  library(rstan)
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  source ("hierarchical_test.R")
  std_upper_bound <- 1000
  standardization <- 1 
  load("friedmanData.RData")
  set.seed(88)
  
  
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
  #to manage weka saved accuracies whihc range between 0 and 100
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
  
  #debug 
  how_many_classifiers <- 2
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
      validDsets <- 42 #for debugging reasons
      #check results are consistently paired as for dataset and foldID
      #data are already properly sorted, so this control pases
      #otherwise we should sort the matrixes
      stopifnot( mean (resultsI[,c(2,4)]==resultsJ[,c(2,4)]) == 1)
      #diffIJ is a matrix, the first column stores the dSetID, the second the 100 differences for that dset
      diffIJ <- cbind (resultsI[,2] , resultsI[,3]-resultsJ[,3])
      
      
      #build matrix of total results
      x<-matrix(ncol  = max(foldID), nrow = length(validDsets) )
      for (idx in 1:length(validDsets)) {
        currentDset <- validDsets[idx] 
        tmp <- diffIJ [diffIJ[,1] == currentDset,] 
        x[idx,]  <- t (tmp [,2])
      }
      
      xTrain<- t(cbind(x[1:10],x[11:20],x[21:30],x[31:40],x[41:50],x[51:60],x[61:70],x[71:80],
                     x[81:90],x[91:100],x[15:24],x[25:34]))
      #xTest  contains all the results (100 for each dset) for the dset in the second half of the permutation
      xTest<- t(cbind(x[11:100],x[21:100],x[c(1:20,31:100)],x[c(1:30,41:100)],
                     x[c(1:40,51:100)],x[c(1:50,61:100)],x[c(1:60,71:100)],
                     x[c(1:70,81:100)],
                     x[c(1:80,91:100)],
                     x[1:90],
                     x[c(1:14,25:100)],
                     x[c(1:24,35:100)]))
      
      #run the hierarchical test
      #we do not provide a simulation ID as this is run locally
      simulationID <- paste(as.character(i*10 + j),as.character(std_upper_bound),sep = "-")
      
      #low number of chains
      currentResults <- hierarchical.test (xTrain,rho,rope_min,rope_max,simulationID,std_upper_bound,samplingType,chains=4)
      #store the results
      stanResults[[counter]] <- currentResults
      hierDeltaEachDset <- c(hierDeltaEachDset, currentResults$delta_each_dset)
      currentMleDeltaEachDset <- apply(xTrain,1,mean)
      mleDeltaEachDset <- c(mleDeltaEachDset, currentMleDeltaEachDset)
      actual <- c (actual, apply(xTest,1,mean))
      classifierI <- c(classifierI,rep(i,dim(xTrain)[1]))
      classifierJ <- c(classifierJ,rep(j,dim(xTrain)[1]))
      dset <- c(dset,1:how_many_dsets)
    }
  }
  
  maeMle <- mean(abs (actual - mleDeltaEachDset)) 
  maeHier <- mean(abs (actual - hierDeltaEachDset)) 
  print (paste ('maeMle', maeMle))
  print (paste ('maeHier', maeHier))
  # return (results)
  
}

