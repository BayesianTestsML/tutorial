analyze_uci_with_rope <- function() {
  
  library(MASS)
  library(matrixStats)
  library(rstan)
  rstan_options(auto_write = TRUE)
  
  options(mc.cores = parallel::detectCores())
  source ("hierarchical_test.R")
  
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
  #to have maximum reliability when comparing real classifiers
  chains = 8 
  
  #output parameters 
  
  #check arguments
  if (std_upper_bound<=1){
    stop("std_std_upper_bound should be larger than 1")
  }
  
  samplingType="student" #default choice
  if (!  ( (samplingType=="normal") | (samplingType=="student") )){
    stop("wrong samplingType")
  }
  
  #for each couple of classifiers: p_value_sign_rank,prob_classANextDelta,prob_ropeNextDelta,prob_classBNextDelta
  
  dsetId<- uci_classification$DatasetID
  classifierId <- uci_classification$ClassifierID
  #we run on accuracy
  score <- uci_classification$Percent.correct
  
  foldID <- uci_classification$Key.Fold
  nFolds <- max (foldID)
  rho=1/nFolds
    
  filename <- paste("uciResults_StdUpper",as.character(std_upper_bound),"samplingType",samplingType, sep="-")
  if (rho==0)  {
    filename <- paste(filename,"noCorrel",sep="-")
  }
  
  Rdata_filename <- paste (filename,"Rdata",sep=".")
  
  
  runID <- uci_classification$Key.Run
  #the fold ID discriminates between different folds and different runs
  foldID <- runID*10+foldID-10
  
  rope_min  <- -0.01
  rope_max  <- 0.01
  
  how_many_classifiers <- max(unique(classifierId))
  
  classifier_names <- c('naive Bayes','aode','hnb','j48','j48_grafted')  #hard-coded
  how_many_comparisons <- how_many_classifiers*(how_many_classifiers-1)/2 #number of pairwise comparisons
  
  #fields to be filled during the multiple comparisons
  p_value_sign_rank <- vector(length = how_many_comparisons, mode = "double")
  p_value_t_test <- vector(length = how_many_comparisons, mode = "double")
  median_difference <- vector(length = how_many_comparisons, mode = "double")
  prob_classANextDelta <- vector(length = how_many_comparisons, mode = "double")
  prob_ropeNextDelta <- vector(length = how_many_comparisons, mode = "double")
  prob_classBNextDelta <- vector(length = how_many_comparisons, mode = "double")
  prob_positiveNextDelta<- vector(length = how_many_comparisons, mode = "double")
  prob_negativeNextDelta<- vector(length = how_many_comparisons, mode = "double")
  classifierI <- vector(length = how_many_comparisons, mode = "integer")
  classifierJ <- vector(length = how_many_comparisons, mode = "integer")
  
  
  dsetsList <- unique(dsetId);
  
  how_many_dsets <- length(dsetsList)
  
  hierarchicalResults <- list()
  counter <- 1
  for (i in 1: (how_many_classifiers-1) ) {
    for (j in (i+1) : how_many_classifiers){
      chains <-4
      show(c(i,j))
      
      classifierI[counter] <- i
      classifierJ[counter] <- j
      
      #run the signed rank
      avgScoreI <- getAvgScore (i,score)         
      avgScoreJ <- getAvgScore (j,score)
      median_difference[counter] <- median(avgScoreI-avgScoreJ)
      wilcoxonStat <-  wilcox.test (avgScoreI-avgScoreJ,alternative = "two.sided");
      p_value_sign_rank[counter] <- wilcoxonStat$p.value
      tTestStat <-  t.test (avgScoreI-avgScoreJ,alternative = "two.sided");
      p_value_t_test[counter] <- tTestStat$p.value
      
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
      
      
      #build matrix of results to be parsed by hierarchical test
      x<-matrix(ncol  = max(foldID), nrow = how_many_dsets )
      for (dsetIdx in 1:how_many_dsets) {
        tmp <- diffIJ [diffIJ[,1] == dsetIdx,] 
        x[dsetIdx,]  <- t (tmp [,2])
      }
      
      
      #run the hierarchical test
      #we do not provide a simulation ID as this is run locally
      startTime<-Sys.time()
      simulationID <- paste(as.character(i*10 + j),as.character(std_upper_bound),samplingType,".dat",sep = "-")
      
      hierarchicalResults[[counter]] <- hierarchical.test (x=x, rope_min=rope_min, rope_max=rope_max, sample_file=simulationID,
                                                           chains = chains,
                                                           samplingType = 'student')
      
      stopTime<-Sys.time()
      show(startTime-stopTime)
      
      #classA is the first classifier, classB the second classifier
      #delta = acc(classA) - acc(classB)
      #thus in the posterior the right tail of delta is the tail in favor
      #of classA
      
      prob_classANextDelta[counter]  <- hierarchicalResults[[counter]]$nextDelta$right
      prob_ropeNextDelta[counter]  <- hierarchicalResults[[counter]]$nextDelta$rope
      prob_classBNextDelta[counter] <- hierarchicalResults[[counter]]$nextDelta$left
      
      prob_positiveNextDelta[counter]  <- hierarchicalResults[[counter]]$nextDelta$positive
      prob_negativeNextDelta[counter]  <- hierarchicalResults[[counter]]$nextDelta$negative
      #waic not useful, we do not track it
      #waic[counter] <- hierarchicalResults[[counter]]$waic
      
      counter <- counter + 1
    }
  }
  
  classifierIString <-vector(length = how_many_comparisons, mode = "character")
  classifierJString <-vector(length = how_many_comparisons, mode = "character")
  
  classifier_names <- c('naive Bayes','aode','hnb','j48','j48_grafted');
  
  results_matrix<-data.frame(
    classifierI=classifierI,
    classifierJ=classifierJ,
    median_difference=median_difference,
    p_value_sign_rank=p_value_sign_rank,
    p_value_t_test=p_value_t_test,
    prob_classANextDelta=prob_classANextDelta,
    prob_ropeNextDelta=prob_ropeNextDelta,
    prob_classBNextDelta=prob_classBNextDelta,
    prob_positiveNextDelta=prob_positiveNextDelta,
    prob_negativeNextDelta=prob_negativeNextDelta,
    rope_min=rope_min,
    rope_max=rope_max
  )
  
  csv_filename <- paste (filename,"csv",sep=".")
  
  write.matrix(results_matrix, file = csv_filename, sep = ",")
  
  results <- list() 
  results[[1]] <- list('classifierI'=classifierI,
                       'classifierJ'=classifierJ,
                       'p_value_sign_rank'=p_value_sign_rank,
                       'median_difference'=median_difference,
                       'prob_classANextDelta'=prob_classANextDelta,
                       'prob_ropeNextDelta'=prob_ropeNextDelta,
                       'prob_classBNextDelta'=prob_classBNextDelta
  )
  
  
  save(hierarchicalResults, file = Rdata_filename)
  
  return (results)
  
}

