#compares the estimates yielded by the hierarchical
#model and the sign test on the next dsets

#todo: 
#1) fix signed-rank 
#2) consider different setd orders

predictiveComparison <- function() {
  
  library(MASS)
  library(matrixStats)
  library(rstan)
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  source ("hierarchical_test.R")
  source ("BayesianSignedRank.R")
  
  
  #utiliy function which gets the average score of the given classifier on a list of data sets
  getAvgScore <- function (currentClassifier,score, dsetList) {
    avg_score <- vector(length = length(dsetList))
    for ( ii in 1:length(dsetList) ){
      avg_score[ii] <- mean ( score [classifierId==currentClassifier & dsetId==dsetsList[ii]] )
    }
    return (avg_score)
  }
  
  
  samplingType="student"  
  std_upper_bound <- 1000
  standardization <- 1 #we run on standardized data
  
  
  load("uci_data.RData")
  
  #output parameters 
  
  #check arguments
  if (std_upper_bound<=1){
    stop("std_std_upper_bound should be larger than 1")
  }
  
  if (!  ( (samplingType=="normal") | (samplingType=="student") )){
    stop("wrong samplingType")
  }
  
  #for each couple of classifiers: p_value_sign_rank,prob_classA_delta0,prob_rope_delta0,prob_classB_delta0
  
  dsetId<- uci_classification$DatasetID
  
  
  
  classifierId <- uci_classification$ClassifierID
  #we run on accuracy
  score <- uci_classification$Percent.correct
  
  foldID <- uci_classification$Key.Fold
  nFolds <- max (foldID)
  rho=1/nFolds
  
  filename <- paste("uciResults_StdUpper",as.character(std_upper_bound),"samplingType",samplingType, sep="-")
  Rdata_filename <- paste (filename,"Rdata",sep=".")
  
  
  runID <- uci_classification$Key.Run
  #the fold ID discriminates between different folds and different runs
  foldID <- runID*10+foldID-10
  rope_min  <- -0.005
  rope_max  <- 0.005
  how_many_classifiers <- max(unique(classifierId))
  classifier_names <- c('naive Bayes','aode','hnb','j48','j48_grafted')  #hard-coded
  how_many_comparisons <- how_many_classifiers*(how_many_classifiers-1)/2 #number of pairwise comparisons
  
  #DEBUG
  # howManyTrainingDsets<-40
  howManyTrainingDsets<-3
  
  #repetitions<-5 
  repetitions<-1 
  
  
  #fields to be filled during the multiple comparisons
  p_value_sign_rank <- matrix(ncol = how_many_comparisons, nrow=repetitions)
  p_value_t_test <- matrix(ncol = how_many_comparisons, nrow=repetitions)
  median_difference <- matrix(ncol = how_many_comparisons, nrow=repetitions)
  classifierI <- vector(length = how_many_comparisons)
  classifierJ <- vector(length = how_many_comparisons)
  logLossHier <- matrix(ncol = how_many_comparisons, nrow=repetitions)
  logLossSignedRank <- matrix(ncol = how_many_comparisons, nrow=repetitions)
  maeHier <- matrix(ncol = how_many_comparisons, nrow=repetitions)
  maeSignedRank <- matrix(ncol = how_many_comparisons, nrow=repetitions)
  probLeftHier <- matrix(ncol = how_many_comparisons, nrow=repetitions)
  probRopeHier <- matrix(ncol = how_many_comparisons, nrow=repetitions)
  probRightHier <- matrix(ncol = how_many_comparisons, nrow=repetitions)
  probLeftSignedRank <- matrix(ncol = how_many_comparisons, nrow=repetitions)
  probRopeSignedRank <- matrix(ncol = how_many_comparisons, nrow=repetitions)
  probRightSignedRank <- matrix(ncol = how_many_comparisons, nrow=repetitions)
  leftCount <- matrix(ncol = how_many_comparisons, nrow=repetitions)
  ropeCount <- matrix(ncol = how_many_comparisons, nrow=repetitions)
  rightCount <- matrix(ncol = how_many_comparisons, nrow=repetitions)
  classifier_names <- c('naive Bayes','aode','hnb','j48','j48_grafted')
  
  dsetsList <- unique(dsetId);
  hierarchicalResults <- list()
  counter <- 1
  for (i in 1: (how_many_classifiers-1) ) {
    for (j in (i+1) : how_many_classifiers){
      
      show(c(i,j))
      
      classifierI[counter] <- i
      classifierJ[counter] <- j
      
      for (reps in 1:repetitions){
        #run the signed rank
        dsetList <- seq (from=1, to = max(dsetId), by =1)
        trainingDsets <- sample (dsetList,howManyTrainingDsets)
        testDsets <- setdiff (dsetList,trainingDsets)
        avgScoreI <- getAvgScore (i,score, trainingDsets)         
        avgScoreJ <- getAvgScore (j,score, trainingDsets)
        median_difference[reps,counter] <- median(avgScoreI-avgScoreJ)
        wilcoxonStat <-  wilcox.test (avgScoreI-avgScoreJ,alternative = "two.sided")
        p_value_sign_rank[reps,counter] <- wilcoxonStat$p.value
        
        #prepare the data for  the hierarchical test
        results <- cbind (classifierId, dsetId, score, foldID)
        resultsI  <- results[classifierId==i,] 
        resultsJ  <- results[classifierId==j,]
        
        
        #check results are consistently paired as for dataset and foldID
        #data are already properly sorted, so this control pases
        #otherwise we should sort the matrixes
        stopifnot( mean (resultsI[,c(2,4)]==resultsJ[,c(2,4)]) == 1)
        diffIJ <- cbind (resultsI[,2] , resultsI[,3]-resultsJ[,3])
        
        
        #build matrix of results to be parsed by hierarchical test
        x <- matrix(ncol  = max(foldID), nrow = howManyTrainingDsets )
        for (dsetIdx in 1:howManyTrainingDsets) {
          tmp <- diffIJ [diffIJ[,1] == dsetIdx,] 
          x[dsetIdx,]  <- t (tmp [,2])
        }
        
        
        #run the hierarchical test
        #we do not provide a simulation ID as this is run locally
        simulationID <- paste(as.character(i*10 + j),as.character(std_upper_bound),samplingType,sep = "-")
        hierarchicalResults[[counter]] <- hierarchical.test (x,rho,rope_min,rope_max,simulationID,std_upper_bound,samplingType)
        
        #hierarchical predictive
        probLeftHier[reps,counter]  <- hierarchicalResults[[counter]]$nextDelta$left
        probRopeHier[reps,counter]  <- hierarchicalResults[[counter]]$nextDelta$rope
        probRightHier[reps,counter] <- hierarchicalResults[[counter]]$nextDelta$right
        
        #sign test
        signedRankResults<-BayesianSignedRank(avgScoreI-avgScoreJ,rope_min,rope_max)
        probRopeSignedRank[reps,counter]  <- signedRankResults$probRope
        probLeftSignedRank[reps,counter]  <- signedRankResults$probLeft
        probRightSignedRank[reps,counter] <- signedRankResults$probRight
        
        #analyze the future dsets
        avgScoreI <- getAvgScore (i,score, (length(trainingDsets)+1):max(dsetId))         
        avgScoreJ <- getAvgScore (j,score, (length(trainingDsets)+1):max(dsetId))
        nextDifferences <- avgScoreI-avgScoreJ
        leftCount[reps,counter]<-  sum(nextDifferences<rope_min)
        rightCount[reps,counter]<- sum(nextDifferences>rope_max)
        ropeCount[reps,counter]<- length(nextDifferences) - leftCount[reps,counter] - rightCount[reps,counter]
        logLossHier[reps,counter]<- -1*(leftCount[reps,counter]*log(probLeftHier[reps,counter]) + rightCount[reps,counter]*log(probRightHier[reps,counter]) + ropeCount[reps,counter]*log(probRopeHier[reps,counter]))
        logLossSignedRank[reps,counter]<- -1*(leftCount[reps,counter]*log(probLeftSignedRank[reps,counter]) + rightCount[reps,counter]*log(probRightSignedRank[reps,counter]) + ropeCount[reps,counter]*log(probRopeSignedRank[reps,counter]))
        maeHier[reps,counter]<- leftCount[reps,counter]*(1-probLeftHier[reps,counter]) + rightCount[reps,counter]*(1-probRightHier[reps,counter]) + ropeCount[reps,counter]*(1-probRopeHier[reps,counter])
        maeSignedRank[reps,counter]<- leftCount[reps,counter]*(1-probLeftSignedRank[reps,counter]) + rightCount[reps,counter]*(1-probRightSignedRank[reps,counter]) + ropeCount[reps,counter]*(1-probRopeSignedRank[reps,counter])
        maeHier[reps,counter]<-maeHier[reps,counter]/length(nextDifferences)
        maeSignedRank[reps,counter]<-maeSignedRank[reps,counter]/length(nextDifferences)
      }
      
      
      #       #write results after having analyzed a pair of classifier
      #       results_matrix<-data.frame(
      #         classifierI=classifierI,
      #         classifierJ=classifierJ,
      #         probLeftHier=mean(probLeftHier[,counter]),
      #         probRopeHier=mean(probRopeHier[,counter]),
      #         probRightHier=mean(probRightHier[,counter]),
      #         probLeftSign=mean(probLeftSign[,counter]),
      #         probRopeSign=mean(probRopeSign[,counter]),
      #         probRightSign=mean(probRightSign[,counter]),
      #         ropeCount=mean(ropeCount[,counter]),
      #         leftCount=mean(leftCount[,counter]),
      #         rightCount=mean(rightCount[,counter]),
      #         logLossHier=mean(logLossHier[,counter]),
      #         logLossSign=mean(logLossSign[,counter]),
      #         maeHier=mean(maeHier[,counter]),
      #         maeSignedRank=mean(maeSignedRank[,counter])
      #       )
      #       filename<-"predictiveSignHier"
      #       csv_filename <- paste (filename,"csv",sep=".")
      #       write.matrix(results_matrix, file = csv_filename, sep = ",")
      
      #update for the next iteration
      counter <- counter + 1
      
    }
  }
  
  results_matrix<-data.frame(
    classifierI=classifierI,
    classifierJ=classifierJ,
    probLeftHier=colMeans(probLeftHier),
    probRopeHier=colMeans(probRopeHier),
    probRightHier=colMeans(probRightHier),
    probLeftSignedRank=colMeans(probLeftSignedRank),
    probRopeSignedRank=colMeans(probRopeSignedRank),
    probRightSignedRank=colMeans(probRightSignedRank),
    ropeCount=colMeans(ropeCount),
    leftCount=colMeans(leftCount),
    rightCount=colMeans(rightCount),
    logLossHier=colMeans(logLossHier),
    logLossSignedRank=colMeans(logLossSignedRank),
    maeHier=colMeans(maeHier),
    maeSignedRank=colMeans(maeSignedRank)
  )
  filename<-"predictiveSignHier"
  csv_filename <- paste (filename,"csv",sep=".")
  write.matrix(results_matrix, file = csv_filename, sep = ",")
  
  
}

