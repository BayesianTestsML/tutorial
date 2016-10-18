sensitivityNormalStudent <- function (class1, class2){
  #class1 and class2 are the two classifier being compared
  #'naive Bayes','aode','hnb','j48','j48_grafted' are coded as 1,2,3,4,5 respectively
  #'e.g. to compare naive Bayes and aode: sensitivityNormalStudent(1,2)
  #in the paper we present the pairs 1-3 and 4-5
  #1) check the inferences of different priors over delta_i: normal, student (Krushcke), student(Juanez), student (GC).
  #2) verify the sensitivity of the shrinkage estimates
  #3) check the log-posterior
  

  library(rstan)
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  source ("hierarchical_test.R")

  #preliminaries
  stdUpperBound <- 1000
  standardization <- 1 
  #to have maximum reliability when comparing real classifiers
  chains = 8
  
  load("uci_data.RData")
  nFolds <- max (uci_classification$Key.Fold)
  rho=1/nFolds
  #this foldID goes between 1 and 100
  foldID <- uci_classification$Key.Run*10+uci_classification$Key.Fold-10
  rope_min  <- -0.01
  rope_max  <- 0.01
  
  #prepare the data for  the hierarchical test
  results <- data.frame (classifierID = uci_classification$ClassifierID, 
                         dsetID=uci_classification$DatasetID, 
                         accuracy=uci_classification$Percent.correct, 
                         fold=foldID)
  
  
  
  diffResults  <- results[results$classifierID==class1,]
  results2  <- results[results$classifierID==class2,]
  stopifnot( mean (diffResults$dsetID==results2$dsetID) ==1)
  diffResults$diff <- diffResults$accuracy - results2$accuracy
  
  
  #build matrix of results to be parsed by hierarchical test
  howManyDsets <- max(diffResults$dsetID)
  x<-matrix(ncol  = max(foldID), nrow =  howManyDsets)
  
  for (dsetIdx in 1:howManyDsets) {
    tmp <- diffResults$diff[diffResults$dsetID == dsetIdx] 
    x[dsetIdx,]  <- t (tmp)
  }
  
  #debug 
  chains <- 4
#debug:currently comment out the  inference of the full model.
  
#   simulationID <- paste('class',class1,'class',class2,"Kruschke",sep ='')
#   hierPosteriorKru <- hierarchical.test (x,rho,rope_min,rope_max,simulationID,stdUpperBound,"studentKruschke",chains)
#     
#   simulationID <- paste('class',class1,'class',class2,"Juanez",sep ='')
#   hierPosteriorJua <- hierarchical.test (x,rho,rope_min,rope_max,simulationID,stdUpperBound,"studentJuanez",chains)
#    
#   simulationID <- paste('class',class1,'class',class2,"GC",sep ='')
#   hierPosterior <- hierarchical.test (x,rho,rope_min,rope_max,simulationID,stdUpperBound,"student",chains)
#  
#   simulationID <- paste('class',class1,'class',class2,"Gaussian",sep ='')
#   hierPosteriorGauss <- hierarchical.test (x,rho,rope_min,rope_max,simulationID,stdUpperBound,"gaussian",chains)
  
  
   
  #now implement log-posterior check 
  #first, infer the model on half the data.
  trainX <- x[,1:(ncol(x)/2)]
  testX  <- x[,(ncol(x)/2 + 1):ncol(x)]
  simulationID <- paste('class',class1,'class',class2,"Kruschke",sep ='')
  halfPosteriorKru <- hierarchical.test (trainX,rho,rope_min,rope_max,simulationID,stdUpperBound,"studentKruschke",chains)
  
  #debug: comment out the other cases
#   simulationID <- paste('class',class1,'class',class2,"Juanez",sep ='')
#   halfPosteriorJua <- hierarchical.test (trainX,rho,rope_min,rope_max,simulationID,stdUpperBound,"studentJuanez",chains)
#   
#   simulationID <- paste('class',class1,'class',class2,"GC",sep ='')
#   halfPosterior <- hierarchical.test (trainX,rho,rope_min,rope_max,simulationID,stdUpperBound,"student",chains)
#   
#   simulationID <- paste('class',class1,'class',class2,"Gaussian",sep ='')
#   halfPosteriorGauss <- hierarchical.test (trainX,rho,rope_min,rope_max,simulationID,stdUpperBound,"gaussian",chains)
  
  hierPosteriorKru$logPredictive <- logPredictive(halfPosteriorKru,testX, rho)
  hierPosteriorJua$logPredictive <- logPredictive(halfPosteriorJua,testX, rho)
  hierPosterior$logPredictive <- logPredictive(halfPosterior,testX, rho)
  hierPosteriorGauss$logPredictive <- logPredictive(halfPosteriorGauss,testX, rho)
  
  #and save to Rdata directory
}