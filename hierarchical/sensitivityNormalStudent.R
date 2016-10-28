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
  source ("logPredictive.R")
  
  #preliminaries
  stdUpperBound <- 1000
  standardization <- 1 
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
  
  #DEBUG: LET'S COMMENT OUT ALL THE OTHER VERSION OF THE HIER MODEL  
  #   simulationID <- paste('class',class1,'class',class2,"Kruschke",sep ='')
  #   hierPosteriorKru <- hierarchical.test (x,rho,rope_min,rope_max,simulationID,stdUpperBound,"studentKruschke",chains)
  #   
  #   simulationID <- paste('class',class1,'class',class2,"Juanez",sep ='')
  #   hierPosteriorJua <- hierarchical.test (x,rho,rope_min,rope_max,simulationID,stdUpperBound,"studentJuanez",chains)
  
  # simulationID <- paste('class',class1,'class',class2,"Gaussian",sep ='')
  # hierPosteriorGauss <- hierarchical.test (x,rho,rope_min,rope_max,simulationID,stdUpperBound,"gaussian",chains)
  
  simulationID <- paste('class',class1,'class',class2,"GC",sep ='')
  hierPosterior <- hierarchical.test (x = x,sample_file = simulationID,samplingType = "student")
  
  simulationID <- paste('class',class1,'class',class2,"GCsens",sep ='')
  alphaBeta = list('lowerAlpha' =0.5,'upperAlpha'= 3,'lowerBeta' = 0.005,'upperBeta' = 0.05)
  hierPosteriorSens <- hierarchical.test (x = x,sample_file = simulationID,samplingType = "student", alphaBeta = alphaBeta)
  
  
  #computation of marginal likelihood on half the data.  
  #   logPredictiveValues <- list()
  #   #now implement log-posterior check 
  #   #first, infer the model on half the data.
  #   trainX <- x[,1:(ncol(x)/2)]
  #   testX  <- x[,(ncol(x)/2 + 1):ncol(x)]
  
  #   simulationID <- paste('class',class1,'class',class2,"GC",sep ='')
  #   halfPosterior <- hierarchical.test (trainX,rho,rope_min,rope_max,simulationID,stdUpperBound,"student",chains)
  #   halfPosterior$logPredictiveEachDset <- logPredictive(halfPosterior,testX, rho)
  #   logPredictiveValues$Gc <- sum (halfPosterior$logPredictiveEachDset)
  #   
  #   simulationID <- paste('class',class1,'class',class2,"Gaussian",sep ='')
  #   halfPosteriorGauss <- hierarchical.test (trainX,rho,rope_min,rope_max,simulationID,stdUpperBound,"gaussian",chains)
  #   halfPosteriorGauss$logPredictiveEachDset <- logPredictive(halfPosteriorGauss,testX, rho)
  #   logPredictiveValues$Gauss <- sum (halfPosteriorGauss$logPredictiveEachDset)
  #   
  #   simulationID <- paste('class',class1,'class',class2,"Kruschke",sep ='')
  #   halfPosteriorKru <- hierarchical.test (trainX,rho,rope_min,rope_max,simulationID,stdUpperBound,"studentKruschke",chains)
  #   halfPosteriorKru$logPredictiveEachDset <- logPredictive(halfPosteriorKru,testX, rho)
  #   logPredictiveValues$Kru <- sum (halfPosteriorKru$logPredictiveEachDset)
  #   
  #   simulationID <- paste('class',class1,'class',class2,"Juanez",sep ='')
  #   halfPosteriorJua <- hierarchical.test (trainX,rho,rope_min,rope_max,simulationID,stdUpperBound,"studentJuanez",chains)
  #   halfPosteriorJua$logPredictiveEachDset <- logPredictive(halfPosteriorJua,testX, rho)
  #   logPredictiveValues$Jua <- sum (halfPosteriorJua$logPredictiveEachDset)
  
  
  
  
  
  
  
  #bayes factor as marginal lik of the proposed model over marg lik of alternative models
  #   bayesFactor <- list()
  #   bayesFactor$Kru <- exp (logPredictiveValues$Gc - logPredictiveValues$Kru)
  #   bayesFactor$Jua <- exp (logPredictiveValues$Gc - logPredictiveValues$Jua)
  #   bayesFactor$Gauss <- exp (logPredictiveValues$Gc - logPredictiveValues$Gauss)
  
  fileName <- paste('Rdata/sensitivityStudentAlphaBeta',class1,class2,'.Rdata', sep='')
  # save (halfPosteriorKru, halfPosteriorJua, halfPosterior, halfPosteriorGauss, bayesFactor, hierPosteriorGauss, hierPosterior, hierPosteriorJua, hierPosteriorKru, file = fileName)  
  save (hierPosterior, hierPosteriorSens,file = fileName)  
  
  
}