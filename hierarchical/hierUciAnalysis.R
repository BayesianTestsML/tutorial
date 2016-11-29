hierUciAnalysis <- function (class1, class2){
  #class1 and class2 are the two classifier being compared
  #'naive Bayes','aode','hnb','j48','j48_grafted' are coded as 1,2,3,4,5 respectively
  #'e.g. to compare naive Bayes and aode: sensitivityNormalStudent(1,2)
  #infers the hierarchical model on the results of classifiers class1 and class2, as loaded from data.
  
  library(rstan)
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  source ("hierarchical_test.R")
  source ("logPredictive.R")
  source ("Utils.R")
  
  #this workspace needs to be there
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
  
  
  
  #those lines if you want to infer the hierarchical model using Gamma prior on the degrees of freedom,
  #as sensitivity analysis
  
  # simulationID <- paste('class',class1,'class',class2,"Kruschke",sep ='')
  # hierPosteriorKru <- hierarchical.test (x=x,rho=rho,samplingType = "studentKruschke",rope_min = rope_min,
  # rope_max = rope_max,std_upper_bound = stdUpperBound,chains = chains,sample_file = simulationID)
  
  # simulationID <- paste('class',class1,'class',class2,"Juanez",sep ='')
  # hierPosteriorJua <- hierarchical.test (x=x,rho=rho,samplingType = "studentJuanez",rope_min = rope_min,
  # rope_max = rope_max,std_upper_bound = stdUpperBound,chains = chains,sample_file = simulationID)
  
  simulationID <- paste('class',class1,'class',class2,"GCsens",sep ='')
  
  
  #novel setup
  #this setup of parameters works well in most cases, providing generally better fit than 
  #a single gamma prior on nu 
  alphaBeta = list('lowerAlpha' =0.5,'upperAlpha'= 5,'lowerBeta' = .05,'upperBeta' = .15)
  hierPosterior <- hierarchical.test (x = x,sample_file = simulationID,samplingType = "student", alphaBeta = alphaBeta)
  
  fileName <- paste('Rdata/hierUciAnalysis',class1,class2,'.Rdata', sep='')
  save (hierPosterior, file = fileName)  
  
}