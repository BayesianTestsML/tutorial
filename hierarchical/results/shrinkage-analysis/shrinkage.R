shrinkage <- function(q=20,generation="mixture"){
  
  library(rstan)
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  
  #generate the data
  std <- 1 #set a large std dev for the groups
  Nsamples <- 10
  totReps <- 10
  rmse_maxLik <- vector(length = totReps)
  rmse_shrinkage <-   vector(length = totReps)
  distance_maxLik <- vector(length = totReps)
  #how much the max lik estimates are pulled towards the general mean
  shrinkage <-   vector(length = totReps)
  
  for (rep in 1:totReps) {
    #generate according to the mixture
    pi <- 0.5
    Msd <- 0.01
    Mmean1 <- 0.05
    Mmean2 <- -0.05
    meanGroups<-c(rnorm(q*pi,mean=Mmean1,sd=Msd),rnorm(q*(1-pi),mean=Mmean2,sd=Msd))
    x <- matrix (ncol = q, nrow = Nsamples)
    mixtureVar <- Msd^2 + pi*Mmean1^2 + (1-pi)*Mmean2^2 - (pi*Mmean1+(1-pi)*Mmean2)^2
    
    
    #overwrite sampling from the normal
    if (generation == "normal"){
      normalMean <- mean(c(Mmean1,Mmean2))
      normalStd  <- sqrt(mixtureVar)
      meanGroups<-rnorm(q,mean=normalMean,sd=normalStd)
    }
    
    for (j in 1:q){
      x[,j] <- rnorm(Nsamples,mean = meanGroups[j], sd = std)
    }
    
    dataList = list(
      Nsamples = Nsamples,
      q = q ,
      x = x,
      std = std
    )
    
    sample_file <- "shrinkage-sample"
    fit <-  stan(file = 'shrinkage.stan', data = dataList, sample_file=sample_file)  
    stanResults<- extract(fit, permuted = TRUE)
    
    shrinkageMeans <- vector(length = q)
    for (i in 1:q){
      shrinkageMeans[i] <- mean(stanResults$mu[,i])
    }
    
    maxLikMeans <- colMeans(x)
    generalMean <- mean(maxLikMeans)
    
    distance_maxLik[rep] <-0
    shrinkage[rep] <- 0
    
    for (i in 1:q){
      distance_maxLik[rep] <- distance_maxLik[rep] + abs (maxLikMeans[i]-Mmean1)
      if (generalMean>maxLikMeans[i]){
        shrinkage[rep] <- shrinkage[rep] + shrinkageMeans[i]-maxLikMeans[i]
      }
      else{
        shrinkage[rep] <- shrinkage[rep] + maxLikMeans[i]-shrinkageMeans[i]
      }
    }
    
    rmse_maxLik[rep] <- sqrt(mean ((maxLikMeans - meanGroups)^2)  )
    rmse_shrinkage[rep] <-  sqrt(mean ((shrinkageMeans - meanGroups)^2)  )
    
    mu0 <- stanResults$mu0
    std0 <- stanResults$std0
    post.samp <-vector(length = length(mu0))
    #plot the posterior
    up <- length(mu0)
    for (k in 1:up) {
      post.samp[k] <- rnorm(1,mu0[k],std0[k]) 
    }
    
    
  }
  
  results <- list('rmse_maxLik'=rmse_maxLik,'rmse_shrinkage'=rmse_shrinkage,"shrinkage"=shrinkage)
  return (results)
}


