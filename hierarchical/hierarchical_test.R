hierarchical.test <- function(x, sample_file, samplingType="student", 
                              alphaBeta = list('lowerAlpha' =1,'upperAlpha'= 2,'lowerBeta' = 0.01,'upperBeta' = 0.1),
                              rho = 0.1, rope_min=-0.01, rope_max=0.01, std_upper_bound=1000,chains=8)
{
                              
  
  # rstan_options(auto_write = TRUE)
  # options(mc.cores = parallel::detectCores())
  library(matrixcalc)
  library(matrixStats)
  library(rstan)
  library(coda) #for computing the HDP interval
  #for sampling from non-standardized topt
  library(metRology)
  
  #------------------------------------------------------------------------------- 
  if ((max(x))>1 & rope_max < 0.02) {
    stop('value of rope_max  not compatible with scale of provided x')
  }
  sample_file <- paste('stanOut/',sample_file,sep='')
  
  #we always standardize the received data 
  standardization=1
  
  Nsamples <- dim(x)[2]
  q <- dim(x)[1]
  sample_file <- paste(sample_file,".StanOut",sep='')
  
  
  
  #data rescaling
  if (standardization==1){
    stdX <- mean(rowSds(x)) #we scale all the data by the mean of the standard deviation of data sets
    x <- x/stdX
    rope_min <- rope_min/stdX
    rope_max <- rope_max/stdX
  }
  
  #search for data sets with 0 variance
  zeroVarIdx <- which (rowSds(x) == 0)
  if ( length(zeroVarIdx) > 0) {
    #to each dset with zero variance we add a gausian noise with mean 0 and  very low std dev
    #this way we preserve the mean of the original data while obtaining a positive std dev
    #minStd <- (rope_max - rope_min)/(6) // old gaussian noise, but we need more variance
    #the noise is uniformly drawn from the rope
    for (i in 1:length(zeroVarIdx)){
      #       noise <- rnorm(Nsamples/2,0,minStd)
      noise <- runif(Nsamples/2,rope_min,rope_max)
      x[zeroVarIdx[i],1:(Nsamples/2)] <- x[zeroVarIdx[i],1:(Nsamples/2)] + noise;
      x[zeroVarIdx[i],(Nsamples/2+1):Nsamples] <- x[zeroVarIdx[i],(Nsamples/2+1):Nsamples] - noise;
    } 
  }
  
  #   #compute the max likelihood variance, to be compared with the estimated one
  #   meanMLE <-rowMeans(x)
  #   beta <- -rho*(1-rho)^(Nsamples-2);
  #   invM <- matrix(beta,nrow=Nsamples,ncol = Nsamples)
  #   for (i in 1:Nsamples){
  #     invM[i,i]<- (1 + (Nsamples-2)*rho) * (1-rho)^(Nsamples-2)
  #   }
  #   detM <- (1+(Nsamples-1)*rho)*(1-rho)^(Nsamples-1);
  #   invM <-invM/detM
  #   
  #   sigmaHat <- vector(length = q)
  #   for (i in 1:q)
  #   {
  #     currentMean<-meanMLE[i]
  #     tmpVec <- as.vector(x[i,]-currentMean)
  #     Z<- tmpVec %*% t(tmpVec)
  #     sigmaHat[i]<-sqrt(matrix.trace(invM %*% Z)/Nsamples)
  #   }
  #   sigmaMatrix<-data.frame(
  #     sigmaHat=sigmaHat
  #   )
  #   meanMatrix<-data.frame(
  #     meanMLE=meanMLE
  #   )
  #save sigmaHat to file
  # filename<-'sigmaHat'
  #   csv_filename <- paste (filename,"csv",sep=".")
  #   write.matrix(sigmaMatrix, file = csv_filename, sep = ",")
  #sigma has been estimated and saved
  
  #   # save meanMle to file
  #   filename<-'meanMLE'
  #   csv_filename <- paste (filename,"csv",sep=".")
  #   write.matrix(meanMatrix, file = csv_filename, sep = ",")
  
  
  
  if (q>1) {
    std_among = sd(rowMeans(x))
  } else {
    #to manage the particular case q=1
    std_among = mean(rowSds(x))
  }
  
  std_within <- mean(rowSds(x))
  
  #notice the lower bound to 0
  dataList = list(
    deltaLow = -max(abs(x)),
    deltaHi = max(abs(x)),
    stdLow = 0,
    stdHi = std_within*std_upper_bound,
    std0Low = 0,
    std0Hi = std_among*std_upper_bound,
    Nsamples = Nsamples,
    q = q ,
    x = x ,
    rho = rho,
    upperAlpha = alphaBeta$upperAlpha,
    lowerAlpha = alphaBeta$lowerAlpha,
    upperBeta = alphaBeta$upperBeta,
    lowerBeta = alphaBeta$lowerBeta
  )
  
  
  
  startTime<-proc.time()
  
  #this calls the Student with learnable dofs 
  if (samplingType=="student") {
    stanfit <-  stan(file = 'stan/hierarchical-t-test.stan', data = dataList,sample_file=sample_file, chains=chains)
  }
  
  #this calls the Gaussian
  else if (samplingType=="studentKruschke") {
    stanfit <-  stan(file = 'stan/hierarchical-t-test_nuKru.stan', data = dataList,sample_file=sample_file, chains=chains)
  }
  #this calls the Gaussian
  else if (samplingType=="studentJuanez") {
    stanfit <-  stan(file = 'stan/hierarchical-t-test_nuJuaSteel.stan', data = dataList,sample_file=sample_file, chains=chains)
  }
  
  #this calls the Gaussian
  else if (samplingType=="gaussian") {
    stanfit <-  stan(file = 'stan/hierarchical-t-testGaussian.stan', data = dataList,sample_file=sample_file, chains=chains)
  }
  
  
  
  #estimate of the posterior variance for comparison purposes
  #     posteriorSigma<-vector(length = q)
  #     for (i in 1:q) {
  #       posteriorSigma[i]<-median(stanResult$sigma[i])
  #     }
  # }
  
  #   if (varianceModel=="fixed"){
  #     dataList$sigmaHat <- sigmaHat
  #     
  #     #this calls the Student with learnable dofs 
  #     if (samplingType=="student") {
  #       stanfit <-  stan(file = 'stan/hierarchical-t-test-fixedSigma.stan', data = dataList,sample_file=sample_file, chains=10)
  #     }
  #     
  #     #this calls the Gaussian, not yet implemented
  #     else if (samplingType=="gaussian") {
  #       stanfit <-  stan(file = 'stan/hierarchical-t-testGaussian-fixedSigma.stan', data = dataList,sample_file=sample_file, chains=4)
  #     }
  #   }
  #   
  
  
  
  stanResults<- extract(stanfit, permuted = TRUE)
  stopTime<-proc.time()
  elapsed=stopTime - startTime
  show(elapsed)
  
  
  
  #get for each data set the probability of left, rope and right
  prob_right_each_dset<-vector(length = q, mode = "double")
  prob_rope_each_dset<-vector(length = q, mode = "double")
  prob_left_each_dset<-vector(length = q, mode = "double")
  
  delta0<-stanResults$delta0
  postSamples <- length(delta0)  
  prob_right_delta0<-mean(delta0>rope_max)
  prob_left_delta0<-mean(delta0<rope_min)
  prob_rope_delta0<-mean(delta0>rope_min & delta0<rope_max)
  prob_positive_delta0 <- mean(delta0>0)
  prob_negative_delta0 <- mean(delta0<0)
  
  delta_each_dset<-vector(length = q, mode = "double")
  #results on non-std data
  sampled_delta_each_dset<-stanResults$delta
  for (j in 1:q){
    prob_right_each_dset[j] <- mean(sampled_delta_each_dset[,j]>rope_max)
    prob_rope_each_dset[j]  <- mean(sampled_delta_each_dset[,j]>rope_min & sampled_delta_each_dset[,j]<rope_max)
    prob_left_each_dset[j]  <- mean(sampled_delta_each_dset[,j]<rope_min)
    delta_each_dset[j] <- mean(sampled_delta_each_dset[,j])*stdX
  }
  
  
  
  #keep small the data to be saved by removing helping variables
  stanResults$diff<-NULL
  stanResults$diagQuad<-NULL
  stanResults$oneOverSigma2<-NULL
  stanResults$nuMinusOne<-NULL
  stanResults$log_lik<-NULL
  
  #compute the probability of delta(q+1) being within the rope, by sampling   
  postSamples <- length(stanResults$delta0)
  sampledRopeWins <- 0 
  sampledLeftWins <- 0 
  sampledRigthWins <- 0 
  sampledPositiveWins <- 0
  sampledNegativeWins <- 0
  cumulativeRope <- vector (length = postSamples)
  cumulativeRight <- vector (length = postSamples)
  cumulativeLeft <- vector (length = postSamples)
  
  std <- stanResults$std0
  mu  <- stanResults$delta0
  #for all the student-based model we sample in the same way
  for (r in 1:postSamples){
    if ( any (samplingType == c('student', 'studentKruschke', 'studentJuanez'))) {
      nu  <- stanResults$nu
      cumulativeRope[r] <- pt.scaled(rope_max, df=nu[r], mean=mu[r], sd=std[r]) - pt.scaled(rope_min, df=nu[r], mean=mu[r], sd=std[r])
      cumulativeLeft[r] <- pt.scaled(rope_min, df=nu[r], mean=mu[r], sd=std[r])
      cumulativeRight[r] <- 1-pt.scaled(rope_max, df=nu[r], mean=mu[r], sd=std[r])
    }
    else if (samplingType=="gaussian") {
      cumulativeRope[r] <- pnorm(rope_max, mean=mu[r], sd=std[r]) - pnorm(rope_min, mean=mu[r], sd=std[r])
      cumulativeLeft[r] <- pnorm(rope_min,  mean=mu[r], sd=std[r])
      cumulativeRight[r] <- 1-pnorm(rope_max, mean=mu[r], sd=std[r])
    }
    
    if (cumulativeRope[r] > cumulativeLeft[r] & cumulativeRope[r] > cumulativeRight[r]){
      sampledRopeWins <- sampledRopeWins + 1
    }
    else if  (cumulativeLeft[r] > cumulativeRope[r] & cumulativeLeft[r] > cumulativeRight[r]){
      sampledLeftWins <- sampledLeftWins + 1
    }
    else {
      sampledRigthWins <- sampledRigthWins +1
    }
    if (mu[r]>0){
      sampledPositiveWins <- sampledPositiveWins + 1
    }
    else {
      sampledNegativeWins <- sampledNegativeWins + 1
    }
  }
  
  
  
  
  probRightNextDelta <- sampledRigthWins/(sampledRigthWins+sampledLeftWins+sampledRopeWins)
  probLeftNextDelta  <- sampledLeftWins/(sampledRigthWins+sampledLeftWins+sampledRopeWins)
  probRopeNextDelta  <- sampledRopeWins/(sampledRigthWins+sampledLeftWins+sampledRopeWins)
  probPositiveNextDelta  <- sampledPositiveWins/(sampledPositiveWins+sampledNegativeWins)
  probNegativeNextDelta  <- sampledNegativeWins /(sampledPositiveWins+sampledNegativeWins)
  
HDP <- coda::HPDinterval(coda::as.mcmc(as.vector(stanResults$delta0 * stdX)))
  
  
  results = list ("delta0"=list("right"=prob_right_delta0, "left"=prob_left_delta0, "rope"=prob_rope_delta0, "positive"=prob_positive_delta0,"negative"=prob_negative_delta0,
                                "HDPlower"= HDP[1], "HDPupper"=HDP[2]),
                  "nextDelta"=list("right"=probRightNextDelta, "left"=probLeftNextDelta, "rope"=probRopeNextDelta, "positive"=probPositiveNextDelta,"negative"=probNegativeNextDelta),
                  "delta_each_dset"=delta_each_dset,
                  "delta"=list("left"=prob_left_each_dset, "rope"=prob_rope_each_dset, 
                               "right"=prob_right_each_dset),"stanResults" = stanResults, "x"=x, "stdX"=stdX,
                  "nu"=stanResults$nu,"rho"=stanResults$rho,"std0"=stanResults$std0)
  
  return (results)
  
}

