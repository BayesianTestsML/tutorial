hierarchical.test <- function(x, sample_file, samplingType="student", 
                              alphaBeta = list('lowerAlpha' =0.5,'upperAlpha'= 5,'lowerBeta' = 0.05,'upperBeta' = .15),
                              rho = 0.1, rope_min=-0.01, rope_max=0.01, std_upper_bound=1000, chains=8)
  

  #usually you want to use the default parameters of the function, providing only x and the sample_file.
  
  
{
                              
  
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  library(matrixcalc)
  library(matrixStats)
  library(rstan)
  #for sampling from non-standardized t distribution
  library(metRology)
  
  #------------------------------------------------------------------------------- 
  if ((max(x))>1 & rope_max < 0.02) {
    stop('value of rope_max  not compatible with scale of provided x')
  }
  sample_file <- paste('stanOut/',sample_file,sep='')
  Nsamples <- dim(x)[2]
  q <- dim(x)[1]
  sample_file <- paste(sample_file,".StanOut",sep='')
  
  
  #data rescaling, to have homogenous scale among all dsets
  stdX <- mean(rowSds(x)) #we scale all the data by the mean of the standard deviation of data sets
  x <- x/stdX
  rope_min <- rope_min/stdX
  rope_max <- rope_max/stdX
  
  #search for data sets with 0 variance, which sometimes are there.
  zeroVarIdx <- which (rowSds(x) == 0)
  if ( length(zeroVarIdx) > 0) {
    #to each dset with zero variance we add a gausian noise with mean 0 and very small std dev
    #this way we preserve the mean of the original data while obtaining a positive std dev
    for (i in 1:length(zeroVarIdx)){
      noise <- runif(Nsamples/2,rope_min,rope_max)
      x[zeroVarIdx[i],1:(Nsamples/2)] <- x[zeroVarIdx[i],1:(Nsamples/2)] + noise;
      x[zeroVarIdx[i],(Nsamples/2+1):Nsamples] <- x[zeroVarIdx[i],(Nsamples/2+1):Nsamples] - noise;
    } 
  }
  
  if (q>1) {
    std_among = sd(rowMeans(x))
  } else {
    #to manage the particular case q=1
    std_among = mean(rowSds(x))
  }
  
  std_within <- mean(rowSds(x))
  
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
  
  
  #this calls the Student with learnable dofs, which is the default
  if (samplingType=="student") {
    stanfit <-  stan(file = 'stan/hierarchical-t-test.stan', data = dataList,sample_file=sample_file, chains=chains)
  }
  
  #this calls the Student with priors on the dof as in Kruschke, Gamma(1,0.0345). This is referred to a shifted exponential in his paper,
  #this Gamma is however very close to  it.
  else if (samplingType=="studentKruschke") {
    stanfit <-  stan(file = 'stan/hierarchical-t-test_nuKru.stan', data = dataList,sample_file=sample_file, chains=chains)
  }
  #this calls the Student with priors Ga (2,0.1) on the dof as in Juanez and Steel
  else if (samplingType=="studentJuanez") {
    stanfit <-  stan(file = 'stan/hierarchical-t-test_nuJuaSteel.stan', data = dataList,sample_file=sample_file, chains=chains)
  }
  
  #this calls the Gaussian
  else if (samplingType=="gaussian") {
    stanfit <-  stan(file = 'stan/hierarchical-t-testGaussian.stan', data = dataList,sample_file=sample_file, chains=chains)
  }
  
  stanResults<- extract(stanfit, permuted = TRUE)
  
  
  #get for each data set the probability of left, rope and right
  prob_right_each_dset<-vector(length = q, mode = "double")
  prob_rope_each_dset<-vector(length = q, mode = "double")
  prob_left_each_dset<-vector(length = q, mode = "double")
  delta_each_dset<-vector(length = q, mode = "double")
  
  #results on non-std data, but the computation is ok because rope_min and rope_max are already scaled.
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
  
  #compute the probability that the most probable outcome for 
  #the next delta is rope, left or right   
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
  
  
  
  results = list (
                  "nextDelta"=list("right"=probRightNextDelta, "left"=probLeftNextDelta, "rope"=probRopeNextDelta, "positive"=probPositiveNextDelta,"negative"=probNegativeNextDelta),
                  "meanDeltaEachDset"=delta_each_dset,
                  "probEachDset"=list("left"=prob_left_each_dset, "rope"=prob_rope_each_dset, 
                               "right"=prob_right_each_dset),"rawResults" = stanResults, "x"=x, "stdX"=stdX)
  
  return (results)
  
}

