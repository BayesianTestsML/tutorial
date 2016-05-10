findPars <- function(){
  
x <- seq (from=1, to=1000, by=1)
KProb <- 1/29 * exp(-(x*1)/29) 
meanKprob <- sum(KProb*x)
varKprob <- sum( KProb * ((x-meanKprob)^2) )
stdKprob <- sqrt (varKprob)

alpha_values <-  seq (from=0.1, to=10, by=0.1)

beta_values <-  alpha_values / meanKprob

varGamma <-  alpha_values / (beta_values^2)
stdGamma <- sqrt (varGamma)

mae = abs(stdKprob-stdGamma)
minIdx <- which (mae == min(mae) )

alpha <- alpha_values[minIdx]
beta <- beta_values[minIdx]



}