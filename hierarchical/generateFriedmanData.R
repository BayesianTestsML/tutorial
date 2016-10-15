generateFriedmanData <- function (friedmanType,settings){
  library(mlbench)
  
  if (friedmanType==1) {
    data <- mlbench.friedman1(settings$sampleSize, sd=settings$friedmanSd)
  }
  else if (friedmanType==2) {
    data <- mlbench.friedman2(settings$sampleSize, sd=settings$friedmanSd)
  }
  else if (friedmanType==3) {
    data <- mlbench.friedman3(settings$sampleSize, sd=settings$friedmanSd)
  }
  
  #add noisy features, generating first random values
  #from a standard normal and then reshaping them
  if (settings$redundantFeats > 0) {
    d <- rnorm(settings$sampleSize * settings$redundantFeats)
    dMatrix <- matrix(d, nrow=settings$sampleSize, ncol = settings$redundantFeats)
    data$x <- cbind( data$x, dMatrix) 
  }

  #generate the class variable by discretizing on the threshold
  data$class <- as.factor ( data$y > settings$threshold)
  
  return ( data )
  
}