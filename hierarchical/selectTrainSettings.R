selectTrainSettings <- function (friedmanTypeVec){
  #selects 2/3 of the provided settings, stratified over the 3 Friedman families
  
  origIdx <- 1:length(friedmanTypeVec)
  
  idx1 <- sample(which(friedmanTypeVec==1))
  idx2 <- sample(which(friedmanTypeVec==2))
  idx3 <- sample(which(friedmanTypeVec==3))
  
  stopifnot(length(idx1)==length(idx2))
  stopifnot(length(idx3)==length(idx2))
  reducedLength= round (2/3 * length(idx1))
  
  idx <- c ( idx1[1:reducedLength], idx2[1:reducedLength], idx3[1:reducedLength])
  return(idx)
}