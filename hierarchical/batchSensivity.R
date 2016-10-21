batchSensitivity <- function (){
  library('MASS')
  source ('sensitivityNormalStudent.R')
  for (i in 1:4) {
    for (j in (i+1) : 5){
      cat (i,j)
      sensitivityNormalStudent(i,j)
#       currentLogPredictiveValues <- sensitivityNormalStudent(i,j)
#       if (i==1 & j==2)
#         logPredictiveValues <- currentLogPredictiveValues
#       else
#         logPredictiveValues <- rbind (logPredictiveValues, currentLogPredictiveValues)
    }
  }
  # csvFilename <- paste('csvResults/logPredictive',".csv",sep='')
  # write.matrix(logPredictiveValues,file=csvFilename, sep=",")
}