batchSensitivity <- function (){
  library('MASS')
  source ('sensitivityNormalStudent.R')
  for (i in 1:4) {
    for (j in (i+1) : 5){
      cat (i,j)
      sensitivityNormalStudent(i,j)
    }
  }
}