compTable <- function() {
  source('Utils.R')
  suffix <- c('12','13','14','15','23','24','25','34','35','45')
  KLmatrix <- matrix(ncol=4,nrow=10)
  colnames(KLmatrix)<-c('klGauss','klJua','klHier','klKru')
  idx <- c(2,3,1)#to pick left, rope, right probabilities from the structures
  
  for (i in 1:length(suffix)){
    # fileName <- paste('sensitivityNormalStudent',suffix[i],'.Rdata',sep='')
    # fileName <- paste('sensitivityStudentAlphaBeta',suffix[i],'.Rdata',sep='')
    fileName <- paste('Rdata/sensitivityNovelAlphaBeta',suffix[i],'.Rdata',sep='')
    load(file = fileName)
    currentHierDataFrame <- as.numeric(hierPosteriorNovel$nextDelta)[idx]
    # currentKruDataFrame <- as.numeric(hierPosteriorKru$nextDelta)[idx]
    # currentJuaDataFrame <- as.numeric(hierPosteriorJua$nextDelta)[idx]
    
    
#     tmp <- matrix(nrow=10,ncol=4)
#     for (j in 1:10){
#       tmp[j,1]=KLPostShrinkage(hierPosteriorGauss,hierPosteriorGauss$delta_each_dset)[1,2] 
#       tmp[j,2]=KLPostShrinkage(hierPosteriorJua,hierPosteriorJua$delta_each_dset)[1,2]
#       tmp[j,3]=KLPostShrinkage(hierPosterior,hierPosterior$delta_each_dset)[1,2]
#       tmp[j,4]=KLPostShrinkage(hierPosteriorKru,hierPosteriorKru$delta_each_dset)[1,2]
#     }
#     KLmatrix[i,] <- apply(tmp,MARGIN = 2, median)
    
    # plotPosteriorGGplot2(hierPosterior, hierPosteriorKru, hierPosteriorJua, suffix[i])
    
    #     currentDataFrame <-
    #       as.data.frame(cbind(halfPosteriorGauss$nextDelta, halfPosteriorJua$nextDelta,halfPosterior$nextDelta,
    #                           halfPosteriorKru$nextDelta))
    #     logPredictive <- c (sum(halfPosteriorGauss$logPredictiveEachDset),sum(halfPosteriorJua$logPredictiveEachDset),
    #                         sum(halfPosterior$logPredictiveEachDset), sum(halfPosteriorKru$logPredictiveEachDset))
    #     currentDataFrame <- rbind(logPredictive,currentDataFrame)
    if (i==1){
      hierDataFrame <- currentHierDataFrame
      # kruDataFrame <- currentKruDataFrame
      # juaDataFrame <- currentJuaDataFrame
      }
    else{
      hierDataFrame <- rbind(hierDataFrame,currentHierDataFrame)
      # kruDataFrame <- rbind(kruDataFrame,currentKruDataFrame)
      # juaDataFrame <- rbind(juaDataFrame,currentJuaDataFrame)
    }
  }
    # colnames (dataFrame) <- c('gauss','jua','gc','kru')
    colnames(hierDataFrame) <- c('left','rope','right')
    # colnames(juaDataFrame) <- c('left','rope','right')
    # colnames(kruDataFrame) <- c('left','rope','right')
    rownames(hierDataFrame)<- suffix
    # rownames(juaDataFrame)<- suffix
    # rownames(kruDataFrame)<- suffix
    
    # df <- as.data.frame(cbind(hierDataFrame,kruDataFrame,juaDataFrame))
    df <- as.data.frame(cbind(hierDataFrame))
     roundedDF <- round(df,digits=2)
    # write.table(roundedDF, file='csvResults/PredictionsHierKruJua.csv', sep=',')
    write.table(roundedDF, file='csvResults/PredictionsHierAlpha05-5-Beta05-015.csv', sep=',')
    # write.table(KLmatrix, file='csvResults/KLNormalJuaKruHier.csv', sep=',', col.names = TRUE)
    
  }