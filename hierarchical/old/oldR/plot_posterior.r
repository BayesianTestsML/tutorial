# setwd("C:/Users/mangili/Dropbox/jmlr-bayesian (1)/code/R/hierarchical")
plot.posterior <- function() {
  dir <- getwd()
  load(paste(dir,"/uci_data.RData",sep=""))
  
  filename <- list(#paste(dir,"/results/uciResults_StdUpper-1000-samplingType-gaussian.Rdata",sep=""),
                   paste(dir,"/results/uciResults_StdUpper-1000-samplingType-hier.Rdata",sep="")
                   )
  
  allResults <- list()
  R <- length(filename) 
  
  for (r in 1:R) {
    load(filename[[r]])
    allResults[[r]]<-hierarchicalResults
  }
  
  col.list <- list("red","blue","black","green")
  stud.list <- c(1,1,0,0)
  
  #utiliy function which gets the average score of the given classifier on each data set
  getAvgScore <- function (currentClassifier,score) {
    avg_score <- vector(length = how_many_dsets, mode = "double");
    for ( ii in 1:how_many_dsets ){
      avg_score[ii] <- mean ( score [classifierId==currentClassifier & dsetId==dsetsList[ii]] )
    }
    return (avg_score)
  }
  
  score <- uci_classification$Percent.correct
  classifierId <- uci_classification$ClassifierID
  dsetId<- uci_classification$DatasetID
  dsetsList <- unique(dsetId);
  how_many_dsets <- length(dsetsList)
  classifierI<-double()
  classifierJ<-double()
  
  counter <-0
  results <- double()

  for (i in 1:4) {
    for (j in (i+1):5 ){

      show(c(i,j))
      counter <- counter+1
      
      # Plot data
      classifierI[counter] <- i
      classifierJ[counter] <- j
      avgScoreI <- getAvgScore (i,score)         
      avgScoreJ <- getAvgScore (j,score)
      avgScoreD <- avgScoreI-avgScoreJ
      h0 <- hist(avgScoreD, plot=FALSE,nclass=25)
      p.range.left <- h0$breaks[1] 
      p.range.right <- h0$breaks[length(h0$breaks)]  
      pdiff <- list(h0)
      ymax <- max(h0$density)
      prob.lr <- double()
      for (r in 1:R) {
        # Plot posterior
        stan <- allResults[[r]][[counter]]$stanResults
        stdx <- allResults[[r]][[counter]]$stdX
        stud <- stud.list[r]
        post.samp <- double()
        mu <- stan$mu0
        prob.lr <- c(prob.lr,mean(mu<0),mean(mu>0))
        std <- stan$std0
        if (stud==1) nu <- stan$nu
        for (k in 1:length(mu)) {
          if (stud==0) {
            post.samp[k] <- rnorm(1,mu[k],std[k])
          } else { 
            post.samp[k] <- rt(1,nu[k])*std[k] + mu[k]; 
          }
        }
        post.samp<-post.samp*stdx
        hdiff <- hist(post.samp[which((post.samp> p.range.left)&(post.samp< p.range.right))],plot=FALSE,nclass=50)
        pdiff[[r]] <- list(mu=hdiff$mids, p=hdiff$density)
        show(c(mean(stan$nu), sd(stan$nu)))
        # lines(pdiff$mu,pdiff$p,col=col.list[[r]])
      }
      ymax <-max(ymax,max(pdiff[[1]]$p))
      plotname <- paste(dir,"/results/posterior_comparison_",i,j,".pdf",sep="")
      pdf(plotname,width=6.5,height=6)
      par(mar=c(5,5,2,2)+0.1)
      plot(h0,freq = FALSE,ylim=c(0,ymax),main=NULL,xlab="x", ylab="P(x)",cex.axis=1.5, cex.lab=1.5)
      lines(pdiff[[1]]$mu,pdiff[[1]]$p,col=col.list[[1]],lwd=2)
      # lines(pdiff[[2]]$mu,pdiff[[2]]$p,lty=2, col=col.list[[2]])
      dev.off()
      results <- rbind(results,prob.lr)
    }
  }
  write.csv(results,paste(dir,"/results/prob_left_right.csv",sep=""))
}

