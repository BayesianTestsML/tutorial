multipleDsetsRopeSignedRank <- function(delta0,std0,sampleSizesSampling,reps,simulation_ID, delta_acc_sampling='cauchy') {
  #let signed-rank and hierarchical test compete on simulated data sets
  #allows for fixed or variable sample sizes
  #allows for mixture or student-generated delta_i
  
  library(MASS)
  library(matrixStats)
  library(rstan)
  library(loo)
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  
  source('cv_nbc.R')
  source('hierarchical_test.R')
  
  if  (( delta_acc_sampling!='cauchy') && (delta_acc_sampling!='mixture') && (delta_acc_sampling!='gaussian') ) { 
    stop('wrong delta_acc_sampling')
  }
  
  #debug
  how_many_dsets=3;
  #how_many_dsets=50;
  
  #check arguments
  std_upper_bound=1000
  
  
  # this control the variance of the simulation.
  # one classifier has accuracy theta_star, the other theta_star=delta_acc.
  # with theta_star=1, we have zero variance.
  # witht theta_star=0.5, we have maximum variance
  theta_star <- 0.9
  n_folds <- 10
  n_runs <- 10
  
  if (sampleSizesSampling=="fixed")
    sample_size=1000
  
  
  file_str0 <- paste('resultsDelta0',delta0,'Std0',std0,'_','sample_sizes',sample_size,'_',delta_acc_sampling,'_simulID_',simulation_ID,sep = "")
  filename <- paste(file_str0,'.csv',sep = "")
  log_filename <- paste('log_',file_str0,'.csv',sep = "")
  
  rope_min <- -0.005;
  rope_max <- 0.005;
  #set the seed for reproducible results
  set.seed(simulation_ID)
  
  
  
  
  #point estimator of the correlation
  rho <- 1/n_folds
  sign_rank_p_value <- matrix(0,reps)
  prob_left_hier_nextDelta <- matrix(0,reps)
  prob_rope_hier_nextDelta <- matrix(0,reps)
  prob_right_hier_nextDelta <- matrix(0,reps)
  logLossHier <- matrix(0,reps)
  mean.diff<- matrix(0,reps)
  var.diff<- matrix(0,reps)
  settings <- list() 
  results <- list() 
  counter <- 1
  futureDsets<-3000
  
  
  current_dset_size <- sample_size;
  test_set_size <- current_dset_size/n_folds;
  #to be appended to various created files
  file_str <- paste('delta0',delta0,'Std0',std0,'_dsets_',how_many_dsets,'_delta_acc_',delta_acc_sampling,'_sampleSizeSampling_',sampleSizesSampling,'simulationID_',simulation_ID,sep="");
  
  
  for (k in 1:reps) {
    
    #if (k%%100==0) cat ('reps',k,'\n')
    
    if (delta_acc_sampling=='cauchy'){
      delta_acc_each_dset <-  rt(how_many_dsets,1)*std0 + delta0
      futureDeltaAcc <- rt(futureDsets,1)*std0 + delta0
    }
    else if (delta_acc_sampling=='gaussian'){
      delta_acc_each_dset <-  rt(how_many_dsets,1)*std0 + delta0
      futureDeltaAcc <- rt(futureDsets,100)*std0 + delta0
    }
    else if (delta_acc_sampling=='mixture') {
      #random number from the uniform and sorted
      idx <- sort(runif(how_many_dsets+futureDsets))
      std <- 0.01
      mean1 <- -0.01
      mean2 <- 0.03
      #initialization of tmp
      tmp <-  idx 
      for (kk in 1:length(idx)) {
        if (runif(1) > 0.5)
          tmp[kk]<-rnorm(1, mean = mean1, sd = std)
        else
          tmp[kk]<-rnorm(1, mean = mean2, sd = std)
      }
      delta_acc_each_dset <- tmp (1:how_many_dsets) 
      futureDeltaAcc <- tmp((how_many_dsets+1):how_many_dsets+futureDsets)
    }
    futureLeftProb=sum(futureDeltaAcc<rope_min)/futureDsets
    futureRightProb=sum(futureDeltaAcc>rope_max)/futureDsets
    futureRopeProb=sum(futureDeltaAcc<rope_max & futureDeltaAcc>rope_min)/futureDsets
    
    
    sample_diff_acc_a_b_each_dset <- matrix(0,how_many_dsets);
    Diff_ab <- matrix(0,n_runs*n_folds,how_many_dsets);
    
    for (i in 1:how_many_dsets){
      
      current_delta_acc <- delta_acc_each_dset[i];
      cv_results <- list('delta'=rep(0,n_runs*n_folds));
      
      while (var(cv_results$delta)==0 && mean(cv_results$delta)==0) { 
        cv_results <- cv_competing_nbc(n_runs,n_folds,current_dset_size,current_delta_acc,theta_star);
        diff_a_b <- cv_results$delta;
      }
      
      sample_diff_acc_a_b_each_dset[i] <- mean(diff_a_b);
      
      
      #equivalent to p-value one-sided right tailed
      Diff_ab[,i]=diff_a_b;
    }
    
    #frequentist Wilcoxon signed rank
    sign_rank_p_value[k]  <-  wilcox.test(sample_diff_acc_a_b_each_dset,alternative = "two.sided",exact=FALSE)$p.value;
    # sign_rank_p_value[k]  <-  t.test(sample_diff_acc_a_b_each_dset)$p.value;
    mean.diff[k] <- mean(delta_acc_each_dset)
    var.diff[k] <- var(delta_acc_each_dset)
    
    stan_prob <- hierarchical.test(t(Diff_ab),rho,rope_min,rope_max,file_str,std_upper_bound,'student')
    prob_left_hier_nextDelta[k]  <- stan_prob$nextDelta$left
    prob_rope_hier_nextDelta[k]  <- stan_prob$nextDelta$rope
    prob_right_hier_nextDelta[k] <- stan_prob$nextDelta$right
    
    #to be implemented
    #logLossSignedRank[k] <- logLossSignedRank[k]/how_many_dsets
    logLossHier[k] <-  futureRightProb * prob_right_hier_nextDelta[k] + futureLeftProb * prob_left_hier_nextDelta[k] +  futureRopeProb * prob_rope_hier_nextDelta[k]
  }
  #close the loop on the repetitions
  
  
  
  results <- list('how_many_dsets'=how_many_dsets,
                  'sample_size'=current_dset_size,
                  'delta0'=delta0,
                  'std0'=std0,
                  'sign_rank_p_value'=sign_rank_p_value,
                  'prob_left_hier_nextDelta'=prob_left_hier_nextDelta,
                  'prob_rope_hier_nextDelta'=prob_rope_hier_nextDelta,
                  'prob_right_hier_nextDelta'=prob_right_hier_nextDelta,
                  'logLossHier'=logLossHier)
  
  save_results <- function(results,filename) {
    
    mystring <- paste('delta_acc_sampling,delta0,std0,num_experiments,',
                      'how_many_dsets,sample_size,signRankPower,signRankPValue,hier_left_95,hier_rope_95,hier_right95,',
                      'hier_left_90,hier_rope_90,hier_right90,median_hier_p_left,median_hier_p_rope,',
                      'median_hier_p_right,logLossHier,',sep="");
    write(mystring, filename, append = FALSE)
    
    num_experiments <- reps;
    h1_sign_rank <- mean(results$sign_rank_p_value<.05)
    h1_sign_rankPvalue <- median(results$sign_rank_p_value)
    hier_left_95 <- mean(results$prob_left_hier_nextDelta>.95)
    hier_rope_95 <- mean(results$prob_rope_hier_nextDelta>.95)
    hier_right_95 <- mean(results$prob_right_hier_nextDelta>.95)
    hier_left_90 <- mean(results$prob_left_hier_nextDelta>.90)
    hier_rope_90 <- mean(results$prob_rope_hier_nextDelta>.90)
    hier_right_90 <- mean(results$prob_right_hier_nextDelta>.90)
    
    tmp_vector <- c(h1_sign_rank, h1_sign_rankPvalue, hier_left_95,hier_rope_95, hier_right_95,hier_left_90,hier_rope_90,hier_right_90,
                    median(results$prob_left_hier_nextDelta),median(results$prob_rope_hier_nextDelta),median(results$prob_right_hier_nextDelta),
                    mean(results$logLossHier) );
    
    mystring <- paste(delta_acc_sampling,',',results$delta0,',',
                      results$std0,',',num_experiments,',',results$how_many_dsets,
                      ',',results$sample_size,',',paste(tmp_vector,collapse=","));
    write(mystring, filename, append = TRUE)
  }
  
  
  #save csv file
  save_results(results, filename)
  
  #save the Rdata file
  rdata_filename <- paste(file_str0,'.Rdata',sep="")
  save(results, file = rdata_filename)
  
  return(results)
}








