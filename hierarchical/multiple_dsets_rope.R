multiple_dsets_rope <- function(delta0=0,std0=0.01,how_many_dsets=50,reps=250,sample_sizes=500,simulation_ID=1, delta_acc_sampling='cauchy') {
  
  library(MASS)
  library(matrixStats)
  library(rstan)
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  
  source('cv_nbc.R')
  source('hierarchical_test.R')
  
  if  (( delta_acc_sampling!='cauchy') && (delta_acc_sampling!='gaussian') && (delta_acc_sampling!='mixture') ) { 
    stop('wrong delta_acc_sampling')
  }
  
  #set the seed for reproducible results
  set.seed(simulation_ID)
  
  # this control the variance of the simulation.
  # the first classifier has accuracy theta_star, the second theta_star+delta_acc.
  # with theta_star=1, we have zero variance.
  # with theta_star=0.5, we have maximum variance
  theta_star <- 0.9
  n_folds <- 10
  n_runs <- 10
  
  file_str0 <- paste('csvResults/Delta0',delta0,'Std0',std0,'_','sample_sizes',sample_sizes,'_',delta_acc_sampling,'_simulID_',simulation_ID,sep = "")
  filename <- paste(file_str0,'.csv',sep = "")
  log_filename <- paste('log_',file_str0,'.csv',sep = "")
  
  rope_min <- -0.01
  rope_max <- 0.01
  
  
  #point estimator of the correlation
  rho <- 1/n_folds;
  sign_rank_p_value <- matrix(0,reps)
  probLeftNextDelta <- matrix(0,reps)
  probRopeNextDelta <- matrix(0,reps)
  probRightNextDelta <- matrix(0,reps)
  rmseMleDelta_i <- matrix(0,reps)
  rmseHierDelta_i <- matrix(0,reps)
  settings <- list() 
  results <- list() 
  counter <- 1
  
  #originally the function was supposed to handle multiple value of parameters.
  #in reality we always used on a cluster running a single setting.
  for (k in 1:length(how_many_dsets)) {
    for (i in 1:length(sample_sizes)) {
      settings[[counter]] <- list('dsets'=how_many_dsets[k],'sample_size'=sample_sizes[i]); ##ok<AGROW>
      counter=counter+1;
    }
  }
  
  #simulation
  for (j in 1:length(settings)) {
    cat('setting',j,'\n');
    current_many_dsets <- settings[[j]]$dsets;
    current_dset_size <- settings[[j]]$sample_size;
    test_set_size <- current_dset_size/n_folds;
    #to be appended to various created files
    file_str <- paste('delta0',delta0,'Std0',std0,'_dsets_',current_many_dsets,'_delta_acc_',delta_acc_sampling,'_samples_sizes_',sample_sizes,'simulationID_',simulation_ID,sep="");
    
    
    for (k in 1:reps) {
      
      if (k%%100==0) cat ('reps',k,'\n');
      
      #Gaussian or cauchy sampling
      if (delta_acc_sampling=='gaussian'){
        #we use the Gaussian only to compare with the mixture, hence the hard coding of the parameters
        delta_mean <- (0.005 + 0.02) / 2
        std <- 0.01
        std_mixture <- sqrt( std^2 + (0.005^2)/2 + (0.02^2)/2 - (0.5*(0.005+0.02))^2 )
        delta_acc_each_dset <- rnorm(current_many_dsets, mean = delta_mean, sd = std_mixture)
      }
      else if (delta_acc_sampling=='cauchy'){
        delta_acc_each_dset <-  rt(current_many_dsets,1)*std0 + delta0; 
      }
      else if (delta_acc_sampling=='mixture'){
        #random number from the uniform and sorted
        idx <- sort(runif(current_many_dsets))
        count1 <- max(1,sum (idx>0.5))
        count1 <- min(count1,current_many_dsets-1)
        count2 <- how_many_dsets - count1 
        std <- 0.01
        mean1 <- 0.005 
        mean2 <- 0.02
        delta_acc_each_dset <- vector(length = current_many_dsets)
        delta_acc_each_dset[1:count1] <- rnorm(count1, mean = mean1, sd = std); 
        delta_acc_each_dset[(count1+1):(count1+count2)] <- rnorm(count2, mean = mean2, sd = std); 
      }
      
      sample_diff_acc_a_b_each_dset <- matrix(0,current_many_dsets);
      diffMatrix <- matrix(0,n_runs*n_folds,current_many_dsets);
      
      for (i in 1:current_many_dsets){
        
        current_delta_acc <- delta_acc_each_dset[i];
        cv_results <- list('delta'=rep(0,n_runs*n_folds))
        while (var(cv_results$delta)==0 && mean(cv_results$delta)==0) { 
          cv_results <- cv_competing_nbc(n_runs,n_folds,current_dset_size,current_delta_acc,theta_star);
          currentDiff <- cv_results$delta;
        }
        sample_diff_acc_a_b_each_dset[i] <- mean(currentDiff);
        diffMatrix[,i]=currentDiff;
      }
      
      #frequentist Wilcoxon signed rank
      sign_rank_p_value[k]  <-  wilcox.test(sample_diff_acc_a_b_each_dset,alternative = "two.sided",exact=FALSE)$p.value;
      
      #running Stan
      #we adopt here a reduced number of chains because artificial data are smooth and convergence is easy
      stanModel <- hierarchical.test(x=t(diffMatrix), sample_file =  file_str,samplingType =  'student',
                                     chains=4)
      
      #probability of left, right and rope being the most probable outcome on the next data set.
      probLeftNextDelta[k]  <- stanModel$nextDelta$left
      probRopeNextDelta[k]  <- stanModel$nextDelta$rope
      probRightNextDelta[k] <- stanModel$nextDelta$right
      
      
      #comparing MLE and hierarchical estimates of Delta_i
      maxLikMeans <- colMeans(diffMatrix)
      rmseMleDelta_i[k] <- sqrt(mean( (delta_acc_each_dset- maxLikMeans)^2 ))
      rmseHierDelta_i[k] <- sqrt(mean( (delta_acc_each_dset- stanModel$meanDeltaEachDset)^2 ))
      
    }#closes the loop on the repetitions
    
    
    
    results[[j]] <- list('how_many_dsets'=current_many_dsets,
                         'sample_size'=current_dset_size,
                         'delta0'=delta0,
                         'std0'=std0,
                         'sign_rank_p_value'=sign_rank_p_value,
                         'probLeftNextDelta'=probLeftNextDelta,
                         'probRopeNextDelta'=probRopeNextDelta,
                         'probRightNextDelta'=probRightNextDelta,
                         'rmseHierDelta_i'=rmseHierDelta_i,
                         'rmseMleDelta_i'=rmseMleDelta_i
                         )  
  }
  
  save_results <- function(results,filename) {
    
    mystring <- paste('delta_acc_sampling,delta0,std0,num_experiments,',
                      'how_many_dsets,sample_size,',
                      'signRankPower,signRankPValue,',
                      'hier_left_95,hier_rope_95,hier_right95,',
                      'medianHierLeft,medianHierRope,medianHierRight,',
                      'rmse_hier,rmse_mle',sep="")
    
    write(mystring, filename, append = FALSE)
    
    for (ii in 1:length(results)){
      
      tmp_vector <- c(
                      results[[ii]]$delta0,
                      results[[ii]]$std0,
                      length(results[[ii]]$sign_rank_p_value),
                      results[[ii]]$how_many_dsets,
                      results[[ii]]$sample_size,
                      mean(results[[ii]]$sign_rank_p_value<.05),
                      median(results[[ii]]$sign_rank_p_value),
                      mean(results[[ii]]$probLeftNextDelta>.95),
                      mean(results[[ii]]$probRopeNextDelta>.95),
                      mean(results[[ii]]$probRightNextDelta>.95),
                      median(results[[ii]]$probLeftNextDelta),
                      median(results[[ii]]$probRopeNextDelta),
                      median(results[[ii]]$probRightNextDelta),
                      mean(results[[ii]]$rmseHierDelta_i),
                      mean(results[[ii]]$rmseMleDelta_i)
                      )
      
      mystring <- paste(delta_acc_sampling,',',paste(tmp_vector,collapse=","))
      write(mystring, filename, append = TRUE)
    }
  }
  
  #save csv file
  save_results(results, filename)
  
  #save the Rdata file
  rdata_filename <- paste(file_str0,'.Rdata',sep="");
  save(results, file = rdata_filename)
  
}

ttest_Bayesian <- function(diff_a_b,rho,rope_min,rope_max) {
  delta <- mean(diff_a_b)
  n <- length(diff_a_b)
  df <- n-1
  stdX <- sd(diff_a_b)
  sp <- sd(diff_a_b)*sqrt(1/n + rho/(1-rho))
  p.left <- pt((rope_min - delta)/sp, df)
  p.rope <- pt((rope_max - delta)/sp, df)-p.left
  indep_p_each_dset <- list('left'=p.left,'rope'=p.rope,'right'=1-p.left-p.rope)
}


