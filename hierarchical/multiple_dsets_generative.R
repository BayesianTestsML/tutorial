
multiple_dsets_generative <- function(delta_min,delta_max,how_many_dsets,sample_sizes,simulation_ID,actual_rho) {
  
  library(MASS)
  library(matrixStats)
  library(rstan)
  library(loo)
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  
  source('cv_nbc.R')
  source('hierarchical_test.R')
  delta_acc_sampling <- "gaussian"
  correlEstimation <- 1
  
  if  (( delta_acc_sampling!='cauchy') && (delta_acc_sampling!='gaussian') && (delta_acc_sampling!='mixture') ) { 
    stop('wrong delta_acc_sampling')
  }
  
  
  #check arguments
  std_upper_bound=1000
  n_runs <- 10
  n_folds <- 10
  #debug
  reps<-20
  
  
  file_str0 <- paste('results/deltaMin',delta_min,'_deltaMax_',delta_max,'_','sample_sizes',sample_sizes,'_',delta_acc_sampling,'_simulID_',simulation_ID,sep = "")
  filename <- paste(file_str0,'.csv',sep = "")
  log_filename <- paste('log_',file_str0,'.csv',sep = "")
  
  rope_min <- -0.01;
  rope_max <- 0.01;
  #set the seed for reproducible results
  set.seed(simulation_ID)
  
  
  if (any(delta_max == 0) || any(delta_max == 0)) "zero encountered"
  
  if (delta_max>0.5 || delta_min< -0.5){
    stop('absolute value of delta_acc should be smaller than 0.5');
  }
  
  if (delta_max < delta_min){
    stop('delta_max larger than delta_min');
  }
  
  #point estimator of the correlation
  rhoHat <- 1/n_folds;
  sign_rank_p_value <- matrix(0,reps);
  prob_left_hier_delta0 <- matrix(0,reps);
  prob_rope_hier_delta0 <- matrix(0,reps);
  prob_right_hier_delta0 <- matrix(0,reps);
  rmse_indep_bay_ttest <- matrix(0,reps);
  rmse_hier_bay_ttest <- matrix(0,reps);
  mae_hier_delta0 <- matrix(0,reps);
  rhoEstimated <- matrix(0,reps);
  stan_elapsed_time <- matrix(0,reps);
  mean.diff<- matrix(0,reps);
  var.diff<- matrix(0,reps);
  settings <- list() 
  results <- list() 
  counter <- 1;
  
  
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
    file_str <- paste('delta_min',delta_min,'_delta_max_',delta_max,'_dsets_',current_many_dsets,'_delta_acc_',delta_acc_sampling,'_samples_sizes_',sample_sizes,'simulationID_',simulation_ID,sep="");
    
    delta_std <- (delta_max-delta_min)/6;
    delta_mean <- (delta_max+delta_min)/2;
    
    initial.time <- proc.time()
    
    for (k in 1:reps) {
      
      if (k%%100==0) cat ('reps',k,'\n');
      
      #Gaussian or cauchy sampling
      if (delta_acc_sampling=='gaussian'){
        delta_acc_each_dset <- rnorm(current_many_dsets, mean = delta_mean, sd = delta_std);
      }
      else if (delta_acc_sampling=='cauchy'){
        delta_acc_each_dset <-  rt(current_many_dsets,3)*delta_std + delta_mean; 
      }
      else if (delta_acc_sampling=='mixture'){
        #random number from the uniform and sorted
        idx <- sort(runif(current_many_dsets))
        count1 <- max(1,sum (idx>0.5))
        count1 <- min(count1,current_many_dsets-1)
        count2 <- how_many_dsets - count1 
        std <- ((delta_max-delta_min)*1/4)/3 #three std fit in a quarter of delta_max - delta_min
        mean1 <- delta_min + (delta_max-delta_min)*1/4
        mean2 <- delta_min + (delta_max-delta_min)*3/4
        delta_acc_each_dset <- vector(length = current_many_dsets)
        delta_acc_each_dset[1:count1] <- rnorm(count1, mean = mean1, sd = std); 
        delta_acc_each_dset[(count1+1):(count1+count2)] <- rnorm(count2, mean = mean2, sd = std); 
      }
      
      indep_p_left_each_dset=matrix(1,current_many_dsets);
      indep_p_rope_each_dset=matrix(1,current_many_dsets);
      indep_p_right_each_dset=matrix(1,current_many_dsets);
      
      sample_diff_acc_a_b_each_dset <- matrix(0,current_many_dsets);
      
      
      Diff_ab <- matrix(0,n_runs*n_folds,current_many_dsets);
      
      #build the covariance matrix, equal for all dsets
      resultsEachDset <- n_folds * n_runs
      variance <- 0.02 * 0.98 /500;
      rho_var <- actual_rho * variance;
      covar_matrix <- matrix(rep.int (rho_var, resultsEachDset), nrow = resultsEachDset, ncol = resultsEachDset)
      for (i in 1:resultsEachDset){
        covar_matrix[i,i]<-variance
      }
      
      
      for (i in 1:current_many_dsets){
        current_delta_acc <- delta_acc_each_dset[i];
        diff_a_b <- rep.int (current_delta_acc, resultsEachDset) 
        diff_a_b <- diff_a_b + mvrnorm(n=1,mu=rep.int (0, resultsEachDset),Sigma=covar_matrix)
        sample_diff_acc_a_b_each_dset[i] <- mean(diff_a_b)
        
        #equivalent to p-value one-sided right tailed
        indep_p_each_dset <- ttest_Bayesian(diff_a_b,rhoHat,rope_min,rope_max)
        Diff_ab[,i]=diff_a_b
        indep_p_left_each_dset[i] <- indep_p_each_dset$left
        indep_p_rope_each_dset[i] <- indep_p_each_dset$rope
        indep_p_right_each_dset[i] <- indep_p_each_dset$right
      }
      
      #frequentist Wilcoxon signed rank
      sign_rank_p_value[k]  <-  wilcox.test(sample_diff_acc_a_b_each_dset,alternative = "two.sided",exact=FALSE)$p.value;
      # sign_rank_p_value[k]  <-  t.test(sample_diff_acc_a_b_each_dset)$p.value;
      mean.diff[k] <- mean(delta_acc_each_dset)
      var.diff[k] <- var(delta_acc_each_dset)
      
      stan.initial.time <- proc.time()
      stan_prob <- hierarchical.test(t(Diff_ab),rhoHat,rope_min,rope_max,file_str,std_upper_bound,correlEstimation)
      stan_elapsed_time[k]<-(proc.time()[3]-stan.initial.time[3]);
      show(stan_elapsed_time[k])
      prob_left_hier_delta0[k]  <- stan_prob$delta0$left
      prob_rope_hier_delta0[k]  <- stan_prob$delta0$rope
      prob_right_hier_delta0[k] <- stan_prob$delta0$right
      rhoEstimated[k] <- mean(stan_prob$rho)
      
      #mae on the delta0. Actual delta0 is halfway delta_min between delta_max
      #the mae on delta0 is computed on the probability. It could be computed
      #on the delta0 vs estimated delta0, for consistency with what we do on dsets.
      actualdelta0 <- mean(c(delta_min,delta_max))
      if ( actualdelta0 < rope_min ) {
        mae_hier_delta0[k] <- 1-prob_left_hier_delta0[k]
      } else if ( actualdelta0 > rope_max ) {
        mae_hier_delta0[k] <- 1-prob_right_hier_delta0[k]
      } else {
        mae_hier_delta0[k] <- 1-prob_rope_hier_delta0[k]
      }
      
      #means a posteriori equals the maxLik mean
      maxLikMeans <- colMeans(Diff_ab)
      rmse_indep_bay_ttest[k] <- mean( (delta_acc_each_dset- maxLikMeans)^2 )
      rmse_hier_bay_ttest[k] <- mean( (delta_acc_each_dset- stan_prob$delta_each_dset)^2 )
    }
    
    elapsed_time <- (proc.time()[3]-initial.time[3]);
    cat('elapsed time: ', elapsed_time)
    
    results[[j]] <- list('how_many_dsets'=current_many_dsets,
                         'sample_size'=current_dset_size,
                         'delta_min'=delta_min,
                         'delta_max'=delta_max,
                         'sign_rank_p_value'=sign_rank_p_value,
                         'prob_left_hier_delta0'=prob_left_hier_delta0,
                         'prob_rope_hier_delta0'=prob_rope_hier_delta0,
                         'prob_right_hier_delta0'=prob_right_hier_delta0,
                         'rmse_hier_each_dset'=rmse_hier_bay_ttest,
                         'rmse_indep_each_dset'=rmse_indep_bay_ttest,
                         'mae_hier_delta0'=mae_hier_delta0,
                         'stan_elapsed_time'=stan_elapsed_time,
                         'rhos'=mean(rhoEstimated)
    )
  }
  
  save_results <- function(results,filename) {
    
    mystring <- paste('delta_acc_sampling,correlEstimation,delta_acc_min,delta_acc_max,num_experiments,sample_size,',
                      'how_many_dsets,sample_size,sign_rank_rejects_null,hier_left_95,hier_rope_95,hier_right95,',
                      'hier_left_90,hier_rope_90,hier_right90,median_hier_p_left,median_hier_p_rope,',
                      'median_hier_p_right,rmse_ttest_hier,rmse_ttest_indep,maedelta0,EstimatedRho,Nfolds,StanElapsedTime',sep="");
    write(mystring, filename, append = FALSE)
    
    for (ii in 1:length(results)){
      num_experiments <- length(results[[ii]]$sign_rank_p_value);
      h1_sign_rank <- mean(results[[ii]]$sign_rank_p_value<.05);
      hier_left_95 <- mean(results[[ii]]$prob_left_hier_delta0>.95);
      hier_rope_95 <- mean(results[[ii]]$prob_rope_hier_delta0>.95);
      hier_right_95 <- mean(results[[ii]]$prob_right_hier_delta0>.95);
      hier_left_90 <- mean(results[[ii]]$prob_left_hier_delta0>.90);
      hier_rope_90 <- mean(results[[ii]]$prob_rope_hier_delta0>.90);
      hier_right_90 <- mean(results[[ii]]$prob_right_hier_delta0>.90);
      
      tmp_vector <- c(h1_sign_rank, hier_left_95,hier_rope_95, hier_right_95,hier_left_90,hier_rope_90,hier_right_90,
                      median(results[[ii]]$prob_left_hier_delta0),median(results[[ii]]$prob_rope_hier_delta0),median(results[[ii]]$prob_right_hier_delta0),
                      mean(results[[ii]]$rmse_hier_each_dset),mean(results[[ii]]$rmse_indep_each_dset),
                      mean(results[[ii]]$mae_hier_delta0), mean(results[[ii]]$rho), n_folds, mean(results[[ii]]$stan_elapsed_time));
      mystring <- paste(delta_acc_sampling,',',correlEstimation,',',results[[ii]]$delta_min,',',
                        results[[ii]]$delta_max,',',num_experiments,',',sample_sizes,',',results[[ii]]$how_many_dsets,
                        ',',results[[ii]]$sample_size,',',paste(tmp_vector,collapse=","));
      write(mystring, filename, append = TRUE)
    }
  }
  
  #save csv file
  save_results(results, filename);
  
  #save the Rdata file
  rdata_filename <- paste(file_str0,'.Rdata',sep="");
  save(results, file = rdata_filename)
  
}

ttest_Bayesian <- function(diff_a_b,rhoHat,rope_min,rope_max) {
  delta <- mean(diff_a_b)
  n <- length(diff_a_b)
  df <- n-1
  stdX <- sd(diff_a_b)
  sp <- sd(diff_a_b)*sqrt(1/n + rhoHat/(1-rhoHat))
  p.left <- pt((rope_min - delta)/sp, df)
  p.rope <- pt((rope_max - delta)/sp, df)-p.left
  indep_p_each_dset <- list('left'=p.left,'rope'=p.rope,'right'=1-p.left-p.rope)
}


