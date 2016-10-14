multiple_dsets_rope <- function(delta0,std0,how_many_dsets,sample_sizes,reps,simulation_ID=1, delta_acc_sampling='cauchy',
                                flag500=FALSE) {
  
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
  
  
  
  #check arguments
  std_upper_bound=1000
  
  
  showPost <- 0
  trackMeans<-0
  
  # this control the variance of the simulation.
  # one classifier has accuracy theta_star, the other theta_star=delta_acc.
  # with theta_star=1, we have zero variance.
  # witht theta_star=0.5, we have maximum variance
  theta_star <- 0.9
  n_folds <- 10
  n_runs <- 10
  
  
  if (trackMeans==1){
    reps <- 1 
  }
  
  file_str0 <- paste('resultsDelta0',delta0,'Std0',std0,'_','sample_sizes',sample_sizes,'_',delta_acc_sampling,'_simulID_',simulation_ID,".dat",sep = "")
  filename <- paste(file_str0,'.csv',sep = "")
  log_filename <- paste('log_',file_str0,'.csv',sep = "")
  
  rope_min <- -0.01;
  rope_max <- 0.01;
  #set the seed for reproducible results
  set.seed(simulation_ID)
  
  
  
  
  #point estimator of the correlation
  rho <- 1/n_folds;
  sign_rank_p_value <- matrix(0,reps)
  prob_left_hier_delta0 <- matrix(0,reps)
  prob_rope_hier_delta0 <- matrix(0,reps)
  prob_right_hier_delta0 <- matrix(0,reps)
  rmse_indep_bay_ttest <- matrix(0,reps)
  rmse_hier_bay_ttest <- matrix(0,reps)
  maeMle_Mle500 <- matrix(0,reps)
  maeHier_Mle500 <-   matrix(0,reps)
  maeDelta_Mle500 <- matrix(0,reps)
  rhoEstimated <- matrix(0,reps)
  stan_elapsed_time <- matrix(0,reps)
  mean.diff<- matrix(0,reps)
  var.diff<- matrix(0,reps)
  settings <- list() 
  results <- list() 
  counter <- 1
  
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
    
    initial.time <- proc.time()
    
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
      if (flag500==TRUE){
        delta_acc_cv500 <- vector(length = current_many_dsets)
      }
      
      indep_p_left_each_dset=matrix(1,current_many_dsets);
      indep_p_rope_each_dset=matrix(1,current_many_dsets);
      indep_p_right_each_dset=matrix(1,current_many_dsets);
      
      sample_diff_acc_a_b_each_dset <- matrix(0,current_many_dsets);
      
      
      Diff_ab <- matrix(0,n_runs*n_folds,current_many_dsets);
      
      for (i in 1:current_many_dsets){
        
        current_delta_acc <- delta_acc_each_dset[i];
        cv_results <- list('delta'=rep(0,n_runs*n_folds))
        if (flag500==FALSE){
          while (var(cv_results$delta)==0 && mean(cv_results$delta)==0) { 
            cv_results <- cv_competing_nbc(n_runs,n_folds,current_dset_size,current_delta_acc,theta_star);
            diff_a_b <- cv_results$delta;
          }
        }
        if (flag500==TRUE){
          for (cvIter in 1:500){
            cv_results <- cv_competing_nbc(n_runs,n_folds,current_dset_size,current_delta_acc,theta_star,seed=cvIter)
            delta_acc_cv500[cvIter] <- mean(cv_results$delta)
#             cat('iter:',cvIter, '\n')
#             cat(mean(cv_results$delta), sd(cv_results$delta),'\n')
          }
          cv_results <- list('delta'=rep(0,n_runs*n_folds));
          while (var(cv_results$delta)==0 && mean(cv_results$delta)==0) { 
            cv_results <- cv_competing_nbc(n_runs,n_folds,current_dset_size,current_delta_acc,theta_star)
            diff_a_b <- cv_results$delta
          }
        }
        
        sample_diff_acc_a_b_each_dset[i] <- mean(diff_a_b);
        
        
        #equivalent to p-value one-sided right tailed
        indep_p_each_dset <- ttest_Bayesian(diff_a_b,rho,rope_min,rope_max);
        Diff_ab[,i]=diff_a_b;
        indep_p_left_each_dset[i] <- indep_p_each_dset$left;
        indep_p_rope_each_dset[i] <- indep_p_each_dset$rope;
        indep_p_right_each_dset[i] <- indep_p_each_dset$right;
        
      }
      
      #frequentist Wilcoxon signed rank
      sign_rank_p_value[k]  <-  wilcox.test(sample_diff_acc_a_b_each_dset,alternative = "two.sided",exact=FALSE)$p.value;
      # sign_rank_p_value[k]  <-  t.test(sample_diff_acc_a_b_each_dset)$p.value;
      mean.diff[k] <- mean(delta_acc_each_dset)
      var.diff[k] <- var(delta_acc_each_dset)
      
      stan.initial.time <- proc.time()
      stan_prob <- hierarchical.test(t(Diff_ab),rho,rope_min,rope_max,file_str,std_upper_bound,'student')
      stan_elapsed_time[k]<-(proc.time()[3]-stan.initial.time[3]);
      show(stan_elapsed_time[k])
      prob_left_hier_delta0[k]  <- stan_prob$delta0$left
      prob_rope_hier_delta0[k]  <- stan_prob$delta0$rope
      prob_right_hier_delta0[k] <- stan_prob$delta0$right
      
      #logLoss, pActual, etc commented out: it refers to the inference on Delta0      
      #       for (aa in 1:current_many_dsets){
      #         if (delta_acc_each_dset[aa]<rope_min) {
      #           logLossIndep[k] <- logLossIndep[k] - log(indep_p_left_each_dset[aa])
      #           logLossHier[k] <- logLossHier[k] - log(stan_prob$delta$left[aa])
      #           pActualHier[k] <-pActualHier[k] + stan_prob$delta$left[aa]
      #           pActualIndep[k] <-pActualIndep[k] + indep_p_left_each_dset[aa]
      #         } else if (delta_acc_each_dset[aa]>rope_max) {
      #           logLossIndep[k] <- logLossIndep[k] - log(indep_p_right_each_dset[aa])
      #           logLossHier[k] <- logLossHier[k] - log(stan_prob$delta$right[aa])
      #           pActualHier[k] <-pActualHier[k] + stan_prob$delta$right[aa]
      #           pActualIndep[k] <-pActualIndep[k] + indep_p_right_each_dset[aa]
      #         } else {
      #           logLossIndep[k] <- logLossIndep[k] - log(indep_p_rope_each_dset[aa])
      #           logLossHier[k] <- logLossHier[k] - log(stan_prob$delta$rope[aa])
      #           pActualHier[k] <-pActualHier[k] + stan_prob$delta$rope[aa] 
      #           pActualIndep[k] <-pActualIndep[k] + indep_p_rope_each_dset[aa]
      #         }
      #       }
      #       logLossIndep[k] <- logLossIndep[k]/current_many_dsets
      #       logLossHier[k] <- logLossHier[k]/current_many_dsets
      #       pActualIndep[k] <-pActualIndep[k]/current_many_dsets
      #       pActualHier[k] <-pActualHier[k]/current_many_dsets
      
      #means a posteriori equals the maxLik mean
      maxLikMeans <- colMeans(Diff_ab)
      rmse_indep_bay_ttest[k] <- sqrt(mean( (delta_acc_each_dset- maxLikMeans)^2 ))
      rmse_hier_bay_ttest[k] <- sqrt(mean( (delta_acc_each_dset- stan_prob$delta_each_dset)^2 ))
      
      maeMle_Mle500[k] <- (mean( abs(delta_acc_cv500- maxLikMeans) ))
      maeHier_Mle500[k] <- (mean( abs(delta_acc_cv500- stan_prob$delta_each_dset) ))
      maeDelta_Mle500[k] <- (mean( abs(delta_acc_cv500 - delta_acc_each_dset ) ))
      
    }#close the loop on the repetitions
    
    
    elapsed_time <- (proc.time()[3]-initial.time[3]);
    cat('elapsed time: ', elapsed_time)
    
    results[[j]] <- list('how_many_dsets'=current_many_dsets,
                         'sample_size'=current_dset_size,
                         'delta0'=delta0,
                         'std0'=std0,
                         'sign_rank_p_value'=sign_rank_p_value,
                         'prob_left_hier_delta0'=prob_left_hier_delta0,
                         'prob_rope_hier_delta0'=prob_rope_hier_delta0,
                         'prob_right_hier_delta0'=prob_right_hier_delta0,
                         'rmse_hier_each_dset'=rmse_hier_bay_ttest,
                         'rmse_indep_each_dset'=rmse_indep_bay_ttest,
                         'stan_elapsed_time'=stan_elapsed_time)  
  }
  
  save_results <- function(results,filename) {
    
    mystring <- paste('delta_acc_sampling,delta0,std0,num_experiments,',
                      'how_many_dsets,sample_size,signRankPower,signRankPValue,hier_left_95,hier_rope_95,hier_right95,',
                      'hier_left_90,hier_rope_90,hier_right90,median_hier_p_left,median_hier_p_rope,',
                      'median_hier_p_right,rmse_hier,rmse_mle,StanElapsedTime',sep="");
    
    #THESE FIELDS ARE NO LONGER TRACKED: logLossHier,logLossIndep,pActualHier,pActualIndep
    write(mystring, filename, append = FALSE)
    
    for (ii in 1:length(results)){
      num_experiments <- length(results[[ii]]$sign_rank_p_value);
      h1_sign_rank <- mean(results[[ii]]$sign_rank_p_value<.05);
      h1_sign_rankPvalue <- median(results[[ii]]$sign_rank_p_value);
      hier_left_95 <- mean(results[[ii]]$prob_left_hier_delta0>.95);
      hier_rope_95 <- mean(results[[ii]]$prob_rope_hier_delta0>.95);
      hier_right_95 <- mean(results[[ii]]$prob_right_hier_delta0>.95);
      hier_left_90 <- mean(results[[ii]]$prob_left_hier_delta0>.90);
      hier_rope_90 <- mean(results[[ii]]$prob_rope_hier_delta0>.90);
      hier_right_90 <- mean(results[[ii]]$prob_right_hier_delta0>.90);
      
      tmp_vector <- c(h1_sign_rank, h1_sign_rankPvalue, hier_left_95,hier_rope_95, hier_right_95,hier_left_90,hier_rope_90,hier_right_90,
                      median(results[[ii]]$prob_left_hier_delta0),median(results[[ii]]$prob_rope_hier_delta0),median(results[[ii]]$prob_right_hier_delta0),
                      mean(results[[ii]]$rmse_hier_each_dset),mean(results[[ii]]$rmse_indep_each_dset),
                      mean(results[[ii]]$stan_elapsed_time));
      mystring <- paste(delta_acc_sampling,',',results[[ii]]$delta0,',',
                        results[[ii]]$std0,',',num_experiments,',',results[[ii]]$how_many_dsets,
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


