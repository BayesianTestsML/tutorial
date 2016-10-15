cv_competing_nbc <- function(nruns,nfolds,n,theta_diff,theta_star,seed=0) {
  # [results]=cv_competing_nbc(nruns,nfolds,n,theta_diff,theta_star)
  #   implements cross-validation, given the number or runs, folds, 
  #   data set size (n) and theta. Compares 2 nbc which uses two different feature for the prediction .
  #   Theta determines the conditional prob of the feature while
  #   generating the data. As a consequence the expected accuracy of naive Bayes is
  #   theta1 and theta2 respectively for the two features. 
  #   theta1=1, namely the first feature is a copty of the class and yields a
  #   perfect classifier.
  #   theta2=1-theta_diff; theta2 is typically high and allows to correclty
  #   learnt the bias.
  #   Both classifiers are random guessers and detected significance are type I
  #   errors.
  #   cv_results contains:
  #     -nbc1 accuracy
  #     -nbc2 accuracy
  #     -delta = nbc1 accuracy - nbc2 accuracy
  #     -var=variance  (nbc1 accuracy - nbc2 accuracy)
  
  
  c_val <- function(nfolds) {
    cv <-list()
    cv$nbc1 <- matrix(1,nfolds);
    cv$nbc2 <- cv$nbc1;
    cv$delta <- cv$nbc1;
    cv$var <- cv$nbc1;
    #column of the feature in the data set for nbc1 and nbc2
    feat_idx_nbc1 <- 2;
    feat_idx_nbc2 <- 3;
    boring_idx <- rep(1:nfolds,ceiling(dim(data)[1]/nfolds));
    boring_idx <- boring_idx[1:dim(data)[1]];
    permuted_idx <- sample(dim(data)[1]);
    permuted_sorted_data <- data[permuted_idx,];
    permuted_sorted_data <- permuted_sorted_data[order(permuted_sorted_data[,1]),];
       
    for (j in 1:nfolds) {
      train <- permuted_sorted_data[which(boring_idx!=j),];
      test <- permuted_sorted_data[which(boring_idx==j),];
      #with the perfect feature
      nbc1 <- learn_nbc(train[,c(1,2)]);
      #with the stochastic feature
      nbc2 <- learn_nbc(train[,c(1,3)]);
      preds_nbc1 <- test_nbc(nbc1,test,feat_idx_nbc1);
      preds_nbc2 <- test_nbc(nbc2,test,feat_idx_nbc2);
      #acc_nbc
      nbc1_correct <- preds_nbc1==test[,1];
      cv$nbc1[j] <- mean(nbc1_correct);
      #acc_zero
      nbc2_correct <- preds_nbc2==test[,1];
      cv$nbc2[j] <- mean(nbc2_correct);
      cv$delta[j] <- mean(nbc1_correct-nbc2_correct);
      #var_delta_acc
      cv$var[j] <- var(nbc1_correct-nbc2_correct);
    }
    cv
  }
  
  
  theta1 <- theta_star;
  theta2 <- theta1-theta_diff;
  data <- generate_data_with_two_feats(n,theta1,theta2);
  
  
  for (i in 1:nruns) {
    if (i==1) {
      cv_results <- c_val(nfolds) 
    }
    else {
      tmp <- c_val(nfolds);
      cv_results$nbc1 <- c(cv_results$nbc1, tmp$nbc1);
      cv_results$nbc2 <- c(cv_results$nbc2, tmp$nbc2);
      cv_results$delta <- c(cv_results$delta, tmp$delta);
      cv_results$var <- c(cv_results$var, tmp$var);
    }
  }
  
 return (cv_results)
  
}

generate_data_with_two_feats <- function(n,theta1,theta2) {
  ##function data=generate_data(n,theta1,theta2)
  #generate n artificial instances with binary class and two binary features.
  #The bias of feature1 and feature2 is theta1 and theta2, namely P(f1|c1)=theta1;
  #P(~f1|~c1)=theta1; P(f2|c1)=theta2, P(~f2|~c1)=theta2.
  
  generate_feature <- function(theta) {
    #feature is obtained by flipping randomly (1-theta)# of class labels
    
    feature <- class-1;
    
    tmp <- runif(n)>theta;
    rand_idx <- tmp>theta;
    feature[tmp] <- 1-feature[tmp];
    
    feature+1
  }
  
  class <- 1*(runif(n)>0.5)+1;
  
  feature1 <- generate_feature(theta1);
  feature2 <- generate_feature(theta2);
  
  data <- cbind(class,feature1,feature2);
  
  data
  
}

learn_nbc <- function(train) {
  ##nbc=learn_nbc(train)
  #learns naive Bayes (BDeu prior) from train data, assuming the class to be in first
  #column and the feature to be in second column.
  class_idx <- 1;
  feat_idx <- 2;
  
  marg <- (sum(train[,class_idx]==1)+.5)/(dim(train)[1]+1);
  marg <- c(marg,1-marg);
  nbc <- list('marg'=marg);
  
  cond <- matrix(0,2,2)
  cond[1,1] <- (sum(train[,class_idx]==1 & train[,feat_idx]==1) + .25) / (sum(train[,feat_idx]==1)+.5);
  cond[2,1] <- 1-cond[1,1];
  cond[1,2] <- (sum(train[,class_idx]==2 & train[,feat_idx]==1) + .25) / ( sum(train[,feat_idx]==2) +.5);
  cond[2,2] <- 1-cond[1,2];
  nbc$cond <- cond;
  
  nbc
}


test_nbc <- function(nbc,test,feat_idx) {
  ##preds_nbc=test_nbc(nbc,test,feat_idx)
  #Returns the most probable class predicted by naiveBayes for each instance of the test set.
  #The test set is assumed to contain a binary class and a
  #a binary feature.
  
  nargin <- length(as.list(match.call()))-1
  if (nargin<3) feat_idx <- 2;
  
  preds_nbc <- matrix(1,dim(test)[1]);  
  for (ii in 1:dim(test)[1]) {
    current_feat <- test[ii,feat_idx];
    prob <- c(nbc$marg[1]*nbc$cond[current_feat,1], nbc$marg[2]*nbc$cond[current_feat,2]);
    if (prob[2]>prob[1]) { 
      preds_nbc[ii] <- 2;
    }
  }
  
  preds_nbc
}
