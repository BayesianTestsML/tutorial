#diff_a_b is a vector of differences between the two classifiers, on each fold of cross-validation.
#If you have done 10 runs of 10-folds cross-validation, you have 100 results for each classifier.
#You should have run cross-validation on the same folds for the two classifiers.
#Then diff_a_b is the difference fold-by-fold.

#rho is the correlation of the cross-validation results: 1/(number of folds)
#rope_min and rope_max are the lower and the upper bound of the rope
 correlatedBayesianTtest <- function(diff_a_b,rho,rope_min,rope_max){
   if (rope_max < rope_min){
     stop("rope_max should be larger than rope_min")
   }
     
  delta <- mean(diff_a_b)
  n <- length(diff_a_b)
  df <- n-1
  stdX <- sd(diff_a_b)
  sp <- sd(diff_a_b)*sqrt(1/n + rho/(1-rho))
  p.left <- pt((rope_min - delta)/sp, df)
  p.rope <- pt((rope_max - delta)/sp, df)-p.left
  results <- list('left'=p.left,'rope'=p.rope,'right'=1-p.left-p.rope)
  return (results)
}


