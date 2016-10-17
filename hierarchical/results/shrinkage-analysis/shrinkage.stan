  data {
    int Nsamples; //number of results for each data set. Typically 100 (10 runs of 10-folds cv)
    int q; //number of data sets. 
    matrix[Nsamples,q] x; //difference of accuracy between the two classifier, on each fold of each data set.
    real std; //known std of the groups
     }



  parameters {
    real mu0; //uniform prior on all pars
    real<lower=0,upper=10> std0;
    vector[q] mu;               // mean  of each group
  }


  model {
  //to generate the data from multivariate normal with independent covariance
    matrix[q,q] Sigma; 
     for (j in 1:q){
       for (i in 1:q){
        if (j==i)
	       Sigma[j,i]<- pow(std,2);
        else
	       Sigma[j,i]<- 0;
       }
    }
    
    
    
  mu ~ normal(mu0, std0); 
        
    for (i in 1:Nsamples) {
      x[i] ~ multi_normal(mu,Sigma);
      }
      
}