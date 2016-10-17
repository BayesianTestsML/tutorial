  data {
    //parameters of the prior for sampling the sigma of each data set
    real unifLo; //bounds of the mean difference of accuracy within each data set
    real unifHi; 
    real mustdLo; //bounds of the mean difference of accuracy across the data sets
    real mustdHi; 

    int stud; // stud = 0/1 gaussian/student prior for mu[i] 

    int kruschke_prior; //Krushke prior for the degrees of freedom or Juarez and Steel prior (Model-Based Clustering of Non-Gaussian Panel Data Based on Skew-t Distributions)

    int<lower=2> Nsamples; //number of results for each data set. Typically 100 (10 runs of 10-folds cv)

    int<lower=1> q; //number of data sets. 

    matrix[q,Nsamples] x; //difference of accuracy between the two classifier, on each fold of each data set.
    

    real rho; //correlation

     }


  transformed data {

    //matrixes which do not depend on the data and which are helpful to compute the likelihood
    vector[Nsamples] H;
    vector[Nsamples] zero_mean_vec;
    matrix[Nsamples,Nsamples] invM;
    real detM;
    
    //parameter of the prior on the degrees of freedom
    real expLambda;
    expLambda <- 1/29.0 ;
  
    /*determinant of M, whose expression is analytically known
    and equal for each data set, as it only depends on rho and Nsamples */
    detM <- (1+(Nsamples-1)*rho)*(1-rho)^(Nsamples-1);

    
    //build H, invM and zero_mean_vec which will not change during sampling
    for (j in 1:Nsamples){
        zero_mean_vec[j]<-0;
        H[j]<-1;
       for (i in 1:Nsamples){
        //notice the mistake in our ML paper,Lemma 2.
        //Those are the entries of adj(M), not of inv(M).
        //inv(M)=1/det(M) adj(M)
        if (j==i)
	       invM[j,i]<- (1 + (Nsamples-2)*rho)*pow((1-rho),Nsamples-2);
        else
	       invM[j,i]<- -rho * pow((1-rho),Nsamples-2);
       }
    }
    //at this point invM contains the adjugate of M.
    //we  divide it by det(M) to obtain the inverse of M.
    invM <-invM/detM;
  }

  parameters {
    real mu0; //mean of the  hyperprior
    //uniform prior on std0 with  upper bound
    real<lower=mustdLo,upper=mustdHi> std0;
    
    
    vector[q] mu;               // mean  of each data set
    //vector<lower=unifLo,upper=unifHi>[q] sigma;   // variance of each data set

        //uniform prior on sigma with upper bound only
    vector<lower=unifLo,upper=unifHi>[q] sigma;   // std dev of each data set
    
    real<lower=0> nuMinusOne ; 
  }

 transformed parameters {
    real<lower=1> nu ;         // actually lower=1
    nu <- nuMinusOne + 1 ;
  }

  model {
    real partial_lik; //lik of the j-th group
    //real det_Sigma; //no longer used
    real log_det_Sigma;
    
    //mu0 does not appear here and it is given uniform priro by default
    //following Gelman we give uniform prior to the scale parameter
    //and we comment out its samplig statement.
    //std0  ~ uniform( mustdLo, mustdHi);
    
    if (kruschke_prior==1) {
       nuMinusOne ~ exponential( expLambda );
    } else {
       nuMinusOne ~ gamma ( 2, 0.1 );
    }
    
    //sample the data given the pars
    for (j in 1:q) {
      if(stud==0) {
        mu[j] ~ normal(mu0, std0); 
      } else {
        mu[j] ~ student_t(nu, mu0, std0); 
      }
      
      //uniform prior on the sigma[j] is implemented by removing the sampling statement.
      //sigma[j] ~ uniform( unifLo, unifHi);

        //naive implementation of det_sigma
        //det_Sigma <- pow(sigma[j],2*Nsamples) * detM;
       //as we need the log(det_sigma) working with logs is more stable
      log_det_Sigma <- 2*Nsamples *log(sigma[j]) + log(detM) ;
        
        
            
      //efficient computation of the likelihood based on the pre-computed inv(M)
       partial_lik <- -0.5 * 1/pow(sigma[j],2) * ( (x[j])' - H*mu[j] )' * invM * ((x[j])' - H*mu[j]);
        
        //partial_lik <- partial_lik - 0.5 * log (det_Sigma) - 0.5*Nsamples*log(6.283);
        partial_lik <- partial_lik - 0.5 * log_det_Sigma - 0.5*Nsamples*log(6.283);

        increment_log_prob(partial_lik);   
       }
  }