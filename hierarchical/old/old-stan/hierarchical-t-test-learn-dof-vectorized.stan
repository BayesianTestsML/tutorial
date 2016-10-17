  data {
    //parameters of the prior for sampling the sigma of each data set
    real unifLo; //bounds of the mean difference of accuracy within each data set
    real unifHi; 
    real mustdLo; //bounds of the mean difference of accuracy across the data sets
    real mustdHi; 
    int samplingType; //0 for student with Krushke prior for the degrees of freedom; 1 for student with Juarez and Steel prior (Model-Based Clustering of Non-Gaussian Panel Data Based on Skew-t Distributions) on the degrees of freedom; 2 for Gaussian

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
    
    //parameter of the prior on the degrees of freedom for Kruschke prior
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
    
    //uniform prior on std0 over the given interval
    real<lower=mustdLo,upper=mustdHi> std0;
    
    
    vector[q] mu;               // mean  of each data set

    //uniform prior on sigma of each data set
    vector<lower=unifLo,upper=unifHi>[q] sigma; 
    
    real<lower=0> nuMinusOne ; 
  }

 transformed parameters {
    real<lower=1> nu ;         
    matrix[q,Nsamples] diff; //difference between the data and the means
    vector[q] diagQuad;
    vector[q] oneOverSigma2; //vector: 1 over the variance of each data set
    
    //actual degrees of freedom
    nu <- nuMinusOne + 1 ;
    
    //the data (x) minus a matrix done as follows:
    //the mu vector (of lenght q) repeated side by side in  Nsamples columns
    diff <- x - rep_matrix(mu,Nsamples); 
    diagQuad <- diagonal (quad_form (invM,diff'));
    oneOverSigma2 <- rep_vector(1, q) ./ sigma;
    oneOverSigma2 <- oneOverSigma2 ./ sigma;
    
  }

  model {
   vector[q] log_det_Sigma;
   vector[q] partial_lik;
    
    //mu0 does not appear here and it is given uniform prior by default
    
    //std is uniformly dranw over its domain
    //std0  ~ uniform( mustdLo, mustdHi);
    
    //kruschke sampling
    if (samplingType==0) {
       nuMinusOne ~ exponential( expLambda );
       mu ~ student_t(nu, mu0, std0);
    } 
    
    //steel sampling
    else if (samplingType==1) {
       nuMinusOne ~ gamma ( 2, 0.1 );
       mu ~ student_t(nu, mu0, std0);
    }
    //gaussian sampling
    else {
              mu ~ normal(mu0, std0); 
    }

    //uniform prior on the sigma[j] is implemented by removing the sampling statement.
    //sigma[j] ~ uniform( unifLo, unifHi);
  
    log_det_Sigma <- 2*Nsamples*log(sigma) + log(detM) ;
    partial_lik <- -0.5 * log_det_Sigma - 0.5*Nsamples*log(6.283);  
    partial_lik <- partial_lik - 0.5 * oneOverSigma2 .* diagQuad;
    
    increment_log_prob(sum(partial_lik));   
 }