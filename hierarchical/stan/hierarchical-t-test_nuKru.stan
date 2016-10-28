/*Hierarchical Bayesian model for the analysis of competing cross-validated classifiers on multiple data sets.
Adopts Kruschke prior as for the degrees of freedom, Gamma(1,0.0345)
*/

  data {
    
    real deltaLow;
    real deltaHi;
    
    //bounds of the sigma of the higher-level distribution
    real std0Low; 
    real std0Hi; 

    //bounds on the domain of the sigma of each data set
    real stdLow; 
    real stdHi; 
    
        
    //number of results for each data set. Typically 100 (10 runs of 10-folds cv)
    int<lower=2> Nsamples; 

    //number of data sets. 
    int<lower=1> q; 

    //difference of accuracy between the two classifier, on each fold of each data set.
    matrix[q,Nsamples] x;
    
    //correlation (1/(number of folds))
    real rho; 
     }


  transformed data {

    //vector of 1s appearing in the likelihood 
    vector[Nsamples] H;
    
    //vector of 0s: the mean of the mvn noise 
    vector[Nsamples] zeroMeanVec;
    
    /* M is the correlation matrix of the mvn noise.
    invM is its inverse, detM its determinant */
    matrix[Nsamples,Nsamples] invM;
    real detM;
    
    //The determinant of M is analytically known
    detM <- (1+(Nsamples-1)*rho)*(1-rho)^(Nsamples-1);

    //build H and invM. They do not depend on the data.
    for (j in 1:Nsamples){
      zeroMeanVec[j]<-0;
      H[j]<-1;
      for (i in 1:Nsamples){
        if (j==i)
          invM[j,i]<- (1 + (Nsamples-2)*rho)*pow((1-rho),Nsamples-2);
        else
          invM[j,i]<- -rho * pow((1-rho),Nsamples-2);
       }
    }
    /*at this point invM contains the adjugate of M.
    we  divide it by det(M) to obtain the inverse of M.*/
    invM <-invM/detM;
  }

  parameters {
    //mean of the  hyperprior from which we sample the delta_i
    real<lower=deltaLow,upper=deltaHi> delta0; 
    
    //std of the hyperprior from which we sample the delta_i
    real<lower=std0Low,upper=std0Hi> std0;
    
    //delta_i of each data set: vector of lenght q.
    vector[q] delta;               

    //sigma of each data set: : vector of lenght q.
    vector<lower=stdLow,upper=stdHi>[q] sigma; 
    
    /* the domain of (nu - 1) starts from 0
    and can be given a gamma prior*/
    real<lower=0> nuMinusOne; 
    
    //parameters of the Gamma prior on nuMinusOne
    real<lower=1,upper=2> gammaAlpha;
    real<lower=0.01, upper=0.1> gammaBeta;
    
  }

 transformed parameters {
    //degrees of freedom
    real<lower=1> nu ;
    
    /*difference between the data (x matrix) and 
    the vector of the q means.*/
    matrix[q,Nsamples] diff; 
    
    vector[q] diagQuad;
    
    /*vector of length q: 
    1 over the variance of each data set*/
    vector[q] oneOverSigma2; 
    
    vector[q] logDetSigma;
    
    vector[q] logLik;
   
    //degrees of freedom
    nu <- nuMinusOne + 1 ;
    
    //1 over the variance of each data set
    oneOverSigma2 <- rep_vector(1, q) ./ sigma;
    oneOverSigma2 <- oneOverSigma2 ./ sigma;

    /*the data (x) minus a matrix done as follows:
    the delta vector (of lenght q) pasted side by side Nsamples times*/
    diff <- x - rep_matrix(delta,Nsamples); 
    
    //efficient matrix computation of the likelihood.
    diagQuad <- diagonal (quad_form (invM,diff'));
    logDetSigma <- 2*Nsamples*log(sigma) + log(detM) ;
    logLik <- -0.5 * logDetSigma - 0.5*Nsamples*log(6.283);  
    logLik <- logLik - 0.5 * oneOverSigma2 .* diagQuad;
    
  }

  model {
    /*mu0 and std0 are not explicitly sampled here.
    Stan automatically samples them: mu0 as uniform and std0 as
    uniform over its domain (std0Low,std0Hi).*/

    //sampling the degrees of freedom according to Kruschke prior
    nuMinusOne ~ gamma ( 1, 0.0345);
    
    
    //vectorial sampling of the delta_i of each data set
    delta ~ student_t(nu, delta0, std0);
    
    //logLik is computed in the previous block 
    increment_log_prob(sum(logLik));   
 }
 
 
