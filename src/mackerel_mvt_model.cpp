////////////////////////////////////////////////////////////////
// Mixture model to analyse the within-year mackerela travel rate
// Author: Kotaro Ono
// Version: 1
// Detail: TMB code
// Applies a change point model for fast vs slower moving mackerel 
// in development (to replace the Stan model that runs quite slowly...)


#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;
  using namespace Eigen;

  // Data for the maturity ogive at length part
  DATA_INTEGER(K);                      //  number of mixture components
  DATA_INTEGER(N);                      //  number of data points
  DATA_INTEGER(Nx);                     //  number of columns in the design matrix
  DATA_MATRIX(X);                  		  //  the design matrix
  DATA_INTEGER(Nthres);
  DATA_VECTOR(thres);
  DATA_SCALAR(mean_diff_tag_area);
  DATA_IVECTOR(is_from_west);
  DATA_VECTOR(Y);                  		  //  response

  // Parameters
  PARAMETER_MATRIX(beta);
  PARAMETER_VECTOR(sigma);

  Type nll=0;
  for (int i=0;i<Nthres;i++){
  nll -= log(Nthres);
  }
  
  matrix<Type> mu(N,K);
  
  // Define mu
  for(int j=0;j<K;j++){
	  mu(,j)=X*beta.col(j);
  }
 
  for (int thr=0;thr<Nthres;thr++) {
    for (int n=0;n<N;n++)
      nll -= dnorm(y[n] | y[n] < (thresh[thr]+is_from_west[n]*mean_diff_tag_area) ? mu[n,1] : mu[n,2], sigma); // to change to TMB code
  }

 
  
  // ============ Outputs =============

  // Model parameters
  REPORT(beta);
  REPORT(sigma);

  //--------------------------------------------------------------------

  return nll;

}
