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
  DATA_MATRIX(X);                  		  //  the design matrix
  DATA_INTEGER(Nthres);
  DATA_VECTOR(thresh);
  DATA_IVECTOR(thresh_start);
  DATA_IVECTOR(thresh_end);
  DATA_SCALAR(mean_diff_tag_area);
  DATA_VECTOR(is_from_west);
  DATA_VECTOR(y);                  		  //  response

  // Parameters
  PARAMETER_MATRIX(beta);
  PARAMETER_VECTOR(log_sigma);

  vector<Type> sigma=exp(log_sigma);

  Type nll=0;
  // nll = -Nthres*log(Nthres);	// the constant for the prior not needed in TMB

  matrix<Type> mu(N,K);

  // Define mu
  for(int j=0;j<K;j++){
	mu.col(j)=X*beta.col(j);
  }

  vector<Type> lp_e(Nthres+1);
  vector<Type> lp_l(Nthres+1);

    for (int n=thresh_start(0);n<thresh_end(0);n++){	// for all observations below the first threshold
      lp_e(0) += dnorm(y(n), mu(n,0), sigma(0), TRUE);  // before thr1
      lp_l(0) += dnorm(y(n), mu(n,1), sigma(1), TRUE);  // after thr2
    }  // Now define the mean travel distance for each mixture component

    for (int t=0;t<Nthres;t++){
      lp_e(t+1) = lp_e(t);
      lp_l(t+1) = lp_l(t);
      for (int n=thresh_start(t+1);n<thresh_end(t+1);n++){ 	  // for all observations between threshold [t,t+1[
        lp_e(t + 1) += dnorm(y(n), mu(n,0), sigma(0), TRUE);  // before thr1
        lp_l(t + 1) += dnorm(y(n), mu(n,1), sigma(1), TRUE);  // after thr2
      }  // Now define the mean travel distance for each mixture component
    }
  
    for (int i=0;i<Nthres;i++){
	nll -= lp_l(Nthres + 1) + lp_e(i) - lp_l(i);
	}
	
  // for (int thr=0;thr<Nthres;thr++) {
    // for (int n=0;n<N;n++){
		// if (y(n) < (thresh(thr)+is_from_west(n)*mean_diff_tag_area)){
			// nll -= dnorm(y(n), mu(n,0), sigma(0), TRUE);
		// }
		// if (y(n) >= (thresh(thr)+is_from_west(n)*mean_diff_tag_area)){
			// nll -= dnorm(y(n), mu(n,1), sigma(1), TRUE);
		// }
    // }
  // }

  // ============ Outputs =============

  // Model parameters
  REPORT(beta);
  REPORT(sigma);
  ADREPORT(mu);

  //--------------------------------------------------------------------

  return nll;

}
