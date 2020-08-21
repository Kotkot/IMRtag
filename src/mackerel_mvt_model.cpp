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
  DATA_SCALAR(mean_diff_tag_area);
  DATA_VECTOR(is_from_west);
  DATA_VECTOR(thres_cov);                  		  //  response
  DATA_VECTOR(y);                  		  //  response
  DATA_INTEGER(Likconfig);

  // Parameters
  PARAMETER_MATRIX(beta);
  PARAMETER_VECTOR(log_sigma);

  vector<Type> sigma=exp(log_sigma);

  Type nll=0;
  // nll = -Nthres*log(Nthres);	// the constant for the prior not needed in TMB

  matrix<Type> mu(N,K);

  // Define mu
  for(int j=0;j<K;j++){
	// for (int i=0;i<N;i++){
	mu.col(j)=X*beta.col(j);
	// vector<Type> temp1=X.row(i);
	// vector<Type> temp2=beta.col(j);
	//  mu(i,j)=(temp1*temp2).sum();
	//  }
	 }

  // Define the gamma parameters (when needed)
  matrix<Type> shape(N,K);
  matrix<Type> scale(N,K);
  for(int j=0;j<K;j++){
	for (int n=0;n<N;n++){
		shape(n,j)= pow(mu(n,j),2)/pow(sigma(j),2);
		scale(n,j)= pow(sigma(j),2)/mu(n,j);
	}
  }

  // vector<Type> lp_e(Nthres+1);
  // vector<Type> lp_l(Nthres+1);
  matrix<Type> LL(N,Nthres);

    // for (int n=thresh_start(0);n<thresh_end(0);n++){	// for all observations below the first threshold
      // lp_e(0) += dnorm(y(n), mu(n,0), sigma(0), TRUE);  // before thr1
      // lp_l(0) += dnorm(y(n), mu(n,1), sigma(1), TRUE);  // after thr2
    // }  // Now define the mean travel distance for each mixture component

    // for (int t=0;t<Nthres;t++){
      // lp_e(t+1) = lp_e(t);
      // lp_l(t+1) = lp_l(t);
      // for (int n=thresh_start(t+1);n<thresh_end(t+1);n++){ 	  // for all observations between threshold [t,t+1[
        // lp_e(t + 1) += dnorm(y(n), mu(n,0), sigma(0), TRUE);  // before thr1
        // lp_l(t + 1) += dnorm(y(n), mu(n,1), sigma(1), TRUE);  // after thr2
      // }  // Now define the mean travel distance for each mixture component
    // }

    // for (int i=0;i<Nthres;i++){
	// nll -= lp_l(Nthres + 1) + lp_e(i) - lp_l(i);
	// }

	for (int thr=0;thr<Nthres;thr++) {
		for (int n=0;n<N;n++){
		// if (thres_cov(n) < (thresh(thr)+is_from_west(n)*mean_diff_tag_area)){
		if (thres_cov(n) < (thresh(thr))){
			if (Likconfig ==0){
				nll -= dnorm(y(n), mu(n,0), sigma(0), TRUE);			// for normal distribution on log scale
			  LL(n,thr)=dnorm(y(n), mu(n,0), sigma(0), TRUE);
			}
			if (Likconfig ==1){
				nll -= dgamma(y(n), shape(n,0), scale(n,0), TRUE);		// for gamma distribution on log scale
			  LL(n,thr)=dgamma(y(n), shape(n,0), scale(n,0), TRUE);
			}
		}
		// if (thres_cov(n) >= (thresh(thr)+is_from_west(n)*mean_diff_tag_area)){
		if (thres_cov(n) >= (thresh(thr))){
			if (Likconfig ==0){
				nll -= dnorm(y(n), mu(n,1), sigma(1), TRUE);
			  LL(n,thr)=dnorm(y(n), mu(n,1), sigma(1), TRUE);
			}
			if (Likconfig ==1){
				nll -= dgamma(y(n), shape(n,1), scale(n,1), TRUE);		// for gamma distribution on log scale
			  LL(n,thr)=dgamma(y(n), shape(n,1), scale(n,1), TRUE);
			}
		}
    }
  }

  // ============ Outputs =============

  // Model parameters
  REPORT(beta);
  REPORT(sigma);
  REPORT(LL);
  ADREPORT(mu);

  //--------------------------------------------------------------------

  return nll;

}
