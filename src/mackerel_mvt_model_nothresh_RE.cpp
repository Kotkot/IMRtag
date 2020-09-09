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
  DATA_INTEGER(N);                      //  number of data points
  DATA_MATRIX(X);                  		  //  the design matrix
  DATA_INTEGER(Nthres);
  DATA_VECTOR(y);                  		  //  response
  DATA_MATRIX(X_pred);                  		  //  the design matrix
  DATA_INTEGER(N_pred);                  		  //  the design matrix
  DATA_INTEGER(N_year);
  DATA_IVECTOR(Year_ID);
  DATA_IVECTOR(Year_ID_pred);
  DATA_INTEGER(Likconfig);

  // Parameters
  PARAMETER_MATRIX(beta);
  PARAMETER(log_sigma);
  PARAMETER_VECTOR(year);
  PARAMETER(log_sigma_year);

  Type sigma=exp(log_sigma);
  Type sigma_year=exp(log_sigma_year);

  Type nll=0;
  // nll = -Nthres*log(Nthres);	// the constant for the prior not needed in TMB

  vector<Type> mu(N);

  // Define mu
	mu=X*beta;

	// add the year random effect to the mus
	for(int n=0;n<N;n++){
	    mu(n) = mu(n) + year(Year_ID(n));
	  }


  // Define the gamma parameters (when needed)
  vector<Type> shape(N);
  vector<Type> scale(N);
	for (int n=0;n<N;n++){
		shape(n)= pow(mu(n),2)/pow(sigma,2);
		scale(n)= pow(sigma,2)/mu(n);
	}

  // vector<Type> lp_e(Nthres+1);
  // vector<Type> lp_l(Nthres+1);
  vector<Type> LL(N);

	for (int thr=0;thr<Nthres;thr++) {
		for (int n=0;n<N;n++){
		// if (thres_cov(n) < (thresh(thr)+is_from_west(n)*mean_diff_tag_area)){
			if (Likconfig ==0){
				nll -= dnorm(y(n), mu(n), sigma, TRUE);			// for normal distribution on log scale
				LL(n)=dnorm(y(n), mu(n), sigma, TRUE);
			}
			if (Likconfig ==1){
				nll -= dgamma(y(n), shape(n), scale(n), TRUE);		// for gamma distribution on log scale
				LL(n)=dgamma(y(n), shape(n), scale(n), TRUE);
			}
		}
	}

	// Add the random effect to the NLL
	  for (int yr=0; yr<N_year; yr++){
	    nll -= dnorm(year(yr), Type(0.0), sigma_year, TRUE);
	  }

  // ============ Outputs =============

  // Prediction
  vector<Type> mu_pred(N_pred);
	mu_pred=X_pred*beta;

	for(int n=0;n<N_pred;n++){
	    mu_pred(n) = mu_pred(n) + year(Year_ID_pred(n));
	  }


  // Model parameters
  REPORT(beta);
  REPORT(sigma);
  REPORT(sigma_year);
  REPORT(LL);
  ADREPORT(mu);
  ADREPORT(mu_pred);

  //--------------------------------------------------------------------

  return nll;

}
