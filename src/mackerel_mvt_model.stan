////////////////////////////////////////////////////////////////
// Mixture model to analyse the within-year mackerela travel rate
// Author: Kotaro Ono
// Version: 1
// Detail:
// Applies a mixture model of fast vs slower moving mackerel (the number of mixture is to be changed to see
// what is best)
//

data {
  int<lower=1> K;          // number of mixture components
  int<lower=1> N;          // number of data points
  int<lower=1> Nx;         // number of columns in the design matrix
  matrix[N,Nx] X;          // the design matrix
  real y[N];               // observations
}

parameters {
  simplex[K] theta;          // mixing proportions
  matrix[Nx,K] beta;             // locations of mixture components
  vector<lower=0>[K] sigma;  // scales of mixture components
}

transformed parameters {
  matrix[N,K] mu;
  // Now define the mean travel distance for each mixture component
    mu = X * beta;
}

model {
  vector[K] log_theta = log(theta);  // cache log calculation
  sigma ~ lognormal(0, 2);


  for (n in 1:N) {
    vector[K] lps = log_theta;
    for (k in 1:K)
      lps[k] += normal_lpdf(y[n] | mu[k], sigma[k]);
    target += log_sum_exp(lps);
  }
}
