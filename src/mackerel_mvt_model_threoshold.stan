////////////////////////////////////////////////////////////////
// Threshold model to analyse the within-year mackerela travel rate
// Author: Kotaro Ono
// Version: 1
// Detail:
// 1 brekpoint analysis (break point is some based on travel distance)
//

data {
  int<lower=1> K;          // number of mixture components
  int<lower=1> N;          // number of data points
  int<lower=1> Nx;         // number of columns in the design matrix
  matrix[N,Nx] X;          // the design matrix
  int<lower=1> Nthres;
  real thresh[Nthres];
  real mean_diff_tag_area;
  real is_from_west[N];
  real y[N];               // observations
}

transformed data {
  real log_unif;
  log_unif = -log(Nthres);
}

parameters {
  matrix[Nx,K] beta;             // locations of mixture components
  real<lower=0> sigma;  // scales of mixture components
}

transformed parameters {
  vector[Nthres] lp;
  matrix[N,K] mu;

  lp = rep_vector(log_unif, Nthres);

  // Now define the mean travel distance for each mixture component
    mu = X * beta;

  for (thr in 1:Nthres) {
    for (n in 1:N)
      lp[thr] += normal_lpdf(y[n] | y[n] < (thresh[thr]+is_from_west[n]*mean_diff_tag_area) ? mu[n,1] : mu[n,2], sigma);
  }

}

model {
  sigma ~ cauchy(0,1);
  target += log_sum_exp(lp);
}

generated quantities {
  matrix[N,Nthres] y_gen;
  for (n in 1:N){
    for (thr in 1:Nthres) {
      y_gen[n,thr] = normal_rng(y[n] < (thresh[thr]+is_from_west[n]*mean_diff_tag_area) ? mu[n,1] : mu[n,2], sigma);
    }
  }
}

