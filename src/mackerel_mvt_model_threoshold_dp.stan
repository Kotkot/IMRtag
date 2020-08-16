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
  int thresh_start[Nthres+1];
  int thresh_end[Nthres+1];
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
  vector<lower=0>[K] sigma;  // scales of mixture components
}

transformed parameters {
  matrix[N,K] mu = X * beta;

  vector[Nthres] lp;
  {
    vector[Nthres + 1] lp_e;
    vector[Nthres + 1] lp_l;
    lp_e[1] = 0;
    lp_l[1] = 0;
    for (n in thresh_start[1]:thresh_end[1]){
      lp_e[1] += normal_lpdf(y[n] | mu[n,1], sigma[1]);  // before thr1
      lp_l[1] += normal_lpdf(y[n] | mu[n,2], sigma[2]);  // after thr2
    }  // Now define the mean travel distance for each mixture component

    for (t in 1:Nthres){
      lp_e[t+1] = lp_e[t];
      lp_l[t+1] = lp_l[t];
      for (n in thresh_start[t+1]:thresh_end[t+1]){
        lp_e[t + 1] += normal_lpdf(y[n] | mu[n,1], sigma[1]);  // before thr1
        lp_l[t + 1] += normal_lpdf(y[n] | mu[n,2], sigma[2]);  // after thr2
      }  // Now define the mean travel distance for each mixture component
    }

    lp = rep_vector(log_unif + lp_l[Nthres + 1], Nthres) + head(lp_e, Nthres) - head(lp_l, Nthres);
  }

  // lp = rep_vector(log_unif, Nthres);
  //
  // // Now define the mean travel distance for each mixture component
  //   mu = X * beta;
  //
  // for (thr in 1:Nthres) {
  //   for (n in 1:N)
  //     lp[thr] += normal_lpdf(y[n] | y[n] < (thresh[thr]+is_from_west[n]*mean_diff_tag_area) ? mu[n,1] : mu[n,2], sigma);
  // }

}

model {
  sigma ~ cauchy(0,1);
  target += log_sum_exp(lp);
}

generated quantities {
  matrix[N,K] y_gen;
  for (n in 1:N){
    for (t in 1:K) {
      y_gen[n,t] = normal_rng(mu[n,t], sigma[t]);
    }
  }
}

