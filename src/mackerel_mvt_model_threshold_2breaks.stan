////////////////////////////////////////////////////////////////
// Threshold model to analyse the within-year mackerela travel rate
// Author: Kotaro Ono
// Version: 1
// Detail:
// 2 brekpoint analysis (break point is some based on travel distance)
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
  int Nthres_tot = Nthres*(Nthres-1)/2;
  log_unif = -log(Nthres_tot);
}

parameters {
  matrix[Nx,K] beta;             // locations of mixture components
  vector<lower=0>[K] sigma;  // scales of mixture components
}

transformed parameters {
  // vector[Nthres_tot] lp;
  matrix[N,K] mu;
  // int counter = 1;

  // lp = rep_vector(log_unif, Nthres_tot);

  // Now define the mean travel distance for each mixture component
    mu = X * beta;

  // for (thr1 in 1:(Nthres-1) {
  //   for (thr2 in (thr1+1):Nthres) {
  //     for (n in 1:N){
  //        lp[counter] += normal_lpdf(y[n] | y[n] < (thresh[thr1]+is_from_west[n]*mean_diff_tag_area) ? mu[n,1] : ((y[n] < thresh[thr2]+is_from_west[n]*mean_diff_tag_area) ? mu[n,2]: mu[n,3]), sigma);
  //     }
  //   counter += counter;
  //   }
  // }

}

model {
  vector[Nthres_tot] lp;
  int counter = 1;

  lp = rep_vector(log_unif, Nthres_tot);

  for (thr1 in 1:(Nthres-1)) {
    for (thr2 in (thr1+1):Nthres) {
      for (n in 1:N){
         lp[counter] += normal_lpdf(y[n] | y[n] < (thresh[thr1]+is_from_west[n]*mean_diff_tag_area) ? mu[n,1] : ((y[n] < thresh[thr2]+is_from_west[n]*mean_diff_tag_area) ? mu[n,2]: mu[n,3]), sigma);
      }
    counter += 1;
    }
  }


  sigma ~ cauchy(0,1);
  target += log_sum_exp(lp);
}

// generated quantities {
  // matrix[N,Nthres_tot] y_gen;

  // int count;
  // for (n in 1:N){
    // count = 1;
    // for (thr1 in 1:(Nthres-1)) {
      // for (thr2 in (thr1+1):Nthres) {
        // y_gen[n,count] = normal_rng(y[n] < (thresh[thr1]+is_from_west[n]*mean_diff_tag_area) ? mu[n,1] : ((y[n] < thresh[thr2]+is_from_west[n]*mean_diff_tag_area) ? mu[n,2]: mu[n,3]), sigma);
        // count += 1;
      // }
    // }
  // }
// }

