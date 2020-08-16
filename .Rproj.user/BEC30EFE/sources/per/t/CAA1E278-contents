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
  int thresh_start[Nthres+1];
  int thresh_end[Nthres+1];
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
  // vector<lower=0>[K] sigma;  // scales of mixture components
  vector<lower=0>[K] sigma;  // scales of mixture components
}

transformed parameters {
  matrix[N,K] mu = X * beta;
  vector[Nthres] lp;
  {
    vector[Nthres + 1] lp_e;
    vector[Nthres + 1] lp_m;
    vector[Nthres + 1] lp_l;
    lp_e[1] = 0;
    lp_m[1] = 0;
    lp_l[1] = 0;
     for (n in thresh_start[1]:thresh_end[1]){
      lp_e[1] += normal_lpdf(y[n] | mu[n,1], sigma[1]);  // before thr1
      lp_m[1] += normal_lpdf(y[n] | mu[n,2], sigma[2]);  // before thr1
      lp_l[1] += normal_lpdf(y[n] | mu[n,3], sigma[3]);  // after thr2
    }  // Now define the mean travel distance for each mixture component


    for (t in 1:Nthres){
      lp_e[t+1] = lp_e[t];
      lp_m[t+1] = lp_m[t];
      lp_l[t+1] = lp_l[t];
      for (n in thresh_start[t+1]:thresh_end[t+1]){
        lp_e[t + 1] += normal_lpdf(y[n] | mu[n,1], sigma[1]);  // before thr1
        lp_m[t + 1] += normal_lpdf(y[n] | mu[n,2], sigma[2]);  // before thr1
        lp_l[t + 1] += normal_lpdf(y[n] | mu[n,3], sigma[3]);  // after thr2
      }  // Now define the mean travel distance for each mixture component
    }

    for (i in 1:Nthres){
    lp[i]=1;
    }

    for (thr2 in 2:(Nthres-1)) {
      for (thr in 1:thr2) {
         lp[1:thr2] += rep_vector(log_unif + lp_m[thr2 + 1], thr2) + head(lp_e, thr2) - head(lp_m, thr2);
      }
      for (thr in (thr2+1):Nthres) {
         lp[(thr2+1):Nthres] += rep_vector(log_unif, Nthres-thr2) + tail(lp_l, Nthres-thr2) - rep_vector(lp_l[thr2+1], Nthres-thr2);
     }
    }
  }
}


model {
  sigma ~ cauchy(0,1);
  target += log_sum_exp(lp);
}

// generated quantities {
//   matrix[N,Nthres_tot] y_gen;
//
//   int count;
//   for (n in 1:N){
//     count = 1;
//     for (thr1 in 1:(Nthres-1)) {
//       for (thr2 in (thr1+1):Nthres) {
//         y_gen[n,count] = normal_rng(y[n] < (thresh[thr1]+is_from_west[n]*mean_diff_tag_area) ? mu[n,1] : ((y[n] < thresh[thr2]+is_from_west[n]*mean_diff_tag_area) ? mu[n,2]: mu[n,3]), sigma);
//         count += 1;
//       }
//     }
//   }
// }
// generated quantities {
//   matrix[N,K] y_gen;
//   for (n in 1:N){
//     for (t in 1:K) {
//       y_gen[n,t] = normal_rng(mu[n,t], sigma[t]);
//     }
//   }
// }


