functions {
  int first_capture(int[] y_i) {
    for (k in 1:size(y_i))
      if (y_i[k])
        return k;
    return 0;
  }
  int last_capture(int[] y_i) {
    for (k_rev in 0:(size(y_i) - 1)) {
      int k;
      k = size(y_i) - k_rev;
      if (y_i[k])
        return k;
    }
    return 0;
  }
  vector prob_uncaptured(int T, vector p, vector phi) {
    vector[T] chi;
    chi[T] = 1.0;
    for (t in 1:(T - 1)) {
      int t_curr;
      int t_next;
      t_curr = T - t;
      t_next = t_curr + 1;
      chi[t_curr] = (1 - phi[t_curr])
                     + phi[t_curr]
                       * (1 - p[t_next])
                       * chi[t_next];
    }
    return chi;
  }
}

data {
  int<lower=2> T;
  int<lower=0> I;
  int<lower=0,upper=1> y[I, T];
  int<lower=1> lh[I];
  int<lower=1> n_lh;
}

transformed data {
  int<lower=0,upper=T> first[I];
  int<lower=0,upper=T> last[I];
  vector<lower=0,upper=I>[T] n_captured;
  for (i in 1:I)
    first[i] = first_capture(y[i]);
  for (i in 1:I)
    last[i] = last_capture(y[i]);
  n_captured = rep_vector(0, T);
  for (t in 1:T)
    for (i in 1:I)
      if (y[i, t])
        n_captured[t] = n_captured[t] + 1;
}

parameters {
  matrix<lower=0,upper=1>[T-1,n_lh] phi;
  matrix<lower=0,upper=1>[T,n_lh] p;
}
transformed parameters {
  matrix<lower=0,upper=1>[T,n_lh] chi;
  for(i in 1:n_lh) {
    chi[,i] = prob_uncaptured(T,p[,i],phi[,i]);
  }
}

model {
  for (i in 1:I) {
    if (first[i] > 0) {
      for (t in (first[i]+1):last[i]) {
        1 ~ bernoulli(phi[t-1,lh[i]]);
        y[i, t] ~ bernoulli(p[t,lh[i]]);
      }
      1 ~ bernoulli(chi[last[i],lh[i]]);
    }
  }
}

generated quantities {
  vector[2] beta;
  real delta_llrtp_grj;
  real delta_grj_boj;
  real delta_llrtp_boj;
  real delta_llrtp_bon;
  real delta_llrtp_gra;
  real delta_sar_bon;
  real delta_sar_gra;

  for(i in 1:n_lh) {
    beta[i] = phi[T-1,i] * p[T,i];
  }

  if(n_lh > 1) {
    delta_llrtp_grj = log(phi[1,1] / phi[1,2]);
    delta_grj_boj = log(phi[2,1] / phi[2,2]);
    delta_llrtp_boj = log((phi[1,1] * phi[2,1]) / (phi[1,2] * phi[2,2]));
    delta_llrtp_bon = log((phi[1,1] * phi[2,1] * phi[3,1]) / (phi[1,2] * phi[2,2] * phi[3,2]));
    delta_llrtp_gra = log((phi[1,1] * phi[2,1] * phi[3,1] * phi[4,1]) / (phi[1,2] * phi[2,2] * phi[3,2] * phi[4,2]));
    delta_sar_bon = log((phi[2,1] * phi[3,1]) / (phi[2,2] * phi[3,2]));
    delta_sar_gra = log((phi[2,1] * phi[3,1] * phi[4,1]) / (phi[2,2] * phi[3,2] * phi[4,2]));
  }


}
