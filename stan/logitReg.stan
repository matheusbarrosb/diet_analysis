data {
  int<lower=0>          N;       // NUMBER OF DATA POINTS
  int<lower=0>          S;       // NUMBER OF SITES
  vector[N]             TL;      // TOTAL LENGTH, CONTINUOUS COVARIATE
  int<lower=0,upper=S>  site[N]; // CATEGORICAL COVARIATE
  int<lower=0,upper=1>  y[N];    // RESPONSE
}
parameters {
  real alpha;         // INTERCEPT
  real beta;          // EFFECT OF TL
  real beta_site[S];  // EFFECT OF SITE
}
model {
  // PRIORS
  alpha     ~ normal(0, 1);
  beta      ~ normal(0, 1);
  beta_site ~ normal(0, 1);
  
  // LIKELIHOOD
  for (i in 1:N) {
      y[i] ~ bernoulli_logit(alpha + beta*TL[i] + beta_site[site[i]]);
  }
}
generated quantities {
  vector<lower=0,upper=1>[S] p;
  
  for (i in 1:N) {
    for (k in 1:S) {
      p[k] = inv_logit(alpha + beta_site[k] + beta*TL[i]);
    }
  }
}

