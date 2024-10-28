data {
  int<lower=0>          N;          // NUMBER OF DATA POINTS
  int<lower=0>          S;          // NUMBER OF SITES
  vector[N]             TL;         // TOTAL LENGTH, CONTINUOUS COVARIATE
  int<lower=0,upper=S>  site[N];    // CATEGORICAL COVARIATE
  int<lower=0,upper=1>  y[N];       // RESPONSE
  real<lower=0>         meanTLs[S]; // FOR POSTERIOR PREDICTION
}
parameters {
  real alpha;         // INTERCEPT, FIXED
  real beta;          // EFFECT OF TL, FIXED
  real beta_site[S];  // EFFECT OF SITE, RANDOM
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
  // POSTERIOR PREDICTIONS, SITE-SPECIFIC PROBABILITIES OF ENCOUNTER
  for (k in 1:S) {
    p[k] = inv_logit(alpha + beta_site[k] + beta*meanTLs[k]); 
  }
}

