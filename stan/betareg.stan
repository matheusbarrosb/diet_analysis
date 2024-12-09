data {
  
  int N;                                // NUMBER OF OBSERVATIONS
  int S;                                // NUMBER OF SITES
  int site[N];                          // CATEGORICAL PREDICTOR
  vector<lower=0,upper=1>[N] GF;        // RESPONSE, GUT FULNESS (%)
  real TL[N];                           // CONTINUOUS PREDICTOR
  
}
parameters {
  
  real alpha;                           // OVERALL GF MEAN, UNTRANSFORMED
  vector[S] beta;                       // DEFLECTIONS, UNTRANSFORMED
  vector<lower=0>[S] sigma_beta;        // DEFLECTION STANDARD DEVIATIONS
  real betaTL;                          // EFFECT OF CONTINUOS COVARIATE
  real<lower=0>phi;                     // DISPERSION PARAMETER
  
}
transformed parameters {
  
  vector[N] A;                          // BETA SHAPE PARAMETER 1
  vector[N] B;                          // BETA SHAPE PARAMETER 2
  vector[N] LP;                         // UNTRANSFORMED LINEAR PREDICTOR
  vector[N] mu;                         // TRANSFORMED LINEAR PREDICTOR
  
  for (i in 1:N) {
    for (s in 1:S) {
      
      LP[i] = alpha + beta[site[i]] + betaTL*TL[i];    // GLM
      
    }
  }
  
  for (i in 1:N) {
    
    mu[i] = inv_logit(LP[i]);           // CONVERT GLM TO 0-1 SCALE
    A[i] = mu[i] * phi;
    B[i] = (1 - mu[i])*phi;
    
  }
}

model {
  
// LIKELIHOOD COMPONENT --------------------------------------------------------
  
  target += beta_lpdf(GF | A , B); # EXPLICIT STATEMENT

// PRIOR COMPONENT -------------------------------------------------------------
  
  phi        ~ cauchy(0,10);           // OTHER OPTION IS A BROAD INV_GAMMA
  alpha      ~ normal(0,2);               
  beta       ~ normal(0,sigma_beta);
  sigma_beta ~ cauchy(0,1);
  betaTL     ~ normal(0,2);
  
}

generated quantities {
  
  vector[S] y_hat;              // SITE MEANS
  vector[S] dflc;               // SCALED DEFLECTION PARAMETERS
  vector[S] cnt[(S-1)];         // CONTRASTS FOR PAIRWISE COMPARISONS
  real      mu_GF;

  for (i in 1:N) {
    for (s in 1:S) {
      
      real mu_hat[S];
      real A_hat[S];
      real B_hat[S];
      mu_hat[s] = inv_logit(alpha + beta[s]);
      A_hat[s]  = mu_hat[s] * phi;
      B_hat[s]  = (1.0 - mu_hat[s]) * phi;
      y_hat[s]  = beta_rng(A_hat[s], B_hat[s]); // POSTERIOR PREDICTIVE SITE MEANS
      
    }
  }
  
  mu_GF = mean(y_hat);      // MEAN OF MODELLED GUT FULLNESS
  
  for (s in 1:S) {

    dflc[s] = y_hat[s] - mu_GF;          // CALCULATE DEFLECTIONS FROM BASELINE
    
  }
  
  for (s in 1:S) {
    for (j in 1:(s-1)) {
      
      cnt[j,s] = y_hat[j] - y_hat[s];      // PAIRWISE COMPARISONS
      // NUMBER OF COMPARISONS = K*(K-1), WHERE K =. No OF GROUPS
    }
  }
  
  // // PRIOR PREDICTIVE CHECK
  real alpha_hat      = normal_rng(0,2);
  real beta_TL_hat    = normal_rng(0,2);
  real phi_hat        = cauchy_rng(0,10);
  
}