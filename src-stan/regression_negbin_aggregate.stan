// https://mc-stan.org/docs/functions-reference/neg-binom-2-log-glm.html
// https://mc-stan.org/docs/stan-users-guide/prediction-forecasting-and-backcasting.html
data {
  int<lower=0> N;            // number of data items
  int<lower=0> N_aggr;       // number of data items for prediction
  int<lower=0> K;            // number of predictors
  matrix[N, K] x;            // predictor matrix
  matrix[N_aggr, K] x_aggr;  // predictor matrix for aggregated x's
  int<lower=0> y[N];         // outcome vector
}

parameters {
  real alpha;             // intercept
  vector[K] beta;         // predictors
  real<lower=0> phi;      // neg. binomial dispersion parameter
}

model {
  // priors:
  phi ~ cauchy(0, 5);
  alpha ~ normal(0, 10);
  for(k in 1:K){
    beta[k] ~ normal(0, 10);
  }
  y ~ neg_binomial_2_log_glm(x, alpha, beta, phi);
}

generated quantities {
    int y_hat[N];           // model posterior predictive distributions
    vector[N_aggr] y_pred;       // model predicted means/expected values for each
    
    y_hat = neg_binomial_2_log_rng(x * beta + alpha, phi);
    y_pred = exp(x_aggr * beta + alpha);
}
