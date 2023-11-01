// adapted from brms 2.20.4
functions {
  /* zero-inflated negative binomial log-PDF of a single response
   * Args:
   *   y: the response value
   *   mu: mean parameter of negative binomial distribution
   *   phi: shape parameter of negative binomial distribution
   *   zi: zero-inflation probability
   * Returns:
   *   a scalar to be added to the log posterior
   */
  real zero_inflated_neg_binomial_lpmf(int y, real mu, real phi,
                                       real zi) {
    if (y == 0) {
      return log_sum_exp(bernoulli_lpmf(1 | zi),
                         bernoulli_lpmf(0 | zi) +
                         neg_binomial_2_lpmf(0 | mu, phi));
    } else {
      return bernoulli_lpmf(0 | zi) +
             neg_binomial_2_lpmf(y | mu, phi);
    }
  }
  /* zero-inflated negative binomial log-PDF of a single response
   * logit parameterization of the zero-inflation part
   * Args:
   *   y: the response value
   *   mu: mean parameter of negative binomial distribution
   *   phi: shape parameter of negative binomial distribution
   *   zi: linear predictor for zero-inflation part
   * Returns:
   *   a scalar to be added to the log posterior
   */
  real zero_inflated_neg_binomial_logit_lpmf(int y, real mu,
                                             real phi, real zi) {
    if (y == 0) {
      return log_sum_exp(bernoulli_logit_lpmf(1 | zi),
                         bernoulli_logit_lpmf(0 | zi) +
                         neg_binomial_2_lpmf(0 | mu, phi));
    } else {
      return bernoulli_logit_lpmf(0 | zi) +
             neg_binomial_2_lpmf(y | mu, phi);
    }
  }
  /* zero-inflated negative binomial log-PDF of a single response
   * log parameterization for the negative binomial part
   * Args:
   *   y: the response value
   *   eta: linear predictor for negative binomial distribution
   *   phi: shape parameter of negative binomial distribution
   *   zi: zero-inflation probability
   * Returns:
   *   a scalar to be added to the log posterior
   */
  real zero_inflated_neg_binomial_log_lpmf(int y, real eta,
                                           real phi, real zi) {
    if (y == 0) {
      return log_sum_exp(bernoulli_lpmf(1 | zi),
                         bernoulli_lpmf(0 | zi) +
                         neg_binomial_2_log_lpmf(0 | eta, phi));
    } else {
      return bernoulli_lpmf(0 | zi) +
             neg_binomial_2_log_lpmf(y | eta, phi);
    }
  }
  /* zero-inflated negative binomial log-PDF of a single response
   * log parameterization for the negative binomial part
   * logit parameterization of the zero-inflation part
   * Args:
   *   y: the response value
   *   eta: linear predictor for negative binomial distribution
   *   phi: shape parameter of negative binomial distribution
   *   zi: linear predictor for zero-inflation part
   * Returns:
   *   a scalar to be added to the log posterior
   */
  real zero_inflated_neg_binomial_log_logit_lpmf(int y, real eta,
                                                 real phi, real zi) {
    if (y == 0) {
      return log_sum_exp(bernoulli_logit_lpmf(1 | zi),
                         bernoulli_logit_lpmf(0 | zi) +
                         neg_binomial_2_log_lpmf(0 | eta, phi));
    } else {
      return bernoulli_logit_lpmf(0 | zi) +
             neg_binomial_2_log_lpmf(y | eta, phi);
    }
  }
  // zero_inflated negative binomial log-CCDF and log-CDF functions
  real zero_inflated_neg_binomial_lccdf(int y, real mu, real phi, real hu) {
    return bernoulli_lpmf(0 | hu) + neg_binomial_2_lccdf(y | mu, phi);
  }
  real zero_inflated_neg_binomial_lcdf(int y, real mu, real phi, real hu) {
    return log1m_exp(zero_inflated_neg_binomial_lccdf(y | mu, phi, hu));
  }
}

data {
  int<lower=0> N;            // number of data items
  int<lower=0> N_aggr;       // number of data items for prediction
  int<lower=0> K;            // number of predictors
  matrix[N, K] x;            // predictor matrix
  matrix[N_aggr, K] x_aggr;  // predictor matrix for aggregated x's
  int<lower=0> y[N];         // outcome vector
  
  int<lower=0> x_end;        // upper bound to compute PMF and CDF values
  vector[N_aggr] ipc_rds_w;  // IPC-RDS weights (sums for each group of individuals)
}

// not necessary since we're not using any continuous coefficients
// transformed data {
//   matrix[N, Kc] Xc;  // centered version of X without an intercept
//   vector[Kc] means_X;  // column means of X before centering
//   for (i in 2:K) {
//     means_X[i - 1] = mean(X[, i]);
//     Xc[, i - 1] = X[, i] - means_X[i - 1];
//   }
// }

parameters {
  real alpha;             // intercept
  vector[K] beta;         // predictors
  // real<lower=0> phi;      // neg. binomial dispersion parameter
  
  // from brms
  // vector[Kc] b;  // regression coefficients
  // real Intercept;  // temporary intercept for centered predictors
  real<lower=0> shape;  // shape parameter
  real<lower=0,upper=1> zi;  // zero-inflation probability
  // ADD ZERO INFLTED UNIFORM PROB
}

// from brms
transformed parameters {
  real lprior = 0;  // prior contributions to the log posterior
  lprior += student_t_lpdf(alpha | 3, 1.1, 2.5);
  lprior += gamma_lpdf(shape | 0.01, 0.01);
  lprior += beta_lpdf(zi | 1, 1);
}

model {
  // priors:
  // for negative binomial model
  // phi ~ cauchy(0, 5);
  alpha ~ normal(0, 10);
  for(k in 1:K){
    beta[k] ~ normal(0, 10);
  }
  
  // likelihood including constants
  // initialize linear predictor term
  vector[N] mu = rep_vector(0.0, N);
  mu += alpha + x * beta;
  for (n in 1:N) {
    target += zero_inflated_neg_binomial_log_lpmf(y[n] | mu[n], shape, zi);
  }
    
  // priors including constants
  target += lprior;
}

generated quantities {
    int y_hat[N];                         // model posterior predictive distributions
    vector[N_aggr] y_pred;                // model predicted means/expected values for each (group of) participant
    vector[x_end+1] pmf;                  // population-level (weighted) PMF

    // y_hat = neg_binomial_2_log_rng(x * beta + alpha, shape);

    // Compute PMF ----
    {
      // (1) first get the model predicted log-mean for each (group of) participants
      y_pred = x_aggr * beta + alpha;

      // (2) then, using that mean, compute the posterior density for each individual
      matrix[N_aggr, x_end+1] pmf_ind;
      matrix[N_aggr, x_end+1] pmf_ind_wt;

      for(i in 1:N_aggr){
        // if(i % 10 == 0){
        //    print("observation: ", i, " of ", N_aggr);
        // }

        for(j in 0:x_end){ // nb: Stan is indexed at 1, we want to compute for x = {0,1,...,x_end}
          // pmf_ind[i, j+1] = neg_binomial_2_log_lpmf(j | y_pred[i], phi);
          // pmf_ind[i, j+1] = zero_inflated_neg_binomial_lpmf(j | y_pred[i], shape, zi);
          pmf_ind[i, j+1] = zero_inflated_neg_binomial_log_lpmf(j | y_pred[i], shape, zi);
        }

        pmf_ind_wt[i, ] = exp(pmf_ind[i, ]) * ipc_rds_w[i];
      }

      // (3) within each iteration, compute the population-wide PMF as the weighted mean of individual densities
      //     by computing P(Y=k) for k = {1,2,...,x_end}
      real sum_wts = sum(ipc_rds_w);

      for(j in 0:x_end){
        pmf[j+1] = sum(pmf_ind_wt[, j+1]) / sum_wts;
      }
    }
    y_pred = exp(y_pred); // exponentiate into mean
}
