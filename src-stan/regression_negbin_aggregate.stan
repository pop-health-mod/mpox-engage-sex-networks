// https://mc-stan.org/docs/functions-reference/neg-binom-2-log-glm.html
// https://mc-stan.org/docs/stan-users-guide/prediction-forecasting-and-backcasting.html
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
    // int y_hat[N];                         // model posterior predictive distributions
    vector[N_aggr] y_pred;                // model predicted means/expected values for each (group of) participant
    vector[x_end+1] pmf;                  // population-level (weighted) PMF
    
    // y_hat = neg_binomial_2_log_rng(x * beta + alpha, phi);
    
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
          pmf_ind[i, j+1] = neg_binomial_2_log_lpmf(j | y_pred[i], phi);
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
