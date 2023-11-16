// https://mc-stan.org/docs/functions-reference/neg-binom-2-log-glm.html
// https://mc-stan.org/docs/stan-users-guide/prediction-forecasting-and-backcasting.html
data {
  int<lower=0> N;            // number of data items at each data_pt
  int<lower=0> n_age;        // number of age groups
  
  int<lower=0> K;            // number of predictors
  matrix[N, K] x;            // predictor matrix (for fitting model)
  int<lower=0> y[N];         // outcome vector

  int<lower=0> N_aggr_age; // number of combinations £¨combos) of covariates (all binary) for each age group (should be the same for each data_pt)
  matrix[N_aggr_age, K] x_aggr_age[n_age];  // an array containing predictor matrix for aggregated x's (for prediction)
  vector[N_aggr_age] ipc_rds_w_age[n_age];  // a matrix containing vectors of IPC-RDS weights (sums for each combo-group of individuals) for each age group
 
  int<lower=0> x_end;        // upper bound to compute PMF and CDF values

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
  y ~ neg_binomial_2_log_glm(x, alpha, beta, phi); // fit model
}

generated quantities {

      vector[N_aggr_age] y_pred[n_age]; // an array of predicted means/expected values for each combo of participant for each age group
      vector[x_end+1] pmf[n_age]; // an array of population-level (weighted) PMF for each age group

      // Compute PMF ----
      {

          for(a in 1:n_age){ // looping over age groups
          // (1) first get the model predicted log-mean for each combo of participants
          y_pred[a] = x_aggr_age[a] * beta + alpha;

          // (2) then, using that mean, compute the posterior density for each combo of participants
          matrix[N_aggr_age, x_end+1] pmf_ind_age;
          matrix[N_aggr_age, x_end+1] pmf_ind_wt_age;

          for(i in 1:N_aggr_age){
             for(j in 0:x_end){ // nb: Stan is indexed at 1, we want to compute for x = {0,1,...,x_end}
                pmf_ind_age[i, j+1] = neg_binomial_2_log_lpmf(j | y_pred[a][i], phi);
              }
              pmf_ind_wt_age[i, ] = exp(pmf_ind_age[i, ]) * ipc_rds_w_age[a][i];
              }
      
           // (3) within each iteration, compute the population-wide PMF as the weighted mean of individual-combo densities
           //     by computing P(Y=k) for k = {1,2,...,x_end}

           real sum_wts_age = sum(ipc_rds_w_age[a]);
      
           for(j in 0:x_end){
               pmf[a][j+1] = sum(pmf_ind_wt_age[, j+1]) / sum_wts_age;
            }
        } // close for age groups looping

      } //close for compute pmf
    y_pred = exp(y_pred); // exponentiate into mean
}
