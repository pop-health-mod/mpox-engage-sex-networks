// authors: Jorge Luis Flores Anato, Fanyu Xiu
// https://mc-stan.org/docs/functions-reference/neg-binom-2-log-glm.html
// https://mc-stan.org/docs/stan-users-guide/prediction-forecasting-and-backcasting.html
data {
  int<lower=0> N;            // number of data items at each data_pt
  int<lower=0> N_aggr;       // number of data items for prediction (city-wide)
  int<lower=0> N_aggr_ah;    // number of combinations of covariates (all binary) for each age-hiv group (should be the same for each data_pt)
  int<lower=0> N_ah;         // number of age-hiv groups
  
  int<lower=0> K;            // number of predictors
  matrix[N, K] x;            // predictor matrix (model fitting)
  int<lower=0> y[N];         // outcome vector

  matrix[N_aggr, K] x_aggr;  // predictor matrix for aggregated x's
  matrix[N_aggr_ah, K] x_aggr_ah[N_ah];  // an array containing predictor matrix for aggregated x's (for prediction)
  vector[N_aggr] ipc_rds_w;              // IPC-RDS weights (sums for each group of individuals)
  vector[N_aggr_ah] ipc_rds_w_ah[N_ah];  // matrix of (IPC-)RDS weights (sums for each combo-group of individuals) for each age-hiv group
 
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
    vector[N_aggr] y_pred;             // model predicted means/expected values for each (group of) participant
    vector[x_end+1] pmf;               // population-level (weighted) PMF
    
    vector[N_aggr_ah] y_pred_ah[N_ah]; // array of predicted means/expected values for each combo of participant for each age-hiv group
    vector[x_end+1] pmf_ah[N_ah];      // array of population-level (weighted) PMF for each age-hiv group
    
    // Compute PMF (city-wide) ----
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
    
    // Compute PMF (by age and hiv group) ----
    for(a in 1:N_ah){ // looping over age-hiv groups
      // (1) first get the model predicted log-mean for each combo of participants
      y_pred_ah[a] = x_aggr_ah[a] * beta + alpha;
      
      // (2) then, using that mean, compute the posterior density for each combinationo of participants
      matrix[N_aggr_ah, x_end+1] pmf_ind_ah;
      matrix[N_aggr_ah, x_end+1] pmf_ind_wt_ah;
      
      for(i in 1:N_aggr_ah){
        for(j in 0:x_end){ // nb: Stan is indexed at 1, we want to compute for x = {0,1,...,x_end}
          pmf_ind_ah[i, j+1] = neg_binomial_2_log_lpmf(j | y_pred_ah[a][i], phi);
        }
        pmf_ind_wt_ah[i, ] = exp(pmf_ind_ah[i, ]) * ipc_rds_w_ah[a][i];
      }
      
      // (3) within each iteration, compute the group-wide PMF as the weighted mean of individual-combo densities
      //     by computing P(Y=k) for k = {1,2,...,x_end}
      real sum_wts_ah = sum(ipc_rds_w_ah[a]);
      
      for(j in 0:x_end){
        pmf_ah[a][j+1] = sum(pmf_ind_wt_ah[, j+1]) / sum_wts_ah;
      }
    } // close for age-hiv groups looping
}
