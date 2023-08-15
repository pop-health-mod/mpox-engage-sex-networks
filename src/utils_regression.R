library(Rcpp)
## TODO CLEAN UP FILE (unused functions, clarify code)
# Process observed data ----
#' compute observed distribution
compute_density_obs <- function(data, var_city = "city",wt="wt_rds_norm"){
  # weighted and unweighted frequencies
  density_observed <- data %>% 
    group_by(get(var_city), nb_part_ttl) %>% 
    summarize(nb_ppl = n(), wt_sum = sum(get(wt)), .groups = "drop_last")
  
  names(density_observed)[names(density_observed) == "get(var_city)"] <- var_city
  
  # weighted and unweighted probabilities/densities
  density_observed <- density_observed %>% 
    mutate(pop = sum(nb_ppl), dens_un = nb_ppl / pop,
           pop_wt = sum(wt_sum), dens_wt = wt_sum / pop_wt) %>% 
    ungroup()
  
  return(density_observed[, c(var_city, "nb_part_ttl", "nb_ppl", "pop", "dens_un",
                              "wt_sum", "pop_wt", "dens_wt")])
}

compute_density_obs_new <- function(data, var_city = "city",wt="wt_rds_norm"){
  # weighted and unweighted frequencies
  density_observed <- data %>% 
    group_by(get(var_city), nb_part_new) %>% 
    summarize(nb_ppl = n(), wt_sum = sum(get(wt)), .groups = "drop_last")
  
  names(density_observed)[names(density_observed) == "get(var_city)"] <- var_city
  
  # weighted and unweighted probabilities/densities
  density_observed <- density_observed %>% 
    mutate(pop = sum(nb_ppl), dens_un = nb_ppl / pop,
           pop_wt = sum(wt_sum), dens_wt = wt_sum / pop_wt) %>% 
    ungroup()
  
  return(density_observed[, c(var_city, "nb_part_new", "nb_ppl", "pop", "dens_un",
                              "wt_sum", "pop_wt", "dens_wt")])
}
# Regression result posterior checks ----
#' turn posterior predictive distributions into a data.frame
#' with a row per iteration and a column for each individual's prediction
extract_yhat_long <- function(stan_fit, row_id){
  # rows: iterations, columns: observations
  stan_pred <- extract(stan_fit, pars = "y_hat")$y_hat
  
  # assign names to each observation
  colnames(stan_pred) <- row_id
  
  # turn into a data.frame in long format
  stan_pred <- as_tibble(stan_pred)
  stan_pred <- stan_pred %>% 
    pivot_longer(cols = everything(), names_to = "part_id", values_to = "y_hat")
  
  # add iteration numbering
  stan_pred <- stan_pred %>% 
    group_by(part_id) %>% 
    mutate(iter = 1:n()) %>% 
    ungroup()
  
  # reorder and return
  return(stan_pred[, c("iter", "part_id", "y_hat")])
}

#' provide a summary of the posterior predictive distributions
summarize_yhat <- function(data, cred_int = .95){
  data %>% 
    group_by(part_id) %>% 
    summarize(mean = mean(y_hat),
              mdn = median(y_hat),
              cr.i_low = quantile(y_hat, .5 - (cred_int/2),na.rm=TRUE),
              cr.i_upp = quantile(y_hat, .5 + (cred_int/2),na.rm=TRUE),
              min = min(y_hat),
              max = max(y_hat))
}
