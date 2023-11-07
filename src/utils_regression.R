# Comparison of CDF between fits ----
# compares whether the specified tail section is bigger in pmf_b compared to pmf_a
compare_pmf_tail <- function(pmf_a, pmf_b, degree_cutoff){
  # first compute the size of the tail for both
  cdf_over_a <- sum(pmf_a[degree_cutoff:length(pmf_a)])
  cdf_over_b <- sum(pmf_b[degree_cutoff:length(pmf_b)])
  
  # then return whether cdf_1 is bigger than cdf_0
  return( as.numeric(cdf_over_b > cdf_over_a) )
}

# compare across the 3 different time points,
# where 1: pre-pandemic time period
#       2: pandemic time period
#       3: post-restrictions time period
compare_timepts <- function(pmf_1, pmf_2, pmf_3,
                            degree_cutoff){
  # first nb is the reference time period and second nb is the comparison
  comparison_1_2 <- compare_pmf_tail(pmf_1, pmf_2, degree_cutoff)
  comparison_2_3 <- compare_pmf_tail(pmf_2, pmf_3, degree_cutoff)
  comparison_1_3 <- compare_pmf_tail(pmf_1, pmf_3, degree_cutoff)
  
  # organize comparisons
  df <- data.frame(
    pre_pand = comparison_1_2,   # did activity decrease during the pandemic?
    pand_post = comparison_2_3,  # did activity increase after restrictions lifting?
    pre_post = comparison_1_3    # did activity go back to pre-pandemic (after lifting)?
  )
  
  return(df)
}

# compare across the 3 cities,
compare_cities <- function(pmf_mtl, pmf_tor, pmf_van,
                           degree_cutoff){
  # first nb is the reference city and second nb is the comparison
  comparison_mtl_tor <- compare_pmf_tail(pmf_mtl, pmf_tor, degree_cutoff)
  comparison_van_tor <- compare_pmf_tail(pmf_van, pmf_tor, degree_cutoff)
  comparison_mtl_van <- compare_pmf_tail(pmf_mtl, pmf_van, degree_cutoff)
  
  # organize comparisons
  df <- data.frame(
    mtl_tor = comparison_mtl_tor,   # was the tail 'fatter' in Toronto vs Montreal?
    van_tor = comparison_van_tor,   #          "            in Toronto vs Vancouver?
    mtl_van = comparison_mtl_van    #          "            in Vancouver vs Montreal?
  )
  
  return(df)
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
