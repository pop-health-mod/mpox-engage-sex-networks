# Libraries ----
library(tidyverse)
library(data.table)
source("./src/plot.R")

# Source NB regression fit ----
# only run if necessary
if(!exists("fit_bayes_ls")){
  source("./21_distr_fit_negbin_coef.R")
}

# Probability mass function ----
## density and PMF computations already performed in Stan

## Collapse PMF into single dataset ----
pmf_iter <- create_city_list(CITIES_DATAPTS)
pmf_wt_by_city <- create_city_list(CITIES_DATAPTS)

if(DO_AGE){
for(cur_city in CITIES_DATAPTS){
  
  for(a in AGES){
  # extract PMF iterations from stan
  index_age <- which(AGES == a)
  pmf_tmp_a <- extract(fit_bayes_ls[[cur_city]], pars = "pmf")$pmf[, index_age, ]
  pmf_iter[[cur_city]][[a]] <- pmf_tmp_a
  
  # get credible intervals
  cred_l_a <- vector("double", ncol(pmf_tmp_a))
  cred_u_a <- vector("double", ncol(pmf_tmp_a))
  
  for(i in 1:ncol(pmf_tmp_a)){
    cred_l_a[i] <- quantile(pmf_tmp_a[ ,i], .025)
    cred_u_a[i] <- quantile(pmf_tmp_a[ ,i], .975)
  }
  
  # create tibble with each city and time period
  pmf_wt_by_city[[cur_city]][[a]] <- tibble(data_pt = cur_city,
                                       age_grp = a,
                                       y_pred = 0:300,
                                       mean = colSums(pmf_tmp_a) / nrow(pmf_tmp_a),
                                       cr.i_low = cred_l_a,
                                       cr.i_upp = cred_u_a)
  
  rm(cred_l_a, cred_u_a)}
  
  pmf_wt_by_city[[cur_city]] <- bind_rows(pmf_wt_by_city[[cur_city]])}
  
  rm(pmf_tmp_a)}else{
  for(cur_city in CITIES_DATAPTS){
    # extract PMF iterations from stan
    pmf_tmp <- extract(fit_bayes_ls[[cur_city]], pars = "pmf")$pmf
    pmf_iter[[cur_city]] <- pmf_tmp
    
    # get credible intervals
    cred_l <- vector("double", ncol(pmf_tmp))
    cred_u <- vector("double", ncol(pmf_tmp))
    
    for(i in 1:ncol(pmf_tmp)){
      cred_l[i] <- quantile(pmf_tmp[, i], .025)
      cred_u[i] <- quantile(pmf_tmp[, i], .975)
    }
    
    # create tibble with each city and time period
    pmf_wt_by_city[[cur_city]] <- tibble(data_pt = cur_city,
                                         y_pred = 0:300,
                                         mean = colSums(pmf_tmp) / nrow(pmf_tmp),
                                         cr.i_low = cred_l,
                                         cr.i_upp = cred_u)
    
    rm(cred_l, cred_u)
  }
  rm(pmf_tmp)
}

# collapse into
pmf_wt_by_city <- bind_rows(pmf_wt_by_city)
pmf_wt_by_city$data_pt <- factor(pmf_wt_by_city$data_pt, levels = CITIES_DATAPTS)

## Verify PMF posterior distributions ----
# verify results by looking at the mean number of partners
ifelse(DO_AGE, group_var <- c("age_grp", "data_pt"), group_var <- "data_pt")
data_mean_nb_partn <- pmf_wt_by_city %>%
  group_by(across(all_of(group_var))) %>%
  summarize(mean_wt = sum(y_pred * mean),
            cr.i_low = sum(y_pred * cr.i_low),
            cr.i_upp = sum(y_pred * cr.i_upp),
            .groups = "drop") %>% 
  mutate(type = "neg. bin.")

data_mean_nb_partn

# verify that weights sum up to 1
pmf_wt_by_city %>% 
  group_by(across(all_of(group_var))) %>% 
  summarize(dens_ttl = sum(mean), 
            .groups = "drop") 

if(DO_AGE){# save full fitted pmf
  write.csv(pmf_wt_by_city,
            sprintf("%s/pmf_weighted_all_partn.csv", out_distr_path),
            row.names = F)
  stop("cdf was not implemented for age-specific output")}

# Cumulative density function ----

# compute fitted CDF in each datapoint (similar procedure as for pmf)
# make lists to hold data for each city
cdf_wt_by_city <- create_city_list(CITIES_DATAPTS)
for(cur_city in CITIES_DATAPTS){
  # extract PMF iterations from data
  pmf_tmp <- pmf_iter[[cur_city]]
  
  # empty CDF matrix
  cdf_tmp <- matrix(0, nrow = nrow(pmf_tmp), ncol = ncol(pmf_tmp))
  
  # sum PMF starting from long tail towards 0
  for(i in ncol(pmf_tmp):1){
    if(i == ncol(pmf_tmp)){
      cdf_tmp[, i] <- pmf_tmp[, i]
    } else {
      cdf_tmp[, i] <- pmf_tmp[, i] + cdf_tmp[, i+1]
    }
  }
  
  # get credible intervals
  cred_l <- vector("double", ncol(cdf_tmp))
  cred_u <- vector("double", ncol(cdf_tmp))
  
  for(i in 1:ncol(cdf_tmp)){
    cred_l[i] <- quantile(cdf_tmp[, i], .025)
    cred_u[i] <- quantile(cdf_tmp[, i], .975)
  }
  
  # create tibble with each city and time period
  cdf_wt_by_city[[cur_city]] <- tibble(data_pt = cur_city,
                                       y_pred = 0:300,
                                       mean = colSums(cdf_tmp) / nrow(cdf_tmp),
                                       cr.i_low = cred_l,
                                       cr.i_upp = cred_u)
  
  rm(cred_l, cred_u)
}
rm(pmf_tmp, cdf_tmp)

# collapse into a single dataframe
cdf_wt_by_city <- bind_rows(cdf_wt_by_city)

cdf_wt_by_city <- cdf_wt_by_city %>% 
  mutate(city = gsub("-Pre-Pandemic|-Pandemic|-Post-Restrictions", "", data_pt),
         time_pt = gsub("Montreal-|Toronto-|Vancouver-", "", data_pt),
         .before = 1)

# factorize time points to order correctly
cdf_wt_by_city$time_pt <- factor(cdf_wt_by_city$time_pt, 
                                 levels = c("Pre-Pandemic", "Pandemic", "Post-Restrictions"))

# Output tables (PMF and CDF) ----
# save iterations (only for main analyses)
saveRDS(pmf_iter, sprintf("./out/pmf_stan_iterations_%s.rds", out_distr_pref))

# save full fitted pmf
write.csv(pmf_wt_by_city,
          sprintf("%s/pmf_weighted_all_partn.csv", out_distr_path),
          row.names = F)

# save full fitted cdf
write.csv(cdf_wt_by_city,
          sprintf("%s/cdf_weighted_all_partn.csv", out_distr_path),
          row.names = F)