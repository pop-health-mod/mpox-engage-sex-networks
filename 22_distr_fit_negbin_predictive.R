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

# collapse into
pmf_wt_by_city <- bind_rows(pmf_wt_by_city)
pmf_wt_by_city$data_pt <- factor(pmf_wt_by_city$data_pt, levels = CITIES_DATAPTS)

## Verify PMF posterior distributions ----
# verify results by looking at the mean number of partners
data_mean_nb_partn <- pmf_wt_by_city %>%
  group_by(data_pt) %>%
  summarize(mean_wt = sum(y_pred * mean),
            cr.i_low = sum(y_pred * cr.i_low),
            cr.i_upp = sum(y_pred * cr.i_upp),
            .groups = "drop") %>% 
  mutate(type = "neg. bin.")

data_mean_nb_partn

# verify that weights sum up to 1
pmf_wt_by_city %>% 
  group_by(data_pt) %>% 
  summarize(dens_ttl = sum(mean), .groups = "drop")

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
saveRDS(pmf_iter, sprintf("./out/pmf_stan_iterations%s.rds", out_distr_pref))

# save full fitted pmf
write.csv(pmf_wt_by_city,
          sprintf("%s/pmf_weighted_all_partn.csv", out_distr_path),
          row.names = F)

# save full fitted cdf
write.csv(cdf_wt_by_city,
          sprintf("%s/cdf_weighted_all_partn.csv", out_distr_path),
          row.names = F)

# Tail comparison ---------------------------------------------------------

# cut-off above which to compute the CDF for the comparisons
# i.e., compute % of population reporting >=x degree of partners
degree <- 100

## within each city, compare
# pandemic vs pre-,
# post- vs pandemic, and
# post- vs pre-
cdf_comparison_by_city <- create_city_list(CITIES)

for(city_name in CITIES){
  # extract PMF iterations from Stan for each time point
  pmf_pre <-  extract(fit_bayes_ls[[paste(city_name, "Pre-Pandemic", sep = "-")]], pars = "pmf")$pmf
  pmf_pand <- extract(fit_bayes_ls[[paste(city_name, "Pandemic", sep = "-")]], pars = "pmf")$pmf
  pmf_post <- extract(fit_bayes_ls[[paste(city_name, "Post-Restrictions", sep = "-")]], pars = "pmf")$pmf
  
  # compare for >=100 partners
  cdf_comparison_by_city[[city_name]] <-  foreach( cur_iter = 1:nrow(pmf_pre) ) %dopar% {
    compare_timepts(
      pmf_1 = pmf_pre,
      pmf_2 = pmf_pand,
      pmf_3 = pmf_post,
      degree_cutoff = 100
    )
  }
  
  # get the % in each comparison that returned TRUE
  cdf_comparison_by_city[[city_name]] <- bind_rows(cdf_comparison_by_city[[city_name]])
  
  cdf_compare_summ <- apply(cdf_comparison_by_city, MARGIN = 2, mean)
  
  cdf_comparison_by_city[[city_name]] <- tibble(
    city = city_name,
    pre_pand = cdf_compare_summ["pre_pand"],
    pand_post = cdf_compare_summ["pand_post"],
    pre_post = cdf_compare_summ["pre_post"]
  )
  
}
rm(pmf_pre, pmf_pand, pmf_post)

## within the 1st and 3rd time periods, compare
# Toronto vs Montreal
# Toronto vs Vancouver
# Vancouver vs Montreal
cdf_comparison_by_timept <- create_city_list(TIMEPTS)

for(time_pt_name in TIMEPTS){
  # extract PMF iterations from Stan for each time point
  pmf_mtl <- extract(fit_bayes_ls[[paste("Montreal", time_pt_name, sep = "-")]], pars = "pmf")$pmf
  pmf_tor <- extract(fit_bayes_ls[[paste("Toronto", time_pt_name, sep = "-")]], pars = "pmf")$pmf
  pmf_van <- extract(fit_bayes_ls[[paste("Vancouver", time_pt_name, sep = "-")]], pars = "pmf")$pmf
  
  # compare for >=100 partners
  cdf_comparison_by_timept[[time_pt_name]] <-  foreach( cur_iter = 1:nrow(pmf_pre) ) %dopar% {
    compare_cities(
      pmf_mtl = pmf_mtl,
      pmf_tor = pmf_tor,
      pmf_van = pmf_van,
      degree_cutoff = 100
    )
  }
  
  # get the % in each comparison that returned TRUE
  cdf_comparison_by_timept[[time_pt_name]] <- bind_rows(cdf_comparison_by_timept[[time_pt_name]])
  
  cdf_compare_summ <- apply(cdf_comparison_by_timept, MARGIN = 2, mean)
  
  # save results in final list
  cdf_comparison_by_timept[[time_pt_name]] <- tibble(
    time_pt = city_name,
    mtl_tor = cdf_compare_summ["mtl_tor"],
    van_tor = cdf_compare_summ["van_tor"],
    mtl_van = cdf_compare_summ["mtl_van"]
  )
  
}
rm(pmf_mtl, pmf_tor, pmf_van)

