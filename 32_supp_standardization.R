# library and data----
library(tidyverse)
library(data.table)
library(rstan)
library(parallel)
library(doParallel)
library(foreach)

source("./src/utils_helper.R")
source("./src/utils_regression.R")
theme_set(theme_bw())

## which variable to use as outcome
outcome_var <- "nb_part_ttl" 
out_distr_path <- "./out/fitted-distr-sens-standard"
out_distr_pref <- "-standard"
fig_path <- "./fig/results-checks-standard"

data_3cities_pre_ipcw <- read_csv("../mpx-engage-params/data-processed/pre_ipcw_3cities.csv")
data_3cities_pand_ipcw <- read_csv("../mpx-engage-params/data-processed/pand_ipcw_3cities.csv")
data_3cities_post_ipcw <- read_csv("../mpx-engage-params/data-processed/post_ipcw_3cities.csv")

# create single dataset with all time periods & cities
data_3cities <- bind_rows(
  data_3cities_pre_ipcw,
  data_3cities_pand_ipcw,
  data_3cities_post_ipcw
) %>%
  mutate(time_pt = factor(time_pt, 
                          levels = c("Pre-Pandemic", "Pandemic", "Post-Restrictions"))) %>%
  mutate(city = recode_factor(city, mtl = "Montreal", trt = "Toronto", van = "Vancouver"))

# create city marker
CITIES <- c("Montreal", "Toronto", "Vancouver")
TIMEPTS <- c("Pre-Pandemic", "Pandemic", "Post-Restrictions")

table(data_3cities$time_pt, 
      data_3cities$city, 
      useNA = "ifany")

data_3cities <- data_3cities %>% 
  mutate(
    data_pt = factor(paste(city, time_pt, sep = "-"),
                     levels = paste(rep(CITIES, each = 3), TIMEPTS, sep = "-"))
  )

CITIES_DATAPTS <- paste(
  rep(CITIES, each = 3), 
  rep(c("Pre-Pandemic", "Pandemic", "Post-Restrictions"), 
      times = 3), 
  sep = "-"
)

# save time by aggregating individuals
# before computing Pr(Y = y) for 0 to 300
nrow(data_3cities) # nb of individuals, could be the same person at different timepoints
data_aggrt <- data_3cities %>% 
  count(data_pt, age_grp, rel_status,
        bath_m, bath_d,
        groupsex_m, groupsex_d, 
        apps_partn_m, apps_partn_d, 
        sex_work_m, sex_work_d)
nrow(data_aggrt) # nb of individual prediction points

# Fit Bayesian models ----

## create dummy variables
# age; reference is 18-29
data_3cities <- make_ind_age(data_3cities)

# partnership status
table(data_3cities$data_pt, data_3cities$rel_status, useNA = "ifany")
table(data_3cities$reg_partn, data_3cities$rel_status, useNA = "ifany")
data_3cities <- make_ind_rel(data_3cities)

# check SPVs and groupsex
table(data_3cities$data_pt, data_3cities$bath_d, useNA = "ifany")
table(data_3cities$data_pt, data_3cities$groupsex_d, useNA = "ifany")

# check apps variable
# NOTE: apps variable not present in FU visit
table(data_3cities$data_pt, data_3cities$apps_partn_m, useNA = "ifany")

# check sex work variable
table(data_3cities$data_pt, data_3cities$sex_work_d, useNA = "ifany")

# hiv status 
table(data_3cities$data_pt, data_3cities$hiv_stat, useNA = "ifany")

# hiv status x age
data_3cities <- data_3cities %>% 
  mutate(
    hiv_age_30_39 = as.integer(age_30_39 == 1 & hiv_stat == 1),
    hiv_age_40_49 = as.integer(age_40_49 == 1 & hiv_stat == 1),
    hiv_age_50_59 = as.integer(age_50_59 == 1 & hiv_stat == 1),
    hiv_age_60_ = as.integer(age_60_ == 1 & hiv_stat == 1),
  )

### fit model
negbin_model <- stan_model(file = "./src-stan/regression_negbin_aggregate.stan",
                           model_name = "negbin_partn")

# variables to use (apps_partn_m and apps_partn_d are removed from FU since it was not asked)
vars_model_base <- c("age_30_39", "age_40_49", "age_50_59", "age_60_",
                     "rel_y_excl", "rel_y_open", "rel_y_uncl", 
                     "hiv_stat",
                     "hiv_age_30_39", "hiv_age_40_49", "hiv_age_50_59", "hiv_age_60_",
                     "bath_m", "bath_d",
                     "groupsex_m", "groupsex_d", 
                     "apps_partn_m", "apps_partn_d",
                     "sex_work_m", "sex_work_d")
vars_model_fu <- setdiff(vars_model_base, 
                         c("apps_partn_m", "apps_partn_d"))

fit_bayes_ls <- create_city_list(CITIES_DATAPTS)

# prepare data_frame for aggregate data points
data_x_aggrt <- data_3cities %>% 
  group_by(data_pt, across(all_of(vars_model_base))) %>% 
  summarize(nb = n(), 
            ipw_rds = sum(ipw_rds), 
            .groups = "drop") %>% 
  select(data_pt, nb, ipw_rds, all_of(vars_model_base))

## select city to standardize to
target_city <- "Montreal"
data_x_aggrt <- data_x_aggrt %>% filter(grepl(target_city, data_pt))

num_cores <- detectCores()
t0 <- Sys.time()

for(cur_city in CITIES_DATAPTS){
  # tracker
  if( grepl("-Pre-Pandemic", cur_city) ){
    print(gsub("-Pre-Pandemic", "", cur_city))
  }
  
  target_pop_stand <- gsub("Vancouver|Montreal|Toronto", target_city, cur_city)
  
  # choose baseline or follow-up variables
  if(grepl("-Pre-Pandemic", cur_city) ){
    vars_model <- vars_model_base
  } else {
    vars_model <- vars_model_fu
  }
  
  fit_bayes_ls[[cur_city]] <- sampling(
    negbin_model,
    data = list(y = filter(data_3cities, data_pt == cur_city)[[outcome_var]],
                
                # data on which model is fit
                x = data_3cities[data_3cities$data_pt == cur_city, vars_model],
                N = sum(data_3cities$data_pt == cur_city),
                
                # data to compute predictions
                x_aggr = data_x_aggrt[data_x_aggrt$data_pt == target_pop_stand, vars_model],
                N_aggr = sum(data_x_aggrt$data_pt == target_pop_stand),
                K = length(vars_model),
                
                # data for PMF and CDF
                x_end = 300,
                ipc_rds_w = data_x_aggrt$ipw_rds[data_x_aggrt$data_pt == target_pop_stand]),
    cores = num_cores,
    chains = 2, iter = 6000
  )
}
t1 <- Sys.time()
t1 - t0

# Inspect model results ----
# model convergence diagnostic
for(cur_city in CITIES_DATAPTS){
  cur_p_trace_plot <- traceplot(fit_bayes_ls[[cur_city]], pars = c("alpha", "beta", "phi"))
  ggsave(sprintf("%s/model-checks-p6m-all-%s-%s.png", fig_path, which(cur_city == CITIES_DATAPTS), cur_city),
         cur_p_trace_plot, device = "png",
         height = 14, width = 30, units = "cm", dpi = 320)
}
rm(cur_p_trace_plot)

# r hat
ess_ls <- create_city_list(CITIES_DATAPTS)
for(cur_city in CITIES_DATAPTS){
  print(cur_city)
  # get model
  cur_model <- summary(fit_bayes_ls[[cur_city]])$summary
  
  # show only intercept and regression coefficients, ignore y_hat and y_pred
  row_param_names <- rownames(cur_model)
  row_param_names <- grep("alpha|beta|phi", row_param_names, value = T)
  
  # output
  print(round(cur_model[row_param_names, ], 3))
  
  # save
  ess_ls[[cur_city]] <- cur_model[row_param_names, c("mean", "se_mean", "2.5%", "50%", "97.5%", "n_eff")]
}
rm(cur_model, row_param_names)

## effective sample size
# Check that all coefficients have >1,000 sample size
ess_ls_tbl <- vector("list", length(ess_ls))
for(i in 1:length(ess_ls)){
  ess_ls_tbl[[i]] <- as_tibble(ess_ls[[i]], rownames = "coeff")
  ess_ls_tbl[[i]] <- mutate(ess_ls_tbl[[i]], city.time = names(ess_ls)[i], .before = 1)
}
ess_tbl <- bind_rows(ess_ls_tbl)
rm(ess_ls_tbl)

# save summary by city & time period
summarize_ess(ess_tbl, beta_only = F)
summarize_ess(ess_tbl, beta_only = T)

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
