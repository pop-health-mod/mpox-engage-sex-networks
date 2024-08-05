# Library and data----
library(tidyverse)
library(data.table)
library(rstan)

source("./src/utils_helper.R")
source("./src/utils_regression.R")
theme_set(theme_bw())

## define paths & output prefixes
fig_path <- "./misc-data-proc/fig-cm"

stan_model_path <- "./src-stan/regression_negbin_aggregate_age_hiv.stan"

out_distr_path <- "./misc-data-proc/outputs-cm"

## load data
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
                          levels = c("Pre-Pandemic", "Pandemic", "Post-Restrictions")))

data_van <- filter(data_3cities, city == "van")

# recode age groups
unique(data_van$age_grp_cat)
data_van <- data_van %>% 
  mutate(
    age_grp_cat = case_when(age <  30 ~ "<30",
                            age <  50 ~ "30-49",
                            age >= 50 ~ "50+"),
    age_grp_cat = factor(age_grp_cat,
                         levels = c("<30", "30-49", "50+"))
  )

table(data_van$age_grp_cat, data_van$age_grp)

# create city marker
TIMEPTS <- c("Pre-Pandemic", "Pandemic", "Post-Restrictions")

AGES        <- paste0("age",sort(unique(data_van$age_grp)))
AGES_AGGR   <- paste0("age",sort(unique(data_van$age_grp_cat)))
HIV         <- paste0("hiv", unique(data_3cities$hiv_stat))
AGEHIV      <- paste(rep(AGES, each = 2), HIV, sep = "_")
AGEHIV_AGGR <- paste(rep(AGES_AGGR, each = 2), HIV, sep = "_")

table(data_van$time_pt, data_van$city, useNA = "ifany")

data_van <- data_van %>% 
  mutate(
    data_pt = factor(paste(city, time_pt, sep = "-"),
                     levels = paste(rep("van", each = 3), TIMEPTS, sep = "-")),
    age_hiv = factor(sprintf("age%s_hiv%s", age_grp, hiv_stat),
                     levels = AGEHIV)
  )

CITIES_DATAPTS <- paste(rep("van", each = 3),
                        c("Pre-Pandemic", "Pandemic", "Post-Restrictions"),
                        sep = "-")

# save time in the PMF step by aggregating individuals
# before computing Pr(Y = y) for 0 to 300
nrow(data_van) # nb of individuals, could be the same person at different timepoints
data_aggrt <- data_van %>% 
  count(data_pt, age_grp, rel_status,
        bath_m, bath_d,
        groupsex_m, groupsex_d, 
        apps_partn_m, apps_partn_d, 
        sex_work_m, sex_work_d)

data_aggrt_ah <- data_van %>% 
  group_by(age_hiv) %>% 
  count(data_pt, rel_status,
        bath_m, bath_d,
        groupsex_m, groupsex_d, 
        apps_partn_m, apps_partn_d, 
        sex_work_m, sex_work_d) %>% 
  ungroup()

# nb of individual prediction points (city-wide and for subgroups)
nrow(data_aggrt)
count(data_aggrt_ah, age_hiv)

# Fit Bayesian model ----

## create dummy variables
# age; reference is 18-29
data_van <- make_ind_age(data_van)

# partnership status
table(data_van$data_pt, data_van$rel_status, useNA = "ifany")
table(data_van$reg_partn, data_van$rel_status, useNA = "ifany")
data_van <- make_ind_rel(data_van)

# check SPVs and groupsex
table(data_van$data_pt, data_van$bath_d, useNA = "ifany")
table(data_van$data_pt, data_van$groupsex_d, useNA = "ifany")

# check apps variable
# NOTE: apps variable not present in FU visit
table(data_van$data_pt, data_van$apps_partn_m, useNA = "ifany")

# check sex work variable
table(data_van$data_pt, data_van$sex_work_d, useNA = "ifany")

# hiv status 
table(data_van$data_pt, data_van$hiv_stat, useNA = "ifany")

# hiv status x age
data_van <- data_van %>% 
  mutate(
    hiv_age_30_39 = as.integer(age_30_39 == 1 & hiv_stat == 1),
    hiv_age_40_49 = as.integer(age_40_49 == 1 & hiv_stat == 1),
    hiv_age_50_59 = as.integer(age_50_59 == 1 & hiv_stat == 1),
    hiv_age_60_ = as.integer(age_60_ == 1 & hiv_stat == 1),
  )

### fit model
negbin_model <- stan_model(file = stan_model_path, model_name = "negbin_partn")

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
data_x_aggrt <- data_van %>% 
  group_by(data_pt, across(all_of(vars_model_base))) %>% 
  summarize(nb = n(), 
            ipw_rds = sum(ipw_rds), 
            .groups = "drop") %>% 
  select(data_pt, nb, ipw_rds, all_of(vars_model_base))

# prepare data_frame for aggregate data points
df_x_aggrt_ah <- split(data_van, data_van$age_hiv)
data_x_aggrt_ah <- list()

for(cur_city in CITIES_DATAPTS){
  data_x_aggrt_ah_cur_city <- list()
  for(ah in AGEHIV){
    cur_city_index <- df_x_aggrt_ah[[ah]]$data_pt == cur_city
    cur_city_data <- df_x_aggrt_ah[[ah]][cur_city_index, ]
    
    # count how many individuals and their total weight in each covariate combination
    ## for the predictions
    data_x_aggrt_ah_cur_city[[ah]] <- cur_city_data %>%
      group_by(across(all_of(vars_model_base))) %>% 
      summarize(nb = n(), 
                ipw_rds = sum(ipw_rds), 
                .groups = "drop") %>% 
      select(nb, ipw_rds, all_of(vars_model_base))
  }
  
  cur_city_max_combo <- max(unname(unlist(lapply(data_x_aggrt_ah_cur_city, nrow))))
  
  data_x_aggrt_ah[[cur_city]] <- array(data = 0,
                                       dim = c(length(AGEHIV), cur_city_max_combo, ncol(data_x_aggrt_ah_cur_city[[ah]])),
                                       dimnames = list(AGEHIV, 1:cur_city_max_combo, colnames(data_x_aggrt_ah_cur_city[[ah]])))
  
  # create the matrices
  for(ah in AGEHIV){
    n_row <- nrow(data_x_aggrt_ah_cur_city[[ah]])
    n_col <- ncol(data_x_aggrt_ah_cur_city[[ah]])
    for(r in 1:n_row){
      for(c in 1:n_col){
        data_x_aggrt_ah[[cur_city]][ah, r, c] <- unname(unlist(data_x_aggrt_ah_cur_city[[ah]][r, c]))
      }
    }
  }
}

num_cores <- parallel::detectCores()
t0 <- Sys.time()

# doing the PMF computations in stan, we get only a [-0.000019  0.000112] difference in estimates
# for Montreal-Pre-Pandemic
for(cur_city in CITIES_DATAPTS){
  # tracker
  print(gsub("van-", "", cur_city))
  
  # choose Pre-Pandemic or follow-up variables
  if(grepl("-Pre-Pandemic", cur_city)){
    vars_model <- vars_model_base
  } else {
    vars_model <- vars_model_fu
  }
  
  fit_bayes_ls[[cur_city]] <- sampling(
    negbin_model,
    data = list(y = filter(data_van, data_pt == cur_city)$nb_part_ttl,
                
                # data on which model is fit
                x = data_van[data_van$data_pt == cur_city, vars_model],
                N = sum(data_van$data_pt == cur_city),
                
                # data to compute city-wide predictions (including prediction targets and weights)
                x_aggr = data_x_aggrt[data_x_aggrt$data_pt == cur_city, vars_model],
                N_aggr = sum(data_x_aggrt$data_pt == cur_city),
                K = length(vars_model),
                x_end = 300,
                ipc_rds_w = data_x_aggrt$ipw_rds[data_x_aggrt$data_pt == cur_city],
                
                # data to compute age x HIV predictions
                N_ah = length(AGEHIV),
                x_aggr_ah = data_x_aggrt_ah[[cur_city]][, , vars_model],
                N_aggr_ah = dim(data_x_aggrt_ah[[cur_city]])[2],
                ipc_rds_w_ah = data_x_aggrt_ah[[cur_city]][, , "ipw_rds"]),
    cores = num_cores,
    chains = 2, iter = 6000
  )
}

t1 <- Sys.time()
t1 - t0 # ~4.1 minutes (regression and PMF)

# save model fits
# saveRDS(fit_bayes_ls, "./misc-data-proc/outputs-cm/stan_model_fit_cm.rds")

## Inspect model convergence diagnostic ----
# convergence of model chains (traceplots)
for(cur_city in CITIES_DATAPTS){
  cur_p_trace_plot <- traceplot(fit_bayes_ls[[cur_city]], pars = c("alpha", "beta", "phi"))
  ggsave(sprintf("%s/model-checks-p6m-all-van-%s.png", fig_path, which(cur_city == CITIES_DATAPTS)),
         cur_p_trace_plot, device = "png",
         height = 14, width = 30, units = "cm", dpi = 320)
}
rm(cur_p_trace_plot)

# r hat and effective sample size
ess_ls <- create_city_list(CITIES_DATAPTS)
for(cur_city in CITIES_DATAPTS){
  print(cur_city)
  # get model
  cur_model <- summary(fit_bayes_ls[[cur_city]])$summary
  
  # show only intercept and regression coefficients, ignore y_hat and y_pred
  row_param_names <- rownames(cur_model)
  row_param_names <- grep("alpha|beta|phi|shape|zi", row_param_names, value = T)
  
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

# City-wide distribution ----
## Probability mass function ----
## density and PMF computations already performed in Stan
### Collapse PMF into single dataset ----
pmf_iter <- create_city_list(CITIES_DATAPTS)
pmf_wt_cty <- create_city_list(CITIES_DATAPTS)

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
  pmf_wt_cty[[cur_city]] <- tibble(data_pt = cur_city,
                                   y_pred = 0:300,
                                   mean = colSums(pmf_tmp) / nrow(pmf_tmp),
                                   cr.i_low = cred_l,
                                   cr.i_upp = cred_u)
  
  rm(cred_l, cred_u)
}
rm(pmf_tmp)

# collapse into single dataframe
pmf_wt_cty <- bind_rows(pmf_wt_cty)
pmf_wt_cty$data_pt <- factor(pmf_wt_cty$data_pt, levels = CITIES_DATAPTS)

pmf_wt_cty <- rename(pmf_wt_cty, time_prd = data_pt)
pmf_wt_cty <- mutate(pmf_wt_cty, time_prd = gsub("van-", "", time_prd))

### Verify PMF posterior distributions ----
# verify results by looking at the mean number of partners
data_mean_nb_partn <- pmf_wt_cty %>%
  mutate(time_prd = factor(time_prd, levels = levels(data_van$time_pt))) %>% 
  group_by(time_prd) %>%
  summarize(mean_wt = sum(y_pred * mean),
            cr.i_low = sum(y_pred * cr.i_low),
            cr.i_upp = sum(y_pred * cr.i_upp),
            .groups = "drop")

data_mean_nb_partn
data_van %>% 
  group_by(time_pt) %>% 
  summarize(mean_nb_partn_p6m = sum(nb_part_ttl * ipw_rds) / sum(ipw_rds))

# verify that weights sum up to 1
pmf_wt_cty %>% 
  group_by(time_prd) %>% 
  summarize(dens_ttl = sum(mean), .groups = "drop")

## Cumulative density function ----
# compute fitted CDF in each datapoint (similar procedure as for pmf)
# make lists to hold data for each city
cdf_wt_cty <- create_city_list(CITIES_DATAPTS)
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
  cdf_wt_cty[[cur_city]] <- tibble(data_pt = cur_city,
                                   y_pred = 0:300,
                                   mean = colSums(cdf_tmp) / nrow(cdf_tmp),
                                   cr.i_low = cred_l,
                                   cr.i_upp = cred_u)
  
  rm(cred_l, cred_u)
}
rm(pmf_tmp, cdf_tmp)

# collapse into a single dataframe
cdf_wt_cty <- bind_rows(cdf_wt_cty)

cdf_wt_cty <- cdf_wt_cty %>% 
  mutate(time_prd = gsub("van-", "", data_pt), .before = 1) %>% 
  select(-data_pt)

# factorize time points to order correctly
cdf_wt_cty$time_prd <- factor(cdf_wt_cty$time_prd,
                              levels = c("Pre-Pandemic", "Pandemic", "Post-Restrictions"))

## Output tables (PMF and CDF) ----
# save iterations (only for main analyses)
saveRDS(pmf_iter, "./misc-data-proc/outputs-cm/pmf_stan_iterations_van.rds")

# save full fitted pmf
write.csv(pmf_wt_cty, "./misc-data-proc/outputs-cm/engage_van_fitted_pmf_nb_partn_p6m.csv",
          row.names = F)

# save full fitted cdf
write.csv(cdf_wt_cty, "./misc-data-proc/outputs-cm/engage_van_fitted_cdf_nb_partn_p6m.csv",
          row.names = F)

# Age x HIV group distributions ----
## Probability mass function ----
## density and PMF computations already performed in Stan

### Collapse PMF into single dataset ----
pmf_iter_ah <- create_city_list(CITIES_DATAPTS)
pmf_wt_by_city_ah <- create_city_list(CITIES_DATAPTS)

# get correspondence between the final groups and the ones used in the regression
ls_ah_corresp <- vector("list", length(AGEHIV_AGGR))
names(ls_ah_corresp) <- AGEHIV_AGGR
{
  ls_ah_corresp[[1]] <- AGEHIV[1]
  ls_ah_corresp[[2]] <- AGEHIV[2]
  ls_ah_corresp[[3]] <- AGEHIV[c(3, 5)]
  ls_ah_corresp[[4]] <- AGEHIV[c(4, 6)]
  ls_ah_corresp[[5]] <- AGEHIV[c(7, 9)]
  ls_ah_corresp[[6]] <- AGEHIV[c(8, 10)]
}
ls_ah_corresp

for(cur_city in CITIES_DATAPTS){
  for(ah in AGEHIV_AGGR){
    # retrieve correspondence
    ah_tgt <- ls_ah_corresp[[ah]]
      
    if(length(ah_tgt) == 1){
      index_ah <- which(AGEHIV == ah_tgt)
      
      # extract PMF iterations from stan
      pmf_tmp_ah <- extract(fit_bayes_ls[[cur_city]], pars = "pmf_ah")$pmf_ah[, index_ah, ]
      pmf_iter_ah[[cur_city]][[ah]] <- pmf_tmp_ah
      
      # get credible intervals
      cred_l_ah <- vector("double", ncol(pmf_tmp_ah))
      cred_u_ah <- vector("double", ncol(pmf_tmp_ah))
      
      for(i in 1:ncol(pmf_tmp_ah)){
        cred_l_ah[i] <- quantile(pmf_tmp_ah[, i], .025)
        cred_u_ah[i] <- quantile(pmf_tmp_ah[, i], .975)
      }
      # get mean pt estimate
      mean_pmf <- colSums(pmf_tmp_ah) / nrow(pmf_tmp_ah)
      
    } else {
      index_ah <- which(AGEHIV %in% ah_tgt)
      
      # extract PMF iterations from stan
      pmf_tmp_ah <- extract(fit_bayes_ls[[cur_city]], pars = "pmf_ah")$pmf_ah[, index_ah, ]
      pmf_tmp_ah_final <- matrix(0,
                                 nrow = nrow(pmf_tmp_ah[, 1,]),
                                 ncol = ncol(pmf_tmp_ah[, 1,]))
      
      # get the total weights from each subgroup
      mtx_ipw_rds <- data_x_aggrt_ah[[cur_city]][AGEHIV[index_ah], , "ipw_rds"]
      vec_ipw_rds <- rowSums(mtx_ipw_rds)
      
      # add up the weighted PMF from each subgroup, and then divide by the total weight
      for(i in 1:length(ah_tgt)){
        pmf_tmp_ah_final <- pmf_tmp_ah_final + (pmf_tmp_ah[, i, ] * vec_ipw_rds[i])
      }
      pmf_tmp_ah_final <- pmf_tmp_ah_final / sum(vec_ipw_rds)
      
      # back to storing the end result and computing the credible intervals
      pmf_iter_ah[[cur_city]][[ah]] <- pmf_tmp_ah_final
      
      # get credible intervals
      cred_l_ah <- vector("double", ncol(pmf_tmp_ah_final))
      cred_u_ah <- vector("double", ncol(pmf_tmp_ah_final))
      
      for(i in 1:ncol(pmf_tmp_ah_final)){
        cred_l_ah[i] <- quantile(pmf_tmp_ah_final[, i], .025)
        cred_u_ah[i] <- quantile(pmf_tmp_ah_final[, i], .975)
      }
      
      # get mean pt estimate
      mean_pmf <- colSums(pmf_tmp_ah_final) / nrow(pmf_tmp_ah_final)
      
      rm(mtx_ipw_rds, vec_ipw_rds)
    }
    
    # create tibble with each city and time period
    pmf_wt_by_city_ah[[cur_city]][[ah]] <- tibble(data_pt = cur_city,
                                                  age_hiv = ah,
                                                  y_pred = 0:300,
                                                  mean = mean_pmf,
                                                  cr.i_low = cred_l_ah,
                                                  cr.i_upp = cred_u_ah) %>%
      mutate(age_hiv = gsub("age|hiv", "", age_hiv),
             age = str_split(age_hiv, "_", simplify = T)[, 1],
             hiv = str_split(age_hiv, "_", simplify = T)[, 2],
             .after = 1) %>% 
      select(-age_hiv)
    
    rm(cred_l_ah, cred_u_ah)
  }
  
  pmf_wt_by_city_ah[[cur_city]] <- bind_rows(pmf_wt_by_city_ah[[cur_city]])
}

rm(pmf_tmp_ah, pmf_tmp_ah_final)

# collapse into single dataframe
pmf_wt_by_city_ah <- bind_rows(pmf_wt_by_city_ah)
pmf_wt_by_city_ah$data_pt <- factor(pmf_wt_by_city_ah$data_pt, levels = CITIES_DATAPTS)

### CONTINUE HERE

### Verify PMF posterior distributions ----
# verify results by looking at the mean number of partners
group_var <- c("age_cats", "hiv_cats", "data_pt")

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
            .groups = "drop") %>% 
  View()

# save full fitted pmf
write.csv(pmf_wt_by_city,
          sprintf("%s/pmf_weighted_%s.csv", out_distr_path, out_distr_pref),
          row.names = F)





