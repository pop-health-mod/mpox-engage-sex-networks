# library and data----
library(tidyverse)
library(data.table)
library(parallel)
library(doParallel)
library(foreach)
library(gridExtra)
library(MCMCvis)
library(rstan)

source("./src/utils_regression_3cities.R")
source("./src/utils_helper.R")
source("./src/plot.R")
set.seed(11)
theme_set(theme_bw())

## which variable to use as outcome
outcome_var <- "nb_part_ttl" 
file_suff <- case_when(outcome_var == "nb_part_ttl" ~ "p6m_all",
                       outcome_var == "nb_part_anal" ~ "p6m_anal")

# restriction analysis
fig_path <- "./figures-3cities/negbin-res"
data_3cities <- read_csv("./data-3cities-feb-2023/restriction/res_fu.csv")
data_3cities <- data_3cities %>%
  mutate(ipw_rds = wt_rds_norm) %>%
  mutate(time_pt = factor(time_pt, levels = c("Pre-Pandemic", "Pandemic", "Post-Restriction"))) %>%
  mutate(city = recode_factor(city, mtl = "Montreal", trt = "Toronto", van = "Vancouver"))

# create city marker
CITIES <- c("Montreal", "Toronto", "Vancouver")
TIMEPTS <- c("Pre-Pandemic", "Pandemic", "Post-Restriction")

table(data_3cities$time_pt, 
      data_3cities$city, 
      useNA = "ifany")

data_3cities <- mutate(data_3cities, 
                       data_pt = paste(city, time_pt, sep = "-"))

CITIES_DATAPTS <- paste(
  rep(CITIES, each = 3), 
  rep(c("Pre-Pandemic", "Pandemic", "Post-Restriction"), 
      times = 3), 
  sep = "-"
)

# could save time by aggregating individuals
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

### fit model
negbin_model <- stan_model(file = "./src-stan/regression_negbin-1_aggregate.stan",
                           model_name = "negbin_partn")

# variables to use (apps_partn_m and apps_partn_d are removed from FU since it was not asked)
vars_model_base <- c("age_30_39", "age_40_49", "age_50_59", "age_60_",
                     "rel_y_excl", "rel_y_open", "rel_y_uncl", 
                     "hiv_stat",
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
  dplyr::select(data_pt, nb, ipw_rds, all_of(vars_model_base))


num_cores <- detectCores()
t0 <- Sys.time()

for(cur_city in CITIES_DATAPTS){
  # tracker
  if( grepl("-Pre-Pandemic", cur_city) ){
    print(gsub("-Pre-Pandemic", "", cur_city))
  }
  
  # choose Pre-Pandemic or follow-up variables
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
                x_aggr = data_x_aggrt[data_x_aggrt$data_pt == cur_city, vars_model],
                N_aggr = sum(data_x_aggrt$data_pt == cur_city),
                K = length(vars_model)),
    cores = num_cores,
    chains = 2, iter = 4000
  )
}

t1 <- Sys.time()
t1 - t0 

# Inspect model results ----
# model convergence diagnostic
for(cur_city in CITIES_DATAPTS){
  MCMCtrace(fit_bayes_ls[[cur_city]], 
            params = c("alpha", "beta", "phi"),
            filename = paste0("traceplot-", cur_city),
            wd = sprintf("./figures-3cities/negbin-res/model-check%s", ifelse(file_suff == "p6m_all", "", "-anal-partn")),
            open_pdf = F)
}

# get posterior predictive distribution ----
# get distribution of expected predictions
# rows: iterations, columns: each individual's prediction (mean of the NB)
data_ypred <- create_city_list(CITIES_DATAPTS)
phi_all <- create_city_list(CITIES_DATAPTS)

for(cur_city in CITIES_DATAPTS){
  data_ypred[[cur_city]] <- extract(fit_bayes_ls[[cur_city]], pars = "y_pred")$y_pred
  phi_all[[cur_city]] <- extract(fit_bayes_ls[[cur_city]], pars = "phi")$phi
}

density_list_by_iter <- create_city_list(CITIES_DATAPTS)

# start a cluster
numCores <- parallel::detectCores()
numCores
numCores <- numCores - 1
cl <- parallel::makeCluster(numCores, type = "PSOCK")

# register
doParallel::registerDoParallel(cl = cl) 
# foreach::getDoParRegistered()

t0 <- Sys.time()
for(cur_city in CITIES_DATAPTS){
  
  print(cur_city)
  
  density_list_by_iter[[cur_city]] <- foreach(i = 1:nrow(data_ypred[[cur_city]])) %dopar% {
    compute_dens_negbin.single(
      data_ypred[[cur_city]][i, ],
      range_x = 0:300,
      phi = phi_all[[cur_city]][i]
    )
  }
}
t1 <- Sys.time()
t1 - t0 

# close cluster
stopCluster(cl)

#   list                          list      matrix
# density_list_by_iter[[city]][[iteration]][0:300, individuals]

# make lists to hold data for each city
tmp_density <- create_city_list(CITIES_DATAPTS) # 4000 iterations per city-data_pt 
dens_wt_by_city <- create_city_list(CITIES_DATAPTS)

t0 <- Sys.time()

for(cur_city in CITIES_DATAPTS){
  cat(sprintf("%s\n", cur_city))
  # compute density in the city and visit
  tmp_density[[cur_city]] <- compute_pmf(
    density_list_by_iter[[cur_city]],
    filter(data_x_aggrt, data_pt == cur_city)$ipw_rds,
    nb_aggregate = filter(data_x_aggrt, data_pt == cur_city)$nb
  )
  
  # # all densities should sum up to 1 (VERIFIED)
  # tmp_density[[cur_city]] %>% group_by(iter) %>% summarize(dens_ttl = sum(dens_wt)) %>%
  #   pull(dens_ttl) %>% unique() %>% print()
  
  dens_wt_by_city[[cur_city]] <- summarize_density(tmp_density[[cur_city]], "dens_wt")
}

t1 <- Sys.time()
t1 - t0 

# collapse into single dataset --------------------------------------------
# collapse each city

for(cur_city in CITIES_DATAPTS){
  dens_wt_by_city[[cur_city]]$data_pt <- cur_city
}


dens_wt_by_city <- bind_rows(dens_wt_by_city)

# reorder to put city and relationship first
dens_wt_by_city <- dens_wt_by_city[, .(data_pt, y_pred, mean, mdn, cr.i_low, cr.i_upp)]

# verification of posterior distributions and densities ----
# verify results by looking at the mean number of partners
data_mean_nb_partn <- dens_wt_by_city %>%
  group_by(data_pt) %>%
  summarize(mean_wt = sum(y_pred * mean),
            cr.i_low = sum(y_pred * cr.i_low),
            cr.i_upp = sum(y_pred * cr.i_upp),
            .groups = "drop") %>% 
  mutate(type = "neg. bin.")

data_mean_nb_partn$data_pt <- factor(data_mean_nb_partn$data_pt, levels = CITIES_DATAPTS)

# verify that weights sum up to 1
dens_wt_by_city %>% 
  group_by(data_pt) %>% 
  summarize(dens_ttl = sum(mean), .groups = "drop")

# output tables (pmf and cdf) ----
# full fitted distribution pmf
write.csv(dens_wt_by_city, sprintf("./output-3cities/negbin-res/full_fitted_pmf_wt_%s.csv", file_suff),
          row.names = F)

# compute fitted CDF in each datapoint (similar procedure as for pmf)
# make lists to hold data for each city
cdf_wt_by_city <- create_city_list(CITIES_DATAPTS)

t0 <- Sys.time()
for(cur_city in CITIES_DATAPTS){
  cat(sprintf("%s\n", cur_city))
  # compute CDF in the city/timepoint
  tmp_density <- compute_cdf(
    density_list_by_iter[[cur_city]],
    filter(data_x_aggrt, data_pt == cur_city)$ipw_rds,
    nb_aggregate = filter(data_x_aggrt, data_pt == cur_city)$nb
  )
  
  # compute weighted densities
  cdf_wt_by_city[[cur_city]] <- summarize_density(tmp_density, "cdf_wt")
}
t1 <- Sys.time()
t1 - t0 

# collapse into single dataset
for(cur_city in CITIES_DATAPTS){
  cdf_wt_by_city[[cur_city]]$data_pt <- cur_city
}

# collapse into a single dataframe
cdf_wt_by_city <- bind_rows(cdf_wt_by_city)

# reorder to put city and relationship first
cdf_wt_by_city <- cdf_wt_by_city[, .(data_pt, y_pred, mean, mdn, cr.i_low, cr.i_upp)]


cdf_wt_by_city <- cdf_wt_by_city %>% 
  mutate(city = gsub("-Pre-Pandemic|-Pandemic|-Post-Restriction", "", data_pt),
         time_pt = gsub("Montreal-|Toronto-|Vancouver-", "", data_pt),
         .before = 1)

# factorize time points to order correctly
cdf_wt_by_city$time_pt <- factor(cdf_wt_by_city$time_pt, 
                                 levels = c("Pre-Pandemic", "Pandemic", "Post-Restriction"))

# save full fitted distribution
write.csv(cdf_wt_by_city, sprintf("./output-3cities/negbin-res/full_fitted_cdf_wt_%s.csv", file_suff),
          row.names = F)
