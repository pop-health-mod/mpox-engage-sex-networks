# Library and data----
library(tidyverse)
library(data.table)
library(rstan)

source("./src/utils_helper.R")
source("./src/utils_regression.R")
theme_set(theme_bw())

## main analyses
outcome_var <- "nb_part_ttl"

## sensitivity analysis with number of *anal* partners as outcome
# outcome_var <- "nb_part_anal"

## sensitivity analysis with zero-inflated negative binomial
DO_ZINF <- FALSE

## define paths & prefixes based on the analysis being done
fig_path <- case_when(outcome_var == "nb_part_ttl" ~ "./fig/results-checks",
                      outcome_var == "nb_part_anal" ~ "./fig/results-checks-anal",
                      DO_ZINF ~ "./fig/results-checks-zinf")

stan_model_path <- ifelse(DO_ZINF, "./src-stan/regression_negbin_zinf_aggregate.stan",
                          "./src-stan/regression_negbin_aggregate.stan")

out_distr_path <- case_when(outcome_var == "nb_part_ttl" ~ "./out/fitted-distributions",
                            outcome_var == "nb_part_anal" ~ "./out/fitted-distr-sens-anal",
                            DO_ZINF ~ "./out/fitted-distr-sens-zinf")
out_distr_pref <- case_when(outcome_var == "nb_part_ttl" ~ "",
                            outcome_var == "nb_part_anal" ~ "-anal",
                            DO_ZINF ~ "-zinf")

## load data
data_3cities_pre_ipcw <- read_csv("../mpx-engage-params/data-3cities-feb-2023/pre_ipcw_3cities.csv")
data_3cities_pand_ipcw <- read_csv("../mpx-engage-params/data-3cities-feb-2023/pand_ipcw_3cities.csv")
data_3cities_post_ipcw <- read_csv("../mpx-engage-params/data-3cities-feb-2023/post_ipcw_3cities.csv")

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
negbin_model <- stan_model(file = stan_model_path,
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

num_cores <- parallel::detectCores()
t0 <- Sys.time()

# doing the PMF computations in stan, we get only a [-0.000019  0.000112] difference in estimates
# for Montreal-Pre-Pandemic
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
                K = length(vars_model),
                
                # data for PMF and CDF
                x_end = 300,
                ipc_rds_w = data_x_aggrt$ipw_rds[data_x_aggrt$data_pt == cur_city]),
    cores = num_cores,
    chains = 2, iter = 4000
  )
}

t1 <- Sys.time()
t1 - t0 # ~9 minutes (regression and PMF)

# save model fits
# saveRDS(fit_bayes_ls, sprintf("./out/text_stan_model_fit%s.rds", out_distr_pref))

## Inspect model convergence diagnostic ----
# convergence of model chains (traceplots)
for(cur_city in CITIES_DATAPTS){
  # MCMCtrace(fit_bayes_ls[[cur_city]], 
  #           params = c("alpha", "beta", "phi"),
  #           filename = paste0("traceplot-", cur_city),
  #           wd = sprintf("%s/model-checks-p6m-all", fig_path),
  #           open_pdf = F)
  cur_p_trace_plot <- traceplot(fit_bayes_ls[[cur_city]], pars = c("alpha", "beta", "phi"))
  ggsave(sprintf("%s/model-checks-p6m-all-%s-%s.png", fig_path, which(cur_city == CITIES_DATAPTS), cur_city),
         cur_p_trace_plot, device = "png",
         height = 14, width = 30, units = "cm", dpi = 320)
}
rm(cur_p_trace_plot)

##### can stop here for sensitivity analyses

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

if(outcome_var == "nb_part_ttl" & !DO_ZINF){
  write_csv(ess_tbl, "./out/stan_model_fit_ess.csv")
}

## Regression coefficients (RR) ----
## extract coefficients (only for post-restrictions period)
coeff_ls <- vector("list", 3)
names(coeff_ls) <- paste(CITIES, "-Post-Restrictions", sep = "")

for(cur_city in names(coeff_ls)){
  coeff_ls[[cur_city]] <- as_tibble(
    summary(fit_bayes_ls[[cur_city]], pars = c("alpha", "beta", "phi"))$summary,
    rownames = "term"
  )
  
  coeff_ls[[cur_city]] <- coeff_ls[[cur_city]] %>% 
    mutate(data_pt = cur_city, .before = 1)
  }

tbl_coeff <- bind_rows(coeff_ls)

# add variable names
vars_model_fu_full <- c("Age 30-39", "Age 40-49", "Age 50-59", "Age ≥60", 
                        "Exclusive Relationship", "Open Relationship", "Unclear Relationship", 
                        "HIV Status",
                        "HIV x Age 30-39", "HIV x Age 40-49", "HIV x Age 50-59", "HIV x Age ≥60",
                        "Bathhouse", "Bathhouse Missing",
                        "Groupsex", "Groupsex Missing",
                        # "Dating Apps", "Dating Apps Missing",
                        "Transactional Sex", "Transactional Sex Missing")

tbl_coeff <- add_column(tbl_coeff,
                        name = rep(c("intercept", vars_model_fu_full, "inverse_overdisp"), times = 3),
                        .after = "term")


# table for posteriors (post-restriction only)
range(exp(tbl_coeff$`97.5%`))

tbl_coeff <- tbl_coeff %>% 
  mutate(coeff = factor(name, levels = unique(name))) %>%
  mutate(city = gsub("-Post-Restrictions", "", data_pt),
         time_pt = gsub("Montreal-|Toronto-|Vancouver-", "", data_pt),
         .before = 1) %>%
  mutate(time_pt = factor(time_pt, 
                          levels = c("Pre-Pandemic", "Pandemic", "Post-Restrictions")))

tbl_coeff_post <- tbl_coeff %>% filter(time_pt == "Post-Restrictions")

coeff_post_tbl <- tbl_coeff_post %>% 
  select(mean, city, time_pt, coeff, mean, se_mean, `2.5%`, `97.5%`) %>%
  mutate(SE = formatC(signif(se_mean,digits = 2), 
                      digits = 2,
                      format = "fg",
                      flag = "#"),
         Mean = paste0(round(exp(mean), 2), " (", round(exp(`2.5%`), 2), ", ", round(exp(`97.5%`), 2), ")")) %>%
  select(city, coeff, Mean, SE)

# turn into one set of columns per city
coeff_post_tbl <- coeff_post_tbl %>% 
  mutate(city_code = case_when(city == "Montreal" ~ "mtl",
                               city == "Toronto" ~ "trt",
                               city == "Vancouver" ~ "van")) %>% 
  select(-city) %>% 
  pivot_wider(names_from = "city_code", values_from = c("Mean", "SE"))

# reorder columns
coeff_post_tbl <- coeff_post_tbl %>% select(coeff, ends_with("mtl"), ends_with("trt"), ends_with("van"))

# only save coefficients of main analyses
if(outcome_var == "nb_part_ttl" & !DO_ZINF){
  write.csv(coeff_post_tbl, "./out/manuscript-tables/table_S4_coef_post.csv")
}
