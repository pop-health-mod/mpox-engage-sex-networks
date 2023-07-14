# library and data----
library(tidyverse)
library(data.table)
library(parallel)
library(doParallel)
library(foreach)
library(gridExtra)
library(MCMCvis)
library(rstan)
library(ggtext)

source("./src/utils_regression_3cities.R")
source("./src/utils_helper.R")
source("./src/plot.R")
set.seed(111)
theme_set(theme_bw())

## which variable to use as outcome
outcome_var <- "nb_part_ttl" 
file_suff <- case_when(outcome_var == "nb_part_ttl" ~ "p6m_all",
                       outcome_var == "nb_part_anal" ~ "p6m_anal")
fig_path <- "./figures-3cities/negbin-ipcw"

data_3cities_pre_ipcw <- read_csv("./data-3cities-feb-2023/pre_ipcw_3cities.csv")
data_3cities_pand_ipcw <- read_csv("./data-3cities-feb-2023/pand_ipcw_3cities.csv")
data_3cities_post_ipcw <- read_csv("./data-3cities-feb-2023/post_ipcw_3cities.csv")

data_3cities <- data_3cities_pre_ipcw %>% 
  bind_rows(data_3cities_pand_ipcw,data_3cities_post_ipcw) %>%
  mutate(time_pt = factor(time_pt, 
                          levels = c("Pre-Pandemic", "Pandemic", "Post-Restriction"))) %>%
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

## Inspect model results ----
# model convergence diagnostic
for(cur_city in CITIES_DATAPTS){
  MCMCtrace(fit_bayes_ls[[cur_city]], 
            params = c("alpha", "beta", "phi"),
            filename = paste0("traceplot-", cur_city),
            wd = sprintf("./figures-3cities/negbin-ipcw/model-check%s", ifelse(file_suff == "p6m_all", "", "-anal-partn")),
            open_pdf = F)
}

## coefficients ------------------------------------------------------------

coeff_ls <- create_city_list(CITIES_DATAPTS)

for(cur_city in CITIES_DATAPTS){
  coeff_ls[[cur_city]] <- as_tibble(
    summary(fit_bayes_ls[[cur_city]], pars = c("alpha", "beta", "phi"))$summary,
    rownames = "term"
  )
  
  coeff_ls[[cur_city]] <- coeff_ls[[cur_city]] %>% 
    mutate(data_pt = cur_city, .before = 1)
  }

tbl_coeff <- bind_rows(coeff_ls)

# add variable names
vars_model_base_full <- c("Age 30-39", "Age 40-49", "Age 50-59", "Age ≥60", 
                         "Exclusive Relationship", "Open Relationship", "Unclear Relationship", 
                         "HIV Status", 
                         "Bathhouse", "Bathhouse Missing",
                         "Groupsex", "Groupsex Missing",
                         "Dating Apps", "Dating Apps Missing",
                         "Transactional Sex", "Transactional Sex Missing")
vars_model_fu_full <- vars_model_base_full[!vars_model_base_full %in% c("Dating Apps", 
                                                                        "Dating Apps Missing")]
tbl_coeff <- add_column(tbl_coeff,
                        name = rep(c("intercept", vars_model_base_full, "inverse_overdisp",
                                     rep(c("intercept", vars_model_fu_full, "inverse_overdisp"),2)), 
                                     times = 3),
                        .after = "term")


# table for posteriors (post-restriction only)
range(exp(tbl_coeff$`97.5%`))

tbl_coeff <- tbl_coeff %>% 
  mutate(name = factor(name, levels = unique(name))) %>%
  mutate(city = gsub("-Pre-Pandemic|-Pandemic|-Post-Restriction", "", data_pt),
         time_pt = gsub("Montreal-|Toronto-|Vancouver-", "", data_pt),
         .before = 1) %>%
  mutate(time_pt = factor(time_pt, 
                        levels = c("Pre-Pandemic", "Pandemic", "Post-Restriction")))

tbl_coeff_post <- tbl_coeff %>% subset(time_pt == "Post-Restriction")

coeff_post_tbl <- tbl_coeff_post %>% 
  dplyr::select(mean, city, time_pt, name, mean, se_mean, `2.5%`, `97.5%`) %>%
  mutate(SE = formatC(signif(se_mean,digits = 2), 
                      digits = 2,
                      format = "fg",
                      flag = "#"),
         Mean = paste0(round(exp(mean), 2), " (", round(exp(`2.5%`), 2), ", ", round(exp(`97.5%`), 2), ")")) %>%
  dplyr::select(city, name, Mean, SE)

coeff_post_tbl <- coeff_post_tbl[!coeff_post_tbl$name %in% c("Bathhouse Missing", "Groupsex Missing", "Dating Apps Missing", "Transactional Sex Missing"), ]
  

# write.csv(coeff_post_tbl, "./output-3cities/negbin-ipcw/table_coef_post.csv")
