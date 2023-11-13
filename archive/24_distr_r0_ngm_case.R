# Libraries and data ----
library(tidyverse)
library(ggtext)
library(parallel)
library(foreach)
library(doParallel)

theme_set(theme_bw())

source("./src/utils_r0_estim.R")
source("./src/utils_helper.R")

## degree distribution for NGM
data_pmf_iter <- readRDS("./out/pmf_stan_iterations.rds")

# for troubleshooting with only point estimate
data_full <- read_csv("./out/fitted-distributions/pmf_weighted_all_partn.csv")

## case data for R0 based on growth rate
case_data <- read.csv("../mpx-engage-params/misc-grant-app/data-case/monkeypox-detailed-2023june13.csv")
case_data$reporting_pt_en[case_data$reporting_pt_en == "Quebec"] <- "Québec"

# only need the 3 concerned provinces
PROVS <- c("Québec", "Ontario", "British Columbia")
case_data <- case_data %>%
  subset(reporting_pt_en %in% PROVS)

CITIES <- c("Montreal", "Toronto", "Vancouver")

# R0 from NGM -----
# range of beta/SAR
#' we use a fine-scale (0.01 around the area where the NGM and case-based methods intersect)
#' to find the equivalence between the two
#' and a coarser-scale (0.05-0.10) for the rest, to complete the plot
beta_range <- unique( c( seq(0, 0.3, 0.1), seq(0.3, 0.5, 0.01), seq(0.5, 1, 0.1) ) )

# list to store R0 results (post-restrictions)
R0_ls_city <- list(grp100 = create_city_list(CITIES), 
                grp248 = create_city_list(CITIES))

# start a cluster
numCores <- parallel::detectCores()
numCores
numCores <- numCores - 1
cl <- parallel::makeCluster(numCores, type = "PSOCK")

# register
doParallel::registerDoParallel(cl = cl)

t0 <- Sys.time()
for(categorization in c("grp100", "grp248")){
  
  for(cur_city in CITIES){
  print(cur_city)
  # create data-pt label and extract city PMF
  # data_pt <- paste(cur_city, "-Post-Restrictions", sep = "")
  # pmf_tmp <- data_pmf_iter[[data_pt]]
  # 
  # R0_ls_city[[categorization]][[cur_city]] <-  foreach( cur_iter = 1:nrow(pmf_tmp) ) %dopar% {
  #   compute_r0_ngm(
  #     pmf_estim_m = pmf_tmp[cur_iter, ],
  #     beta_range = beta_range,
  #     run_par = TRUE,
  #     finer = (categorization == "248"))
  #   )
  # }
  
  # troubleshooting code, only uses the point estimates
  data_pt <- paste(cur_city, "-Post-Restrictions", sep = "")
  R0_ls_city[[categorization]][[cur_city]] <-  compute_r0_ngm(
    pmf_estim_m = data_full[data_full$data_pt == data_pt, ]$mean,
    beta_range = beta_range,
    run_par = FALSE,
    finer = (categorization == "grp248"))
}}
t1 <- Sys.time()
t1 - t0 # 1.03 hours

R0_ls_city

rm(data_pt, pmf_tmp)

stopCluster(cl)

# save iterations
# saveRDS(R0_ls_city, "./out/r0-estim-iterations.rds")

# R0 from cases ----
## format data, to get time from first case in each province and ln(cumulative cases)
case_data <- case_data %>% 
  group_by(reporting_pt_en) %>%
  mutate(
    earliest_date = as.Date(min(date)),
    time_conti = as.integer(as.Date(date) - earliest_date),
    ln_cumu_cases = log(num_confirmedcases_cumulative),
    ln_incid = log(num_confirmedcases_delta),
    reporting_pt_en = factor(reporting_pt_en, level = PROVS)
  )

# cutoff data to use in growth rate estimation for R0
cutoff <- 50
case_data_cutoff <- case_data %>% 
  filter(time_conti <= cutoff & reporting_pt_en == "British Columbia" |
           time_conti <= cutoff & reporting_pt_en == "Ontario" |
           time_conti <= cutoff & reporting_pt_en == "Québec")

## plot cumulative incidence and fitted log-growth rate
png("./fig/fig_S2_cumul_incidence.png",
    width = 15, height = 7, units = "cm", res = 600)
ggplot(case_data, aes(x = time_conti,
                      y = ln_cumu_cases,
                      group = reporting_pt_en,
                      col = reporting_pt_en)) +
  geom_point(size = 0.6, alpha = 0.8) + 
  # geom_line(aes(y = ln_incid), linetype = "dashed") + # can use to compare growth in cumulative cases vs incidence
  # plot a straight line for the data used in the regression
  # this is equivalent to doing lm(ln_cumu_cases ~ time_conti) for each province
  geom_smooth(data = case_data_cutoff, method = "lm", se = FALSE, linewidth = 0.5) +
  # change axis titles
  
  coord_cartesian(xlim = c(0, 200)) +
  scale_colour_viridis_d(option = "C", end = 0.8)  +
  labs(x = "Days since first case",
       y = "ln(Cumulative confirmed mpox cases)",
       col = "Province") +
  
  theme(axis.title = element_text(size = 9),
        axis.text = element_text(size = 7),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 7))
dev.off()

lambda <- vector(mode = "numeric", 
                 length = 3)
names(lambda) <- PROVS
for (province in PROVS){
  case_data_province <- case_data_cutoff %>% 
    subset(reporting_pt_en == province)
  linear_model <- lm(formula = ln_cumu_cases ~ time_conti, case_data_province)
  # linear_model <- lm(formula = incid ~ time_conti, case_data_province)
  lambda[[province]] <- linear_model$coefficients[2] %>%
    unname()
}

# formula is (1 + lambda*D)(1 + lambda*D') where D: infectious period & D': 
SI = 8.99
D_apostr = (7.12 - 2.5)
D = SI - D_apostr

r0_case <- (1 + lambda * D) * (1 + lambda * D_apostr)
round(r0_case, 2)

# simplified formula, should roughly correspond but use above in main paper
# round( (1 + lambda * SI) , 2)

write.csv(r0_case, "./out/r0-estim-cases.csv", row.names = FALSE)

# Estimate SAR and R0 given pre-pandemic sexual activity level -------
SAR_by_city <- list(grp100 = create_city_list(CITIES),
                    grp248 = create_city_list(CITIES))
R0_by_city <- list(grp100 = create_city_list(CITIES),
                   grp248 = create_city_list(CITIES))

for(categorization in c("grp100", "grp248")){
  for(cur_city in CITIES){
  # SAR estimation and R0 CrI's ----
  
  ## Match the case R0 to the NGM R0 ----
  # arrange R0 results as [iterations, SAR] for the specified city
  # r0_mtx <- do.call(rbind, R0_ls_city_248[[categorization]][[cur_city]])
  r0_mtx <- R0_ls_city[[categorization]][[cur_city]] # for troubleshooting with only point esimate

  ### find closest point to match to case-based R0
  # get R0 from cases, keep only relevant city
  prov_keep <- case_when(grepl("Montreal", cur_city) ~ "Québec",
                         grepl("Toronto", cur_city) ~ "Ontario",
                         grepl("Vancouver", cur_city) ~ "British Columbia")
  r0_case_city <- r0_case[prov_keep]
  
  ## R0 from NGM
  # first get difference between NGM and case to find closest value
  r0_mtx_diff <- r0_mtx - r0_case_city
  
  # find minimum in each iteration (row)
  # r0_mtx_mins <- apply(abs(r0_mtx_diff), 1, which.min) # comment out for troubleshooting w/ only pt estimate
  
  # for troubleshooting
  r0_mtx_mins <- which.min(abs(r0_mtx_diff))
  
  # get the corresponding SARs
  SAR_by_city[[categorization]][[cur_city]] <- tibble( 
                                    ngrp = categorization,
                                    city = cur_city,
                                    # iter = 1:nrow(r0_mtx), # comment out for troubleshooting w/ only pt estimate
                                    SAR = beta_range[r0_mtx_mins]) 
  
  ## Get R0 values for plotting ----
  # get credible intervals
  # cred_l <- vector("double", ncol(r0_mtx))
  # cred_u <- vector("double", ncol(r0_mtx))
  
  # troubleshooting
  cred_l <- r0_mtx 
  cred_u <- r0_mtx
  
  # for(i in 1:ncol(r0_mtx)){  # comment out for troubleshooting w/ only pt estimate
  #   cred_l[i] <- quantile(r0_mtx[, i], .025)
  #   cred_u[i] <- quantile(r0_mtx[, i], .975)
  # }
  
  # create tibble with each city and time period
  R0_by_city[[categorization]][[cur_city]] <- tibble(
                                   ngrp = categorization,
                                   city = cur_city,
                                   SAR = beta_range,
                                   # mean_r0 = colSums(r0_mtx) / nrow(r0_mtx), 
                                   mean_r0 = r0_mtx, # troubleshooting
                                   r0.l = cred_l,
                                   r0.u = cred_u)
    
}}
rm(r0_mtx, r0_mtx_diff, r0_mtx_mins)

## Get CrI's for SAR ----
# df_sar <- bind_rows(SAR_by_city)
df_sar <- map_df(SAR_by_city, bind_rows) # troubleshoot

# compute city-specific SAR
df_sar_city <- df_sar %>% 
  group_by(ngrp, city) %>% 
  summarize(
    mean = mean(SAR),
    cr.i_low = quantile(SAR, .025),
    cr.i_upp = quantile(SAR, .975)
  )

# compute average SAR
df_sar_avg <- df_sar %>% 
  group_by(ngrp) %>%
  summarize(
    city = "avg",
    mean = mean(SAR),
    cr.i_low = quantile(SAR, .025),
    cr.i_upp = quantile(SAR, .975)
  )

df_sar_all <- bind_rows(df_sar_avg, df_sar_city)
df_sar_all

write.csv(df_sar_all, "./out/r0-estim-ngm/SAR_estimates.csv", row.names = FALSE)

# Project R0 using estimated SAR ----
# use point estimate
SAR_avg <- df_sar_avg$mean
names(SAR_avg) <- c("grp100", "grp248")
SAR_avg

# list to store R0 results (pre-pandemic)
R0_projected_ls <-  list(grp100 = create_city_list(CITIES), 
                         grp248 = create_city_list(CITIES))

# start a cluster
numCores <- parallel::detectCores()
numCores
numCores <- numCores - 1
cl <- parallel::makeCluster(numCores, type = "PSOCK")

# register
doParallel::registerDoParallel(cl = cl)

t0 <- Sys.time()
for(categorization in c("grp100", "grp248")){
for(cur_city in CITIES){
  print(cur_city)
  # create data-pt label and extract city PMF
  data_pt <- paste(cur_city, "-Pre-Pandemic", sep = "")
  # pmf_tmp <- data_pmf_iter[[data_pt]]
  # 
  # R0_projected_ls[[categorization]][[cur_city]] <-  foreach( cur_iter = 1:nrow(pmf_tmp) ) %dopar% {
  #   compute_r0_ngm(
  #     pmf_estim_m = pmf_tmp[cur_iter, ],
  #     beta_range = SAR_avg,
  #     run_par = TRUE)
      
  # troubleshooting code, only uses the point estimates
  R0_projected_ls[[categorization]][[cur_city]] <-  compute_r0_ngm(
        pmf_estim_m = data_full[data_full$data_pt == data_pt, ]$mean,
        beta_range = SAR_avg[categorization],
        run_par = TRUE,
        finer = (categorization == "grp248"))
    
  }
}
t1 <- Sys.time()
t1 - t0 # 4.8 minutes
rm(data_pt, pmf_tmp)

stopCluster(cl)

## Get point estim. and CrI ----
df_r0_proj <- list(grp100 = create_city_list(CITIES), 
                   grp248 = create_city_list(CITIES))
for(categorization in c("grp100", "grp248")){
for(cur_city in CITIES){
  # r0_mtx <- do.call(rbind, R0_projected_ls[[categorization]][[cur_city]])
  r0_mtx <- df_r0_proj[[categorization]][[cur_city]] # for troubleshooting with only point esimate
  
  ## Get R0 values for plotting ----
  # # get credible intervals 
  # cred_l <- vector("double", ncol(r0_mtx))
  # cred_u <- vector("double", ncol(r0_mtx))
  
  # troubleshooting
  cred_l <- r0_mtx 
  cred_u <- r0_mtx
  
  # for(i in 1:ncol(r0_mtx)){ # comment out for troubleshooting w/ only pt estimate
  #   cred_l[i] <- quantile(r0_mtx[, i], .025)
  #   cred_u[i] <- quantile(r0_mtx[, i], .975)
  # }
  
  # create tibble with each city and time period
  df_r0_proj[[categorization]][[cur_city]] <- tibble(ngrp = categorization,
                                                     city = cur_city,
                                   SAR = SAR_avg,
                                   # mean_r0 = colSums(r0_mtx) / nrow(r0_mtx), 
                                   mean_r0 = r0_mtx, # troubleshooting
                                   r0.l = cred_l,
                                   r0.u = cred_u)
  
}}

## R0 dataframes for plotting ----
# df_r0_ngm <- bind_rows(R0_by_city)
# df_r0_proj <- bind_rows(df_r0_proj)
df_r0_ngm <- map_df(R0_by_city, bind_rows)
df_r0_proj <- map_df(df_r0_proj, bind_rows)

# save estimated NGM R0 and projected R0 using average SAR
write.csv(df_r0_ngm, "./out/r0-estim-ngm/R0_estim_NGM.csv",
          row.names = FALSE)

write.csv(df_r0_proj, "./out/r0-estim-ngm/R0_projected.csv",
          row.names = FALSE)
