# Libraries & data ----
library(tidyverse)

source("./src/seir_fns.R")

theme_set(theme_bw())

# partnership contact data
sex_  <- read.csv("./out/fitted-distributions-grouped/pmf_weighted_partn_p6m_20grps.csv")
data_pmf_prepand <- read_csv("./out/fitted-distributions-grouped/pmf_weighted_partn_p6m_20grps_prepandemic.csv")
data_case <- read.csv("./data-public/mpox_case_data_phac.csv")

PROV <- c("QC", "ON", "BC")
CITIES <- c("mtl", "tor", "van")

# Model running and calibration -----
# single province for troubleshooting
# province <- "ON"
# cty <- case_when(province == "QC" ~ "mtl",
#                  province == "ON" ~ "tor",
#                  province == "BC" ~ "van")

# format incidence data
data_case$date <- as.Date(data_case$date)

## Fit for each city ----
fit_ls <- list(mtl = NULL,
               tor = NULL,
               van = NULL)
theta_ls <- fit_ls
par_ls <- fit_ls
opt_ls <- fit_ls
incid_ls <- fit_ls
fixed_par_ls <- fit_ls
likdat_ls <- fit_ls

for(province in PROV){
  cty <- case_when(province == "QC" ~ "mtl",
                   province == "ON" ~ "tor",
                   province == "BC" ~ "van")
  
  print(paste("Running model for", cty, "=============="))
  
  # up to when to fit data
  date_end <- case_when(province == "QC" ~ as.Date("2022-06-15"),
                        province == "ON" ~ as.Date("2022-07-10"),
                        province == "BC" ~ as.Date("2022-07-10"))
  
  city_dat <- subset(data_case, prov == province & date < date_end); nrow(city_dat)
  print(mean(city_dat$incidence)); sd(city_dat$incidence)
  
  # population size and contact rates for model
  sex <- subset(sex_, city == cty)
  N <- ifelse(cty == "van", 26100, ifelse(cty == "tor", 78000, 54000))
  k <- nrow(sex)
  k_size <- sex$prop; if (sum(k_size) != 1) { print("error, pop size") }
  contact <- sex$mean_rate / 180
  avg_contact <- sum(contact * k_size)
  
  # initialize model population
  init_pop <- array(data = 0, dim = c(6, k))
  init_pop[1, ] <- N * k_size
  # imported_cases <- ifelse(cty == "van", 1, ifelse(cty == "mtl", 3, 4))
  imported_low <- ifelse(cty == "van", 1, ifelse(cty == "mtl", 2, 3))
  imported_upp <- ifelse(cty == "van", 3, ifelse(cty == "mtl", 6, 6))
  imported_start_val <- (imported_low + imported_upp) / 2
  
  # we import cases in the 5% highest activity group
  ind <- which(cumsum(k_size) >= 0.95)
  # # distributed equally within the E1, E2, I1, I2 compartments
  # init_pop[2, ind] <- imported_cases * k_size[ind] / sum(k_size[ind]) / 4
  # init_pop[3, ind] <- imported_cases * k_size[ind] / sum(k_size[ind]) / 4
  # init_pop[4, ind] <- imported_cases * k_size[ind] / sum(k_size[ind]) / 4
  # init_pop[5, ind] <- imported_cases * k_size[ind] / sum(k_size[ind]) / 4
  
  fixed_par <- list(init_pop = init_pop,
                    k_size = k_size,
                    seed_indx = ind,
                    seed_nb_init = imported_start_val,
                    seed_lbound = imported_low,
                    seed_ubound = imported_upp,
                    seed_diff = NA,
                    contact = contact,
                    alpha = 1 / 5.1,
                    # gamma = 1 / 4,
                    report_delay = 1 / 2,
                    start = 0, end = 150, dt = 0.25,
                    lag_introduction = ifelse(cty == "mtl", 21, 10))
  fixed_par$seed_diff <- (fixed_par$seed_ubound - fixed_par$seed_lbound)
  city_dat$time_intro <- city_dat$time_conti + fixed_par$lag_introduction # x days after introduction
  
  ### calibration
  theta0 <- c(qlogis(0.7), qlogis(0.5), qlogis(0.8), log(1/5),
              # to transform to logit scale from bounded scale do (val - lbound) / (ubound - lbound)
              qlogis( (fixed_par$seed_nb_init - fixed_par$seed_lbound) / fixed_par$seed_diff )
              )
  
  likdat <- city_dat
  llk(theta0, fixed_par, likdat)
  
  # Nelder-Mead better to get good starting values
  opt_nelder_mead <- optim(theta0, llk, fixed_par = fixed_par, likdat = likdat, method = "Nelder-Mead", 
                           control = list(fnscale = -1, trace = 0, maxit = 250), hessian = FALSE)
  # BFGS better to get the mode if supplied with good starting values
  opt <- optim(opt_nelder_mead$par, llk, fixed_par = fixed_par, likdat = likdat, method = "BFGS", 
               control = list(fnscale = -1, trace = 4, REPORT = 1, maxit = 250), hessian = TRUE)
  theta <- opt$par
  fit_par <- data.frame(
    parameter = c("transmission parameter (beta)", "assortativity (omega)", 
                  "reporting fraction", "duration infectiousness (1/gamma)",
                  "imported cases"),
    # to backtransform from logit scale do (lbound + val) / (ubound - lbound)
    value = round(
      c( plogis(opt$par[1:3]), 1 / exp(opt$par[4]), (fixed_par$seed_lbound + plogis(opt$par[5]) * fixed_par$seed_diff) ),
      2
    )
  ); print(fit_par)
  
  ## using the calibrated parameters get the point estimate for outputs
  fit <- mpox_mod(init_pop = fixed_par$init_pop,
                  contact = fixed_par$contact,
                  k_size = fixed_par$k_size,
                  seed_indx = fixed_par$seed_indx,
                  seed_nb = fixed_par$seed_lbound + plogis(opt$par[5]) * fixed_par$seed_diff,
                  beta = plogis(theta[1]), alpha = fixed_par$alpha, gamma = exp(theta[4]), omega = plogis(theta[2]),
                  report_frac = plogis(theta[3]), report_delay = fixed_par$report_delay,
                  start = fixed_par$start, end = fixed_par$end, dt = fixed_par$dt)
  
  fit$date <- fit$time - fixed_par$lag_introduction + min(city_dat$date)
  
  ## store fitted results
  # model parameters
  theta_ls[[cty]] <- theta
  par_ls[[cty]] <- fit_par
  opt_ls[[cty]] <- opt
  
  # model fit
  fit_ls[[cty]] <- fit
  incid_ls[[cty]] <- city_dat
  
  # model inputs
  fixed_par_ls[[cty]] <- fixed_par
  likdat_ls[[cty]] <- likdat
}

## check fit results
province <- "ON"
cty <- CITIES[province == PROV]
plot(fit_ls[[cty]]$cases ~ fit_ls[[cty]]$date, type = 'l', lwd = 3, col = "firebrick4",
     main = paste0("mpox in ", province), ylim = c(0, max(incid_ls[[cty]]$incidence * 1.2)), 
     xlab = "time", ylab = "reported mpox cases (day)")
  lines(incid_ls[[cty]]$incidence ~ incid_ls[[cty]]$date, type = "s", col = "steelblue4", lwd = 2)

# saveRDS(fit_ls, sprintf("./out-seir/%s_model_fit_out.rds", Sys.Date()))
# saveRDS(par_ls, sprintf("./out-seir/%s_model_fit_pars.rds", Sys.Date()))

## Credible interval simulations ----
# simulate CrI
# suggest nsir = 5000 and sim = 1000 when SIR is TRUE. Can reduce when exploring results
ci_ls <- list(mtl = NULL,
              tor = NULL,
              van = NULL)

t0 <- Sys.time()
for(cty in CITIES){
  print(cty)
  ci_ls[[cty]] <- simul_fun(opt_ls[[cty]]$hessian, theta_ls[[cty]], fixed_par_ls[[cty]], sim = 1000,
                            SIR = TRUE, likdat = likdat_ls[[cty]], nsir = 15000, with_replacement = TRUE,
                            parallel = TRUE)
  print(ci_ls[[cty]]$posterior_ci)
  
  ci_ls[[cty]]$result$date <- ci_ls[[cty]]$result$time - fixed_par_ls[[cty]]$lag_introduction + min(incid_ls[[cty]]$date)
}
t1 <- Sys.time()
t1 - t0 # ~1hr with nsir=15000

# saveRDS(ci_ls, sprintf("./out-seir/%s_CIs.rds", Sys.Date()))
# ci_ls <- readRDS(sprintf("./out-seir/%s_CIs.rds", Sys.Date()))

### Verification of results ----
## plotting results with CIs
province <- "BC"
cty <- CITIES[province == PROV]
plot(fit_ls[[cty]]$cases ~ fit_ls[[cty]]$date , type = 'l', lwd = 3, col = "firebrick4",
     main = paste0("mpox in ", province), ylim = c(0, max(incid_ls[[cty]]$incidence * 1.2)),
     xlab = "time", ylab = "reported mpox cases (day)")
  lines(incid_ls[[cty]]$incidence ~ incid_ls[[cty]]$date, type = "s", col = "steelblue4", lwd = 2)
  polygon(x = c(ci_ls[[cty]]$result$date, rev(ci_ls[[cty]]$result$date)),
          y = c(ci_ls[[cty]]$result$lci,  rev(ci_ls[[cty]]$result$uci)),
          border = NA, col = rgb(245, 160, 142, 125, max = 255))
  lines(fit_ls[[cty]]$cases ~ fit_ls[[cty]]$date, lwd = 3, col = "firebrick4")
  legend("topleft", col = c("steelblue4", "firebrick4", rgb(245, 160, 142, 125, max = 255)),
         legend = c("Reported mpox cases", "Modeled cases", "95% credible intervals"),
         lwd = c(3, 3, NA), pch = c(NA, NA, 15), pt.cex = c(NA, NA, 3), bty = "n")

## plot proportion recovered by risk groups
province <- "ON"
cty <- CITIES[province == PROV]

prp_rec <- t(t(fit_ls[[cty]]$X[, 6, ]) / fit_ls[[cty]]$X[1, 1, ]) * 100

plot(prp_rec[ , ncol(prp_rec)] ~ fit$date, type = 'l', lwd = 0, xlab = "time (days)", ylab = "Proportion immune (%)", ylim = c(0, 100),
     main = paste0("Proportion with natural immunity in the ", ncol(prp_rec), " risk groups (", province, ")"))
  pal <- wesanderson::wes_palette("Zissou1", ncol(prp_rec), type = "continuous")
  for (i in 1:ncol(prp_rec)) {
    lines(prp_rec[ , i] ~ fit_ls[[cty]]$date, lwd = 2, col = pal[i])
  }

# total cumulative risk
(sum(fit_ls[[cty]]$X[1, , ]) - sum(fit_ls[[cty]]$X[length(fit_ls[[cty]]$time), 1, ])) / sum(fit_ls[[cty]]$X[1, , ])
# total cumulative case
(sum(fit_ls[[cty]]$X[1, , ]) - sum(fit_ls[[cty]]$X[length(fit_ls[[cty]]$time), 1, ]))

## Compute R0 ----
### Estimated R0 (post-restrictions) ----
df_r0_ls <- list(mtl = NULL,
                 tor = NULL,
                 van = NULL)

prp_imm <- cumsum(rev(k_size))

for(cty in CITIES){
  df_r0 <- data.frame()
  
  for(grp_rm in 0:which(prp_imm == 0.05)){
    if(grp_rm == 0){
      # cur_r0 <- r0_fun(fit_ls[[cty]], theta = theta_ls[[cty]], fixed_par = fixed_par_ls[[cty]], compute_ci = FALSE) # pt estimate
      cur_r0 <- r0_fun(fit_ls[[cty]], theta = theta_ls[[cty]], fixed_par = fixed_par_ls[[cty]],
                       post_samples = ci_ls[[cty]]$posterior_samples)
      cur_r0 <- mutate(cur_r0, prop_imm = 0, .before = 1)
    } else {
      # cur_r0 <- r0_fun(fit_ls[[cty]], theta = theta_ls[[cty]], fixed_par = fixed_par_ls[[cty]], grp_immunity = grp_rm, compute_ci = FALSE) # pt estimate
      cur_r0 <- r0_fun(fit_ls[[cty]], theta = theta_ls[[cty]], fixed_par = fixed_par_ls[[cty]],
                       post_samples = ci_ls[[cty]]$posterior_samples, grp_immunity = grp_rm)
      cur_r0 <- mutate(cur_r0, prop_imm = prp_imm[grp_rm], .before = 1)
    }
    
    # add to df
    df_r0 <- bind_rows(df_r0, cur_r0)
    df_r0 <- mutate(df_r0, city = cty, .before = 1)
  }
  
  df_r0_ls[[cty]] <- df_r0
}
df_r0 <- bind_rows(df_r0_ls)

head(df_r0)
rm(df_r0_ls)

### Projected R0 (pre-pandemic) ----
## need to change contact rates in fixed_par
# compute and save R0
df_r0_pre_ls <- list(mtl = NULL,
                     tor = NULL,
                     van = NULL)

for(cty in CITIES){
  # set contact rates to pre-pandemic levels
  fixed_par_tmp <- list()
  fixed_par_tmp$contact <- subset(data_pmf_prepand, city == cty)$mean_rate / 180
  
  df_r0_pre <- data.frame()
  
  for(grp_rm in 0:which(prp_imm == 0.06)){
    if(grp_rm == 0){
      # cur_r0_pre <- r0_fun(fit_ls[[cty]], theta_ls[[cty]], fixed_par_tmp, compute_ci = FALSE) # pt estimate
      cur_r0_pre <- r0_fun(fit_ls[[cty]], theta_ls[[cty]], fixed_par_tmp,
                           post_samples = ci_ls[[cty]]$posterior_samples)
      cur_r0_pre <- mutate(cur_r0_pre, prop_imm = 0, .before = 1)
    } else {
      # cur_r0_pre <- r0_fun(fit_ls[[cty]], theta_ls[[cty]], fixed_par_tmp, grp_immunity = grp_rm, compute_ci = FALSE) # pt estimate
      cur_r0_pre <- r0_fun(fit_ls[[cty]], theta_ls[[cty]], fixed_par_tmp,
                           post_samples = ci_ls[[cty]]$posterior_samples, grp_immunity = grp_rm)
      cur_r0_pre <- mutate(cur_r0_pre, prop_imm = prp_imm[grp_rm], .before = 1)
    }
    
    # add to df
    df_r0_pre <- bind_rows(df_r0_pre, cur_r0_pre)
    df_r0_pre <- mutate(df_r0_pre, city = cty, .before = 1)
  }
  
  df_r0_pre_ls[[cty]] <- df_r0_pre
}
rm(fixed_par_tmp)

df_r0_pre <- bind_rows(df_r0_pre_ls)
head(df_r0_pre)

# Format model outputs for paper ----
## Incidence ----
df_incid <- data.frame()

for(cty in CITIES){
  df_incid <- bind_rows(
    df_incid,
    data.frame(city = cty,
               time = fit_ls[[cty]]$time,
               date = fit_ls[[cty]]$date,
               cases = fit_ls[[cty]]$cases,
               cases_lci = ci_ls[[cty]]$result$lci,
               cases_uci = ci_ls[[cty]]$result$uci,
               incidence = fit_ls[[cty]]$inc)
  )
}

head(df_incid)
write.csv(df_incid, "./out-seir/fit_incidence.csv", row.names = FALSE)

## Parameters ----
df_pars <- data.frame()

for(cty in CITIES){
  # point estimate and CrIs
  df_par_tmp <- full_join(
    par_ls[[cty]],
    ci_ls[[cty]]$posterior_ci,
    by = c("parameter" = "names")
  )
  
  df_par_tmp <- mutate(df_par_tmp, city = cty, .before = 1)
  
  # save in final df
  df_pars <- bind_rows(df_pars, df_par_tmp)
}

# get average
df_pars_avg <- df_pars %>% 
  mutate(parameter = factor(parameter, levels = unique(parameter))) %>% 
  group_by(parameter) %>% 
  summarize(value = mean(value)) %>% 
  mutate(city = "all", .before = 1)

df_pars <- bind_rows(df_pars_avg, df_pars); rm(df_pars_avg)

df_pars
write.csv(df_pars, "./out-seir/fit_pars.csv", row.names = FALSE)

### Format parameters ----
# parameter order for final output
df_pars$parameter <- factor(df_pars$parameter,
                            levels = c("duration infectiousness (1/gamma)",
                                       "transmission parameter (beta)",
                                       "assortativity (omega)",
                                       "reporting fraction",
                                        "imported cases"))
df_pars <- arrange(df_pars, parameter, city)

# put 'all' mean at the start
df_pars <- bind_rows(filter(df_pars, city == "all"), filter(df_pars, city != "all"))

# create CrI column
df_pars[, 3:5] <- round(df_pars[, 3:5], 2)
df_pars$cri <- sprintf("(%s\u2013%s)", df_pars$lci, df_pars$uci)
df_pars$cri[is.na(df_pars$lci)] <- NA_character_

# save table
write.csv(df_pars, "./out/manuscript-tables/table_S2_seir_params.csv", row.names = FALSE)

## R0 and R_eff ----
### post-restrictions R0
head(df_r0)
write.csv(df_r0, "./out-seir/fit_r0.csv", row.names = FALSE)

### pre-pandemic projected R0
head(df_r0_pre)
write.csv(df_r0_pre, "./out-seir/fit_r0_prepand.csv", row.names = FALSE)
