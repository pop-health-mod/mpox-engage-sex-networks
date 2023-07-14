# library and data--------------------------------------------------------------
library(tidyverse)
source("./src/utils_helper.R")
library(dplyr)
library(ggtext)

data_full <- read_csv("./output-3cities/negbin-ipcw/full_fitted_pmf_wt_p6m_all.csv")
# TODO: MOVE SUPPLEMENTARY TO DIFFERENT SCRIPT
data_restr <- read_csv("./output-3cities/negbin-res/full_fitted_cdf_wt_p6m_all.csv")
data_age <- read_csv("./output-3cities/negbin-age/full_fitted_cdf_wt_p6m_all.csv")
data_anal <- read_csv("./output-3cities/negbin-ipcw/full_fitted_cdf_wt_p6m_anal.csv")
case_data <- read.csv("./misc-grant-app/data-case/monkeypox-detailed-2023june13.csv")

CITIES <- c("Montreal", "Toronto", "Vancouver")
# DATASET <- c("data_full", "data_restr", "data_age", "data_anal")
DATASET <- c("data_full")
df_draw <- create_city_list(DATASET)
df_point = df_draw
SAR_est = df_draw
SAR_est.lb = df_draw
SAR_est.ub = df_draw

# R0 cases ----------------------------------------------------------------

PROVS <- c( "Quebec","Ontario","British Columbia")

# plot cumulative cases vs. t, by province
case_data_pro <- case_data %>%
  subset(reporting_pt_en %in% PROVS)

## x-axis == date
ggplot(case_data_pro,aes(x = as.Date(date),
                         num_confirmedcases_cumulative, 
                         group = reporting_pt_en,
                         col = reporting_pt_en)) +
  geom_line() + 
  xlab("Date") + 
  ylab("Number of Cumulative Confirmed Cases") + 
  guides(col = guide_legend(title = "Province")) +
  scale_x_date(date_labels = "%b %y", date_breaks  = "1 month")

## x-axis == time from (first case in each province), y-axis == ln(cumulative cases)
case_data_pro <- case_data_pro %>% 
  group_by(reporting_pt_en) %>%
  mutate(earliest_date = as.Date(min(date)),
         time_conti = as.integer(as.Date(date) - earliest_date),
         ln_cumu_cases = log(num_confirmedcases_cumulative),
         reporting_pt_en = factor(reporting_pt_en, level = c("Quebec", "Ontario", "British Columbia")))

# plot log cumulative cases with a fitted line for the data used 
case_data_cutoff <- case_data_pro %>% 
  filter(time_conti <= 100 & reporting_pt_en == "British Columbia" |
           time_conti <= 100 & reporting_pt_en == "Ontario" |
           time_conti <= 100 & reporting_pt_en == "Quebec") 


png("./figures-3cities/negbin-ipcw/figS_r0_case.png",
    width = 12, height = 12, units = "cm", res = 600)

ggplot(case_data_pro, aes(x = time_conti,
                         y = ln_cumu_cases,
                         group = reporting_pt_en,
                         col = reporting_pt_en)) +
  geom_point() + 
  # plot a straight line for the data used in the regression
  # this is equivalent to doing lm(ln_cumu_cases ~ time_conti) for each province
  # geom_smooth(data = case_data_cutoff, method = "lm", se = FALSE) +
  # change axis titles
  xlab("Days since first case") + 
  ylab("ln(Number of cumulative confirmed mpox cases)") + 
  guides(col = guide_legend(title = "Province")) +
  scale_colour_viridis_d(option = "C", end = 0.8)  +
  theme(plot.title.position = "plot",
        plot.caption.position = "plot",
        plot.caption = element_markdown(hjust = 0, size = 5.5))
  # labs(caption = c("**Figure S1. Natural logarithm of cumulative comfirmed mpox cases in the province of Quebec, Ontario, and British Columbia<br>against time based on reported mpox cases updated on June 13th, 2023.** R0 was estimated as the gradient of the simple linear<br>regression fitted line during the first 100 days after the first mpox case was reported in each province."))

dev.off()

lambda <- vector(mode = "numeric", 
                 length = 3)
names(lambda) <- PROVS
for (province in PROVS){
  case_data_province <- case_data_cutoff %>% 
    subset(reporting_pt_en == province)
  linear_model <- lm(formula = ln_cumu_cases ~ time_conti, case_data_province)
  lambda[[province]] <- linear_model$coefficients[2] %>%
    unname()
}

r0_case <- 1 + lambda * 17.3

# R0 NGM ------------------------------------------------------------------

data_pt_name = unique(data_full$data_pt)

# DATASET <- "data_full" # paper result

for(dataset_name in DATASET){
  # contact rates and group density (dij = Nij/Nj where i is a sexual activity group, j is a city-timept)
  data_pmf_summ <- get(dataset_name) %>% 
    dplyr::select(data_pt, y_pred, mean, cr.i_low, cr.i_upp)
  print(dataset_name)
  # 100 sexual activity groups
  for (j in 1:99) {
    i = 99 - j + 2
    data_pmf_summ <- data_pmf_summ %>%
    mutate(quant_grp = ifelse(y_pred <= 3 * i, paste0("grp_", i), quant_grp))
  }
  
  data_pmf_summ <- data_pmf_summ %>%
  mutate(quant_grp = case_when(
    y_pred == 0         ~ "grp_1",
    TRUE ~ quant_grp))
  
  # add group density
  data_pmf_summ <- data_pmf_summ %>% 
    group_by(data_pt, quant_grp) %>% 
    summarize(y_mean = mean(y_pred), 
              qt = sum(mean), 
              qt.l = sum(cr.i_low),
              qt.u = sum(cr.i_upp),
              .groups = "drop")
  
  n = length(unique(data_pmf_summ$quant_grp)) #number of sexual activity groups
  beta_range = seq(0, 1, 0.01)
  
  # compute the contact rates
  # c: number of sexual contacts per day (estimate from P6M data)
  # d: proportion of the population in that sexual activity group
  contact <- data_pmf_summ %>% 
    transmute(data_pt, 
              c = y_mean / (6*30), 
              d = qt, 
              d.l = qt.l, 
              d.u = qt.u)
  
  # Delta matrix (a diagonal matirx with d's on diagonal)
  Delta <- create_city_list(data_pt_name)
  Delta.l <- create_city_list(data_pt_name)
  Delta.u <- create_city_list(data_pt_name)
  for(pt in data_pt_name){
    Delta[[pt]] <- matrix(0, n, n)
    Delta.l[[pt]] <- matrix(0, n, n)
    Delta.u[[pt]] <- matrix(0, n, n)
    contact_pt <- contact %>% subset(data_pt == pt)
    for(i in 1:n){
      Delta[[pt]][i, i] <- contact_pt$d[i]
      Delta.l[[pt]][i, i] <- contact_pt$d.l[i]
      Delta.u[[pt]][i, i] <- contact_pt$d.u[i]
    }
  }
  
  # TODO: this should be a function like r0_estim_ngm
  # or multiple functions, first one being a file like r0_var_setup to
  # instantiate all the variables, then 
  # matrices to store the R0
  R0 <- matrix(0, length(beta_range), length(data_pt_name))
  R0.l <- matrix(0, length(beta_range), length(data_pt_name))
  R0.u <- matrix(0, length(beta_range), length(data_pt_name))
  
  colnames(R0) <- data_pt_name
  colnames(R0.l) <- data_pt_name
  colnames(R0.u) <- data_pt_name
  
  # incubation period
  v <- 1/(7.9)
  
  # disease duration
  gamma <- 1/(17.3)
  
  # matrix 
  
  m = 2 * n # need E and I compartments for each group
  Sigma = matrix(0, m, m)
  E = matrix(0, m, n)
  
  
  # matrice A to calculate in final size estimation
  # TODO: delete?
  A <- create_city_list(data_pt_name)
  A.l <- create_city_list(data_pt_name)
  A.u <- create_city_list(data_pt_name)
  
  # Transition matrix sigma
  for(i in 1:n){
    Sigma[2*i - 1, 2*i - 1] <- -v
    Sigma[2*i, 2*i - 1] <- v
    Sigma[2*i, 2*i] <- -gamma
  }
  
  # Auxiliary matrix E
  for(i in 1:n){
    E[2*i - 1, i] <- 1
  }
  
  for(pt in data_pt_name){
    print(pt)
    
    # Transmission matrix T
    T_pt = matrix(0, m, m)
    T_pt.l = matrix(0, m, m)
    T_pt.u = matrix(0, m, m)
    
    # c: city and time-specific sexual activity rate
    # d: proportion in sexual activity group (instead of using population size N)
    c_pt = contact[contact$data_pt == pt, ]$c
    d_pt = contact[contact$data_pt == pt, ]$d
    d_pt.l = contact[contact$data_pt == pt, ]$d.l
    d_pt.u = contact[contact$data_pt == pt, ]$d.u
    
    # c_ttl:    total sexual activity rate
    #           c_1 + c_2 + c_3 + ... + c_n
    # denom_pt: total number of partnerships generated in the population
    #           c_1*d_1 + c_2*d_2 + ... c_n*d_n
    c_ttl = sum(c_pt)
    
    denom_pt = sum(c_pt * d_pt)
    denom_pt.l = sum(c_pt * d_pt.l)
    denom_pt.u = sum(c_pt * d_pt.u)
    
    #' the columns are in the order of E_i -> I_i -> E_i+1 ....
    #' where i is the sexual activity group
    for(beta in beta_range){
      for(i in 1:n){ 
        for(j in 1:n){
          T_pt[2 * i - 1, 2 * j] <- c_pt[i] * d_pt[i] * c_pt[j] * beta / denom_pt
          T_pt.l[2 * i - 1, 2 * j] <- c_pt[i] * d_pt.l[i] * c_pt[j] * beta / denom_pt.l
          T_pt.u[2 * i - 1, 2 * j] <- c_pt[i] * d_pt.u[i] * c_pt[j] * beta / denom_pt.u
        }
      }
      K_pt = -t(E) %*% T_pt %*% solve(Sigma) %*% E
      K_pt.l = -t(E) %*% T_pt.l %*% solve(Sigma) %*% E
      K_pt.u = -t(E) %*% T_pt.u %*% solve(Sigma) %*% E
      
      # verified that using NGM with large domain K_l = -T *sigma^-1
      # gives same result as below
      # K_alt <- -T_pt %*% solve(Sigma)
      # Re(eigen(K_alt)$values[1])
      
      R0_pt <- Re(eigen(K_pt)$values[1])
      R0_pt.l <- Re(eigen(K_pt.l)$values[1])
      R0_pt.u <- Re(eigen(K_pt.u)$values[1])
      
      R0[match(beta, beta_range), pt] <- R0_pt
      R0.l[match(beta, beta_range), pt] <- R0_pt.l
      R0.u[match(beta, beta_range), pt] <- R0_pt.u
      
      
      A[[pt]] <- solve(Delta[[pt]]) %*% K_pt %*% Delta[[pt]]
      A.l[[pt]] <- solve(Delta.l[[pt]]) %*% K_pt.l %*% Delta.l[[pt]]
      A.u[[pt]] <- solve(Delta.u[[pt]]) %*% K_pt.u %*% Delta.u[[pt]]
      
    }
  }
  
  R0 <- data.frame(SAR = beta_range, R0)
  R0.l <- data.frame(SAR = beta_range, R0.l)
  R0.u <- data.frame(SAR = beta_range, R0.u)
  
  # long format for easier plotting in ggplot
  # one row per R0 estimate, columns for city and time point
  dat_r0 <- R0
  dat_r0.l <- R0.l
  dat_r0.u <- R0.u
  
  names(dat_r0) <- gsub("Pre\\.Pandemic", "Pre-Pandemic", names(dat_r0))
  names(dat_r0) <- gsub("Post\\.Restriction", "Post-Restriction", names(dat_r0))
  
  names(dat_r0.l) <- gsub("Pre\\.Pandemic", "Pre-Pandemic", names(dat_r0.l))
  names(dat_r0.l) <- gsub("Post\\.Restriction", "Post-Restriction", names(dat_r0.l))
  
  names(dat_r0.u) <- gsub("Pre\\.Pandemic", "Pre-Pandemic", names(dat_r0.u))
  names(dat_r0.u) <- gsub("Post\\.Restriction", "Post-Restriction", names(dat_r0.u))
  
  dat_r0 <- dat_r0 %>% 
    pivot_longer(cols = `Montreal.Pre-Pandemic`:`Vancouver.Post-Restriction`,
                 names_to = c("city", "time_pt"), names_sep = "\\.", values_to = "r0") 
  dat_r0.l <- dat_r0.l %>% 
    pivot_longer(cols = `Montreal.Pre-Pandemic`:`Vancouver.Post-Restriction`,
                 names_to = c("city", "time_pt"), names_sep = "\\.", values_to = "r0") 
  dat_r0.u <- dat_r0.u %>% 
    pivot_longer(cols = `Montreal.Pre-Pandemic`:`Vancouver.Post-Restriction`,
                 names_to = c("city", "time_pt"), names_sep = "\\.", values_to = "r0") 
  
  dat_r0$time_pt <- factor(dat_r0$time_pt,
                           levels = c("Pre-Pandemic", "Pandemic", "Post-Restriction"))
  dat_r0.l$time_pt <- factor(dat_r0.l$time_pt,
                             levels = c("Pre-Pandemic", "Pandemic", "Post-Restriction"))
  dat_r0.u$time_pt <- factor(dat_r0.u$time_pt,
                             levels = c("Pre-Pandemic", "Pandemic", "Post-Restriction"))
  
  df_list <- list(dat_r0, dat_r0.l, dat_r0.u)      
  dat_r0_fu <- df_list %>% reduce(full_join, by=c('SAR', "city", "time_pt")) %>% rename(r0 = r0.x,
                                                                                        r0.l = r0.y,
                                                                                        r0.u = r0)
  # one time point per panel
  
  ## create a df with R0 and city name
  
  dat_r0_draw <- dat_r0_fu %>% subset(time_pt != "Pandemic")
  dat_r0_draw$city[dat_r0_draw$city == "Montreal"] <- "Montréal"
  
  # R0 from cases using 100-day growth period 
  df_case <- data.frame(city = c("Montréal", "Toronto", "Vancouver"),
                        r0_case = r0_case)
  R0_case <- df_case$r0_case
  
  # R0_sexual activity, when SAR = 1
  R0_sex_post <- dat_r0_fu[dat_r0_fu$SAR == 1 & dat_r0_fu$time_pt == "Post-Restriction", ]$r0
  R0_sex_post.lb <- dat_r0_fu[dat_r0_fu$SAR == 1 & dat_r0_fu$time_pt == "Post-Restriction", ]$r0.l
  R0_sex_post.ub <- dat_r0_fu[dat_r0_fu$SAR == 1 & dat_r0_fu$time_pt == "Post-Restriction", ]$r0.u
  
  SAR_est[[dataset_name]] <- mean(R0_case / R0_sex_post)
  SAR_est.lb[[dataset_name]] <- mean(R0_case / R0_sex_post.ub)
  SAR_est.ub[[dataset_name]] <- mean(R0_case / R0_sex_post.lb)
  
  R0_sex_pre <- dat_r0_fu[dat_r0_fu$SAR == 1 & dat_r0_fu$time_pt == "Pre-Pandemic", ]$r0
  R0_sex_pre.lb <- dat_r0_fu[dat_r0_fu$SAR == 1 & dat_r0_fu$time_pt == "Pre-Pandemic", ]$r0.l
  R0_sex_pre.ub <- dat_r0_fu[dat_r0_fu$SAR == 1 & dat_r0_fu$time_pt == "Pre-Pandemic", ]$r0.u

  df_r0_pre <- data.frame(city = c("Montréal", "Toronto", "Vancouver"),
                          r0_est = SAR_est[[dataset_name]] * R0_sex_pre, # estimated R0 if sexual activities were at the pre-pandemic level
                          r0_pre_post_prop = SAR_est[[dataset_name]] * R0_sex_pre / R0_case,
                          r0_pre_post_prop.lb = SAR_est.lb[[dataset_name]] * R0_sex_pre.lb / R0_case,
                          r0_pre_post_prop.ub = SAR_est.ub[[dataset_name]] * R0_sex_pre.ub / R0_case)
  
  ## add case-based R0 estimate to df with sex network-based R0
  # note: this should replicate the same R0 for each row by city
  dat_r0_draw <- left_join(dat_r0_draw, df_case, by = "city")
  dat_r0_draw <- left_join(dat_r0_draw, df_r0_pre, by ="city")
  
  # get R0_ntwrk with closest value
  dat_r0_draw$r0_case_diff <- abs(dat_r0_draw$r0_case - dat_r0_draw$r0)
  dat_r0_draw$r0_est_diff <- abs(dat_r0_draw$r0_est - dat_r0_draw$r0)
  
  dat_r0_draw <- dat_r0_draw %>% 
    group_by(city, time_pt) %>% 
    mutate(is_closest_case = (r0_case_diff == min(r0_case_diff)),
           is_closest_est = (r0_est_diff == min(r0_est_diff))) %>% 
    ungroup()
  
  df_draw[[dataset_name]] <- dat_r0_draw
  
  # make df to plot the points
  dat_points <- subset(dat_r0_draw, (is_closest_case & time_pt == "Post-Restriction") | (is_closest_est & time_pt == "Pre-Pandemic") )
  df_point[[dataset_name]] <- dat_points
}

# estimate SAR and R0 given pre-pandemic sexual activity level -------

SAR_df <- do.call(cbind.data.frame, SAR_est) %>%
  pivot_longer(cols = all_of(DATASET),
               names_to = "dataset",
               values_to = "value")

draw_df <- do.call(rbind.data.frame, df_draw) %>% 
  mutate(dataset = as.vector(sapply(DATASET, function (x) rep(x, 606)))) %>% 
  mutate(Dataset = case_when(dataset == "data_full" ~ "Full",
                             dataset == "data_restr" ~ "Restricted",
                             dataset == "data_age" ~ "Age-Standardized",
                             dataset == "data_anal" ~ "Anal Sexual Partners"),
                             Dataset = factor(Dataset, levels = c("Full", "Restricted", "Age-Standardized", "Anal Sexual Partners")))

point_df <- do.call(rbind.data.frame, df_point) %>% 
  mutate(dataset = as.vector(sapply(DATASET, function (x) rep(x, 6)))) %>% 
  mutate(Dataset = case_when(dataset == "data_full" ~ "Full",
                                                  dataset == "data_restr" ~ "Restricted",
                                                  dataset == "data_age" ~ "Age-Standardized",
                                                  dataset == "data_anal" ~ "Anal Sexual Partners"),
                              Dataset = factor(Dataset, levels = c("Full", "Restricted", "Age-Standardized", "Anal Sexual Partners")))

# plot R0 estimate -----
# comparison across diff datasets 
# STILL EDITING
# TODO: MOVE INTO SEPARATE SCRIPT
p_r0_compare <- ggplot(draw_df, aes(x = SAR, y = r0, col = city)) +
  # line to sh where R0=1 is
  geom_hline(yintercept = 1, linetype = "dotted") +
  # draw the 
  geom_line() +
  geom_point(data = point_df,
             aes(x = SAR, y = r0, col = city),
             size = 1) +
  facet_grid(Dataset ~ time_pt) +
  geom_segment(data = point_df,
               aes(x = SAR, xend = SAR, y = 0, yend = r0, col = city),
               linetype = "dashed",
               size = 0.6) +
  geom_segment(data = point_df,
               aes(x = 0, xend = SAR, y = r0, yend = r0, col = city),
               linetype = "dashed",
               size = 0.6) +
  coord_cartesian(ylim = c(0, 3.5)) +
  geom_ribbon(aes(ymin = r0.l, ymax = r0.u, fill = city), alpha = 0.2, colour = NA) +
  scale_colour_viridis_d(option = "C", end = 0.8) +
  scale_fill_viridis_d(option = "C", end = 0.8) +
  xlab("Secondary Attack Rate (SAR)") + 
  ylab(expression(italic('R'[0])))  +
  theme(plot.title.position = "plot",
        plot.caption.position = "plot",
        plot.caption = element_markdown(hjust = 0, size = 5.5))

png("./figures-3cities/negbin-ipcw/figS_r0.png",
    width = 18, height = 18, units = "cm", res = 600)

dev.off()

# TODO: harmonize sizing with fig1
# only using full dataset
draw_df_main <- filter(draw_df, dataset == "data_full" & time_pt == "Post-Restriction")
p_r0_main <- ggplot(draw_df_main, aes(x = SAR, y = r0, col = city)) +
  # line to show where R0=1 is
  geom_hline(yintercept = 1, linetype = "dotted") +
  
  # draw the R0 NGM estimate
  geom_ribbon(aes(ymin = r0.l, ymax = r0.u, fill = city), alpha = 0.2, colour = NA) +
  geom_line() +
  # draw the R0 estimates from growth rate
  geom_point(data = filter(point_df, dataset == "data_full" & time_pt == "Post-Restriction"),
             aes(x = SAR, y = r0, col = city),
             size = 1) +
  geom_segment(data = filter(point_df, dataset == "data_full" & time_pt == "Post-Restriction"),
               aes(x = SAR, xend = SAR, y = 0, yend = r0, col = city),
               linetype = "dashed",
               size = 0.6) +
  geom_segment(data = filter(point_df, dataset == "data_full" & time_pt == "Post-Restriction"),
               aes(x = 0, xend = SAR, y = r0, yend = r0, col = city),
               linetype = "dashed",
               size = 0.6) +
  coord_cartesian(ylim = c(0, 3.5)) +
  
  # facet_wrap( ~ time_pt) +
  
  # aesthetics
  scale_colour_viridis_d(option = "C", end = 0.8) +
  scale_fill_viridis_d(option = "C", end = 0.8) +
  
  # axis labels
  labs(title = "A)",
       x = "Secondary Attack Rate", 
       y = expression(italic("R"[0])*" (estimated)"),
       # y = bquote(atop(R["0"])),
       col = "City", fill = "City") +
  theme_default

draw_r0_prepand <- filter(draw_df, dataset == "data_full" & time_pt == "Pre-Pandemic")

# plot only the SAR that was estimated
draw_r0_prepand <- draw_r0_prepand %>% 
  filter(SAR == round(SAR_df$value, 2))

p_r0_rel <- ggplot(draw_r0_prepand, aes(x = city, y = r0, col = city)) +
  # line to show where R0=1 is
  geom_hline(yintercept = 1, linetype = "dotted") +
  
  # draw the R0 NGM estimate
  geom_pointrange(aes(ymin = r0.l, ymax = r0.u), fatten = 0.5) +
  coord_cartesian(ylim = c(0, 3.5)) +
  
  # aesthetics
  scale_colour_viridis_d(option = "C", end = 0.8) +
  
  # axis labels
  labs(title = "B)", x = "City",
       y = expression(italic("R"[0])*" (projected)")) +
       # y = bquote(atop(R["0"]*" assuming", "pre-pandemic sexual activity"))) +
  theme(legend.position = "none") +
  theme_default

# get legend
plt_legend <- cowplot::get_legend(p_r0_main)

png("./figures-3cities/negbin-ipcw/fig1_r0.png",
    width = 12, height = 5.5, units = "cm", res = 600)
grid.arrange(
  p_r0_main + theme(legend.position = "none"),
  p_r0_rel,
  plt_legend,
  nrow = 1, widths = unit(c(5, 5, 2), "cm")
)
dev.off()
