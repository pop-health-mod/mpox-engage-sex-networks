# Libraries and data -----
library(tidyverse)
library(gridExtra)

theme_set(theme_bw())

source("./src/plot.R")
source("./src/plot_config.R")

# whether to use log-log plots for the sensitivity analysis figures
# or use natural scale for Y axis (X axis always in log scale)
plot_sens_loglog <- TRUE

# main analyses
data_ipcw <- read_csv("./out/fitted-distributions/cdf_weighted_all_partn.csv")

# sensitivity analyses
data_restr <- read_csv("./out/fitted-distr-sens-restrict/cdf_weighted_all_partn.csv")
data_stand <- read_csv("./out/fitted-distr-sens-standard/cdf_weighted_all_partn.csv")
data_anal <- read_csv("../mpx-engage-params/output-3cities/negbin-ipcw/full_fitted_cdf_wt_p6m_anal.csv")

data_anal$time_pt <- gsub("Restriction", "Restrictions", data_anal$time_pt)

# put into single data frame 
data_fit_cdf <- bind_rows(
  mutate(data_ipcw, dataset = "Main"),
  mutate(data_restr, dataset = "Restriction"),
  mutate(data_stand, dataset = "Standardization"),
  mutate(data_anal, dataset = "Anal")
)

# create city and time point variables
data_fit_cdf <- data_fit_cdf %>% 
  mutate(city = case_when(grepl("mtl|Montreal", data_pt) ~ "Montréal",
                          grepl("trt|Toronto", data_pt) ~ "Toronto",
                          grepl("van|Vancouver", data_pt) ~ "Vancouver"),
         time_pt = case_when(grepl("baseline|Pre-Pandemic", data_pt) ~ "Pre-Pandemic",
                             grepl("pandemic|Pandemic", data_pt) ~ "Pandemic",
                             grepl("mpox|Post-Restriction", data_pt) ~ "Post-Restrictions"),
         time_pt = factor(time_pt,
                          levels = c("Pre-Pandemic", "Pandemic", "Post-Restrictions"))
  )

## Transformation for log-log plots ----
# if plotting in log-log, pre-transform y values into log
if(plot_sens_loglog){
  # transform fitted CDF
  data_fit_cdf <- data_fit_cdf %>% 
    mutate(
      mean = log(mean),
      cr.i_low = log(cr.i_low),
      cr.i_upp = log(cr.i_upp)
    )
  
  # breaks for the Y axis
  y_log_breaks <- c(1, .1, .01, 1e-4, 1e-6)
  
  # min and max to use in the Y axis
  y_minmax <- log(c(10^-6.5, 1))
}

# Plot: restriction ----
# scales::show_col(scales::viridis_pal(option = "B", end = 0.8)(6))
col_pal <- scales::viridis_pal(option = "B", end = 0.8)(6)[c(3, 6)]

p_restr <- plot_cdf_fit(
  filter(data_fit_cdf, dataset %in% c("Main", "Restriction")),
  col_var = "dataset", outcome_type = "all"
) +
  facet_grid(time_pt ~ city) +
  scale_x_continuous(trans = scales::pseudo_log_trans(),
                     breaks = c(0, 10, 100, 300)) +
  
  scale_fill_manual(values = col_pal) +
  scale_colour_manual(values = col_pal) +

  labs(col = "Analysis", fill = "Analysis") +
  
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 12),    # axis titles size
    axis.text = element_text(size = 10),     # axis text size
  )

if(plot_sens_loglog){
  p_restr <- p_restr +
    coord_cartesian(ylim = y_minmax) +
    scale_y_continuous(breaks = log(y_log_breaks),
                       labels = y_log_breaks)
}

p_restr

ggsave("./fig/fig_S3_cdf_main_vs_restriction.png",
       device = "png",
       width = 15, height = 13, units = "cm", dpi = 600)

# Plot: standardization ----
p_stand <- plot_cdf_fit(
  filter(data_fit_cdf, dataset %in% c("Main", "Standardization")),
  col_var = "dataset", outcome_type = "all"
) +
  facet_grid(time_pt ~ city) +
  scale_x_continuous(trans = scales::pseudo_log_trans(),
                     breaks = c(0, 10, 100, 300)) +
  
  scale_fill_manual(values = col_pal) +
  scale_colour_manual(values = col_pal) +
  
  labs(col = "Analysis", fill = "Analysis") +
  
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 12),    # axis titles size
    axis.text = element_text(size = 10),     # axis text size
  )

if(plot_sens_loglog){
  p_stand <- p_stand +
    coord_cartesian(ylim = y_minmax) +
    scale_y_continuous(breaks = log(y_log_breaks),
                       labels = y_log_breaks)
}

p_stand

ggsave("./fig/fig_S4_cdf_main_vs_standardization.png",
       device = "png",
       width = 15, height = 13, units = "cm", dpi = 600)

# Plot: anal partners ----
# for anal partners outcome, need to relabel the 'dataset' variable
data_fit_anal <- data_fit_cdf %>% 
  filter(dataset %in% c("Main", "Anal")) %>% 
  mutate(dataset = factor(dataset,
                          levels = c("Main", "Anal"),
                          labels = c("All partners", "Anal partners")))

p_anal <- plot_cdf_fit(
  data_fit_anal,
  col_var = "dataset", outcome_type = "all"
) +
  facet_grid(time_pt ~ city) +
  scale_x_continuous(trans = scales::pseudo_log_trans(),
                     breaks = c(0, 10, 100, 300)) +
  
  scale_fill_manual(values = col_pal) +
  scale_colour_manual(values = col_pal) +
  
  labs(col = "Outcome", fill = "Outcome") +
  
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 12),    # axis titles size
    axis.text = element_text(size = 10),     # axis text size
  )

if(plot_sens_loglog){
  p_anal <- p_anal +
    coord_cartesian(ylim = y_minmax) +
    scale_y_continuous(breaks = log(y_log_breaks),
                       labels = y_log_breaks)
}

p_anal

ggsave("./fig/fig_S5_cdf_main_all_vs_anal.png",
       device = "png",
       width = 15, height = 13, units = "cm", dpi = 600)
