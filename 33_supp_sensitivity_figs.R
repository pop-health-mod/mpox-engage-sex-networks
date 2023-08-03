# Libraries and data -----
library(tidyverse)
library(gridExtra)

theme_set(theme_bw())

source("./src/plot.R")
source("./src/plot_config.R")

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
) %>%
  mutate(dataset = factor(dataset, level = c("Main", "Restriction", "Standardization", "Anal")))

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

# Plot: restriction ----
# scales::show_col(scales::viridis_pal(option = "B", end = 0.8)(6))
col_pal <- scales::viridis_pal(option = "B", end = 0.8)(6)[c(3, 6)]

plot_cdf_fit(
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

ggsave("./fig/fig_S3_cdf_main_vs_restriction.png",
       device = "png",
       width = 15, height = 13, units = "cm", dpi = 600)

# Plot: standardization ----
plot_cdf_fit(
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

ggsave("./fig/fig_S4_cdf_main_vs_standardization.png",
       device = "png",
       width = 15, height = 13, units = "cm", dpi = 600)

# Plot: anal partners ----
# TODO write properly
plot_cdf_fit(
  filter(data_fit_cdf, dataset %in% c("Main", "Anal")) %>% 
    mutate(dataset = factor(dataset,
                            levels = c("Main", "Anal"),
                            labels = c("All partners", "Anal partners"))),
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

ggsave("./fig/fig_S5_cdf_main_all_vs_anal.png",
       device = "png",
       width = 15, height = 13, units = "cm", dpi = 600)
