# dataset and library ----
library(tidyverse)
library(data.table)
library(gridExtra)
library(ggtext)
source("./src/utils_regression_3cities.R")
source("./src/plot.R")

theme_set(theme_bw())

## which variable to use as outcome
file_suff <- "p6m_all"
fig_path <- "./figures-3cities/negbin-ipcw"

# load main fitted results
data_fit <- read_csv("./output-3cities/negbin-ipcw/full_fitted_cdf_wt_p6m_all.csv")

# create city and time point variables
data_fit <- data_fit %>% 
  mutate(city = case_when(grepl("mtl|Montreal", data_pt) ~ "Montréal",
                          grepl("trt|Toronto", data_pt) ~ "Toronto",
                          grepl("van|Vancouver", data_pt) ~ "Vancouver"),
         time_pt = case_when(grepl("baseline|Pre-Pandemic", data_pt) ~ "Pre-Pandemic",
                             grepl("pandemic|Pandemic", data_pt) ~ "Pandemic",
                             grepl("mpox|Post-Restriction", data_pt) ~ "Post-Restrictions"),
         time_pt = factor(time_pt,
                          levels = c("Pre-Pandemic", "Pandemic", "Post-Restrictions"))
  )

# load raw data
ipcw_pre <- read.csv("./data-3cities-feb-2023/pre_ipcw_3cities.csv")
ipcw_pand <- read.csv("./data-3cities-feb-2023/pand_ipcw_3cities.csv")
ipcw_post <- read.csv("./data-3cities-feb-2023/post_ipcw_3cities.csv")

data_3cities <- ipcw_pre %>% 
  bind_rows(ipcw_pand, ipcw_post) %>%
  mutate(time_pt = factor(time_pt, 
                          levels = c("Pre-Pandemic", "Pandemic", "Post-Restriction"),
                          labels = c("Pre-Pandemic", "Pandemic", "Post-Restrictions"))) %>%
  mutate(city = recode_factor(city, mtl = "Montréal", trt = "Toronto", van = "Vancouver")) %>% 
  as_tibble()

# Plot empirical and fitted CDF ----
# compare CDF negative binomial fit vs observed
## frequency each one was reported
freq_nb_partn <- data_3cities %>% 
  group_by(city, time_pt, nb_part_ttl) %>% 
  summarize(n = n(), wt = sum(ipw_rds), .groups = "drop_last")

# compute frequencies and check that population and total density are OK
freq_nb_partn <- freq_nb_partn %>% 
  mutate(prop_raw = n / sum(n),
         prop_wt = wt / sum(wt))

# verify sample size and proportions
tapply(freq_nb_partn$n, paste(freq_nb_partn$city, freq_nb_partn$time_pt), sum)
tapply(freq_nb_partn$prop_raw, paste(freq_nb_partn$city, freq_nb_partn$time_pt), sum)
tapply(freq_nb_partn$prop_wt, paste(freq_nb_partn$city, freq_nb_partn$time_pt), sum)

# compute inverse CDF
freq_nb_partn <- freq_nb_partn %>% 
  arrange(city, time_pt, -nb_part_ttl) %>% 
  mutate(inv_cdf_raw = cumsum(prop_raw),
         inv_cdf_wt = cumsum(prop_wt)) %>% 
  arrange(city, time_pt, nb_part_ttl)

freq_nb_partn <- freq_nb_partn %>% ungroup()

### plot
plot_cdf_fit(data_fit, "city",
             p_title = NULL,
             outcome_type = "all") +
  
  # plot  observed data
  geom_step(data = freq_nb_partn,
            aes(x = nb_part_ttl, y = inv_cdf_wt, col = city)) +
  geom_point(data = freq_nb_partn,
             aes(x = nb_part_ttl, y = inv_cdf_wt, col = city),
             size = 0.7) +
  
  facet_grid(time_pt ~ city) +
  
  coord_cartesian(ylim = c(0, 1)) +
  scale_x_continuous(trans = scales::pseudo_log_trans(),
                     breaks = c(0, 10, 100, 300)) +
  
  labs(x = "Number of reported sexual partners in the P6M",
       y = "Proportion reporting at least x partners",
       col = NULL) +
  
  # fix proportions of plot elements
  # adjust bold face, and caption position
  theme(
    legend.position = "none",                               # not needed because of wrapping
    
    axis.title = element_text(size = 12),    # axis titles size
    axis.text = element_text(size = 10),     # axis text size
  )

ggsave("./figures-3cities/negbin-ipcw/figS_model_fit_to_data.png",
       device = "png",
       width = 15, height = 11, units = "cm", dpi = 600)

