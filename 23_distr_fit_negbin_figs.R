# Library and data----
library(tidyverse)
library(gridExtra)
library(ggtext)

source("./src/plot.R")
source("./src/plot_config.R")

theme_set(theme_bw())

## load main fitted results
data_ipcw <- read_csv("./out/fitted-distributions/cdf_weighted_all_partn.csv")
data_ipcw$time_pt[is.na(data_ipcw$time_pt)] <- "Post-Restrictions" ## RERUN FITTING AND DELETE
data_ipcw$city <- gsub("s", "", data_ipcw$city)

# format city and time point variables (for plotting)
data_ipcw <- data_ipcw %>% 
  mutate(
    city = factor(city,
                  levels = c("Montreal", "Toronto", "Vancouver"),
                  labels = c("Montréal", "Toronto", "Vancouver")),
    time_pt = factor(time_pt,
                     levels = c("Pre-Pandemic", "Pandemic", "Post-Restrictions"))
  )

## load Engage data
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
  mutate(city = recode_factor(city, mtl = "Montréal", trt = "Toronto", van = "Vancouver"))

# Plot full CDF and select points (main figures) ----
## plots of full CDF
# both plots use log scale for x axis
# one facet per timepoint
p_by_time <- plot_cdf_fit(data_ipcw, col_var = "city",
                          p_title = NULL,
                          outcome_type = gsub("p6m_", "", file_suff)) +
  facet_wrap(~time_pt) +
  scale_x_continuous(trans = scales::pseudo_log_trans(),
                     breaks = c(0, 10, 100, 300)) +
  
  labs(title = "A)", y = NULL, 
       col = "City", fill = "City")

# one facet per city
p_by_city <- plot_cdf_fit(data_ipcw, col_var = "time_pt",
                          p_title = NULL,
                          outcome_type = gsub("p6m_", "", file_suff)) +
  facet_wrap(~city) +
  scale_x_continuous(trans = scales::pseudo_log_trans(),
                     breaks = c(0, 10, 100, 300)) +
  
  scale_fill_viridis_d(option = "B", end = 0.8) +
  scale_colour_viridis_d(option = "B", end = 0.8) +
  labs(title = "C)", y = NULL, 
       col = "Time period", fill = "Time period")

# y label
yaxis <- grid::textGrob(
  "Proportion reporting at least x partners",
  gp = grid::gpar(fontsize = 9, lineheight = 0.8), rot = 90
)

## select cut-offs to show and turn into proportion
# TODO harmonize sizing of legends between time points and city
prop_at_least_x <- data_ipcw %>% 
  filter(y_pred %in% c(25, 50, 100, 150)) %>% 
  arrange(y_pred, city, time_pt) 

# one facet per timepoint
p_by_time_sel <- plot_cdf_sel(filter(prop_at_least_x), "city", file_suff) +
  facet_wrap(~ time_pt) +
  coord_cartesian(ylim = c(0, .20)) +
  labs(title = "B)",
       x = "Number of reported sexual partners in the P6M", y = NULL,
       col = "City")

# plot full figure
png("./fig/fig_1_cdf_main_model.png",
    width = 15, height = 16, units = "cm", res = 600)
grid.arrange(
  p_by_time + theme_cdf,
  p_by_time_sel + theme_cdf,
  p_by_city + theme_cdf,
  yaxis,
  layout_matrix = matrix(c(4, 4, 4, 1, 2, 3), nrow = 3, byrow = F),
  widths = unit(c(0.4, 14.6), "cm")
)
dev.off()

# Check goodness of fit (CDF from NB vs weighted observed) ----
## compute frequency each one was reported
data_cdf_obs <- data_3cities %>% 
  group_by(city, time_pt, nb_part_ttl) %>% 
  summarize(n = n(), wt = sum(ipw_rds), .groups = "drop_last")

# compute frequencies and check that population and total density are OK
data_cdf_obs <- data_cdf_obs %>% 
  mutate(prop_raw = n / sum(n),
         prop_wt = wt / sum(wt))

# verify sample size and proportions
tapply(data_cdf_obs$n, paste(data_cdf_obs$city, data_cdf_obs$time_pt), sum)
tapply(data_cdf_obs$prop_raw, paste(data_cdf_obs$city, data_cdf_obs$time_pt), sum)
tapply(data_cdf_obs$prop_wt, paste(data_cdf_obs$city, data_cdf_obs$time_pt), sum)

# compute inverse CDF
data_cdf_obs <- data_cdf_obs %>% 
  arrange(city, time_pt, -nb_part_ttl) %>% 
  mutate(inv_cdf_raw = cumsum(prop_raw),
         inv_cdf_wt = cumsum(prop_wt)) %>% 
  arrange(city, time_pt, nb_part_ttl)

data_cdf_obs <- data_cdf_obs %>% ungroup() # note: needs to be grouped up to now for sums to work as intended

### plot
plot_cdf_fit(data_ipcw, "city", # plots fitted data
             p_title = NULL,
             outcome_type = "all") +
  
  # plot  observed data
  geom_step(data = data_cdf_obs,
            aes(x = nb_part_ttl, y = inv_cdf_wt, col = city),
            linewidth = 0.3) +
  geom_point(data = data_cdf_obs,
             aes(x = nb_part_ttl, y = inv_cdf_wt, col = city),
             size = 0.7) +
  
  facet_grid(time_pt ~ city) +
  
  # formt
  coord_cartesian(ylim = c(0, 1)) +
  scale_x_continuous(trans = scales::pseudo_log_trans(),
                     breaks = c(0, 10, 100, 300)) +
  
  labs(col = NULL) +
  
  # fix proportions of plot elements
  # adjust bold face, and caption position
  theme(
    legend.position = "none",                               # not needed because of wrapping
    
    axis.title = element_text(size = 12),    # axis titles size
    axis.text = element_text(size = 10),     # axis text size
  )

ggsave("./fig/fig_S1_cdf_model_fit_to_data.png",
       device = "png",
       width = 15, height = 11, units = "cm", dpi = 600)
