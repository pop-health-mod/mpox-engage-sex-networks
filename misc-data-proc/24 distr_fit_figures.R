# Library and data----
library(tidyverse)

source("./src/plot.R")
source("./src/plot_config.R")
theme_set(theme_bw())

## define output paths
fig_path <- "./misc-data-proc/fig-cm"

### load data
## empirical
dat_distr_empir_cty <- read.csv("./misc-data-proc/outputs-cm/engage_van_empiric_nb_partn_p6m.csv")
dat_distr_empir_cty <- subset(dat_distr_empir_cty, n > 0)

# compute empirical CDF
dat_distr_empir_cty <- dat_distr_empir_cty %>% 
  group_by(time_pt) %>% 
  arrange(-x) %>% 
  mutate(cdf = cumsum(p), cdf.rds = cumsum(p.rds)) %>% 
  arrange(x)
dat_distr_empir_cty <- ungroup(dat_distr_empir_cty)

## fitted
dat_distr_fit_cty <- read.csv("./misc-data-proc/outputs-cm/engage_van_fitted_cdf_nb_partn_p6m.csv")
# dat_distr_fit_cty <- subset(dat_distr_fit_cty, time_pt == "Pre-Pandemic")
# dat_distr_fit_cty <- mutate(dat_distr_fit_cty, city = "van", .before = 1)

# Plot city-wide ----
# whether to use log scale on both x and y axis
plot_loglog <- TRUE

# if plotting in log-log, pre-transform y values into log
if(plot_loglog){
  # transform fitted CDF
  dat_distr_fit_cty <- dat_distr_fit_cty %>% 
    mutate(
      mean = log(mean),
      cr.i_low = log(cr.i_low),
      cr.i_upp = log(cr.i_upp)
    )
  
  # transform observed data
  dat_distr_empir_cty <- dat_distr_empir_cty %>% 
    mutate(
      cdf.rds = log(cdf.rds)
    )
}

# breaks, min and max to use in Y axis
y_log_breaks <- c(1, .1, .01, 1e-4, 1e-6)

if(plot_loglog){
  y_minmax <- log(c(10^-6.5, 1))
} else {
  y_minmax <- c(0, 1)
}

# plot
p_cdf <- plot_cdf_fit(dat_distr_fit_cty, "time_pt", # plots fitted data
                      p_title = NULL,
                      outcome_type = "all") +
  
  # plot  observed data
  geom_step(data = dat_distr_empir_cty,
            aes(x = x, y = cdf.rds, col = time_pt),
            linewidth = 0.3) +
  geom_point(data = dat_distr_empir_cty,
             aes(x = x, y = cdf.rds, col = time_pt),
             size = 0.7) +
  
  # format
  coord_cartesian(ylim = y_minmax) +
  scale_x_continuous(trans = scales::pseudo_log_trans(),
                     breaks = c(0, 10, 100, 300)) +
  scale_fill_viridis_d(option = "B", begin = 0.4, end = 0.8) +
  scale_colour_viridis_d(option = "B", begin = 0.4, end = 0.8) +
  
  facet_grid(~time_pt) +
  
  # fix proportions of plot elements
  labs(x = "Number of reported sexual partners\nin the P6M",
       y = "Proportion reporting\n>=x partners") +
  theme(
    legend.position = "none",                # not needed, only Vancouver
    
    axis.title = element_text(size = 12),    # axis titles size
    axis.text = element_text(size = 10),     # axis text size
  )

if(plot_loglog){
  p_cdf <- p_cdf +
    scale_y_continuous(breaks = log(y_log_breaks),
                       labels = y_log_breaks)
}

p_cdf

ggsave("./misc-data-proc/fig-cm/nb_partn_distr_fit_check.png",
       device = "png",
       width = 21, height = 8, units = "cm", dpi = 600)


