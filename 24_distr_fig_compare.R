# dataset and library ----
library(tidyverse)
library(data.table)
library(gridExtra)
library(ggtext)
source("./src/utils_regression_3cities.R")
source("./src/plot.R")

theme_set(theme_bw())

## which variable to use as outcome
outcome_var <- "nb_part_ttl" 
file_suff <- case_when(outcome_var == "nb_part_ttl" ~ "p6m_all",
                       outcome_var == "nb_part_anal" ~ "p6m_anal")
fig_path <- "./figures-3cities/negbin-ipcw"

# load main fitted results
data_ipcw <- read_csv("./output-3cities/negbin-ipcw/full_fitted_cdf_wt_p6m_all.csv")
# TODO: MOVE THESE INTO SEPARATE SCRIPT FOR SUPPLEMENTARY FIGS
data_restr <- read_csv("./output-3cities/negbin-res/full_fitted_cdf_wt_p6m_all.csv")
data_age <- read_csv("./output-3cities/negbin-age/full_fitted_cdf_wt_p6m_all.csv")
data_anal <- read_csv("./output-3cities/negbin-ipcw/full_fitted_cdf_wt_p6m_anal.csv")

# # fix  labels (check that order matches) (verified)
# unique(data_restr$data_pt)
# unique(data_ipcw$data_pt)
# unique(data_age$data_pt)
# unique(data_anal$data_pt)

data_fit_cdf <- bind_rows(
  mutate(data_ipcw, dataset = "Full"),
  mutate(data_restr, dataset = "Restriction"),
  mutate(data_restr, dataset = "Age-Standardized"),
  mutate(data_anal, dataset = "Anal-Sex-Partner")
) %>%
  mutate(dataset = factor(dataset, level = c("Full", "Restriction", "Age-Standardized", "Anal-Sex-Partner")))
  # mutate(across(mean:cr.i_upp, prop_to_perc))

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


# plot full model --------------------------------------------------------------------

cdf_wt_by_city <- data_fit_cdf %>% subset(dataset == "Full")
## theme for plot
# TODO move into setup file
theme_default <- theme(axis.title = element_text(size = 9),
                       axis.text = element_text(size = 7),
                       strip.text = element_text(size = 7),
                       legend.title = element_text(size = 9),
                       legend.text = element_text(size = 7),
                       plot.title = element_text(size = 10),
                       plot.title.position = "plot",
                       plot.caption.position = "plot")

# one facet per timepoint
# TODO: move x label into plot function
p_wt_by_time <- plot_cdf_fit(cdf_wt_by_city, col_var = "city",
                             p_title = NULL,
                             outcome_type = gsub("p6m_", "", file_suff)) +
  facet_wrap(~time_pt) +
  labs(title = "A)",
       x = "Number of reported sexual partners in the P6M", y = NULL, 
       col = "City", fill = "City")

# one facet per city
p_wt_by_city <- plot_cdf_fit(cdf_wt_by_city, col_var = "time_pt",
                             p_title = NULL,
                             outcome_type = gsub("p6m_", "", file_suff)) +
  facet_wrap(~city) +
  labs(title = "A)",
       x = "Number of reported sexual partners in the P6M", y = NULL, 
       col = "Time period", fill = "Time period")

# log scale
p_out_by_time <- p_wt_by_time +
  scale_x_continuous(trans = scales::pseudo_log_trans(),
                     breaks = c(0, 10, 100, 300))
  # labs(caption = "Displayed in pseudo-log transformation") 

p_out_by_city <- p_wt_by_city +
  scale_x_continuous(trans = scales::pseudo_log_trans(),
                     breaks = c(0, 10, 100, 300))
  # labs(caption = "Displayed in pseudo-log transformation")


# theme for text size
# theme_txt <- theme(axis.title = element_text(size = 9),
                   # axis.text = element_text(size = 7),
                   # strip.text = element_text(size = 7),
                   # legend.text = element_text(size = 9),
                   # plot.title = element_text(size = 7))

# y label
yaxis <- grid::textGrob(
  "Proportion reporting at least x partners",
  gp = grid::gpar(fontsize = 9, lineheight = 0.8), rot = 90
)

# select cut-offs to show and turn into proportion
# TODO: CREATE SINGLE OUTPUT FOR COLOUR LEGEND
prop_at_least_x <- cdf_wt_by_city %>% 
  filter(y_pred %in% c(25, 50, 100, 150)) %>% 
  arrange(y_pred, city, time_pt) 

p_out_by_time_sel <- plot_cdf_sel(filter(prop_at_least_x), "city", file_suff) +
  facet_wrap(~ time_pt) +
  coord_cartesian(ylim = c(0, .20)) +
  labs(title = "B)",
       x = "Number of reported sexual partners in the P6M", y = NULL,
       col = "City")

p_out_by_city_sel <- plot_cdf_sel(filter(prop_at_least_x), "time_pt", file_suff) +
  facet_wrap(~ city) +
  coord_cartesian(ylim = c(0, .20)) +
  labs(title = "B)",
       x = "Number of reported sexual partners in the P6M", y = NULL,
       col = "Time period")


png("./figures-3cities/negbin-ipcw/figS_full_model_cdf_strat_by_time.png",
    width = 15, height = 11, units = "cm", res = 600)
grid.arrange(
  p_out_by_time + theme_default,
  p_out_by_time_sel + theme_default,
  yaxis, layout_matrix = matrix(c(3, 3, 1, 2), nrow = 2, byrow = F),
  widths = unit(c(0.4, 14.6), "cm")
)
dev.off()

png("./figures-3cities/negbin-ipcw/figS_full_model_cdf_strat_by_city.png",
    width = 15, height = 11, units = "cm", res = 600)
grid.arrange(
  p_out_by_city + labs(y = NULL, title = "A)") + theme_default,
  p_out_by_city_sel + labs(y = NULL, title = "B)") + theme_default,
  yaxis, layout_matrix = matrix(c(3, 3, 1, 2), nrow = 2, byrow = F),
  widths = unit(c(0.4, 14.6), "cm")
)
dev.off()

# TODO: MOVE INTO SEPARATE SCRIPT
# plot comparison between models --------------------------------------------------------------------

png("./figures-3cities/negbin-ipcw/figS_mega_comp.png",
    width = 15, height = 12, units = "cm", res = 600)

plot_cdf_fit(data_fit_cdf, col_var = "dataset", outcome_type = "all") +
  facet_grid(time_pt~city) +
  scale_x_continuous(trans = scales::pseudo_log_trans(),
                     breaks = c(0, 10, 25, 50, 100, 300)) +
  scale_colour_viridis_d(option = "F", end = 0.9) +
  scale_fill_viridis_d(option = "F", end = 0.9) + 
  labs(col = "dataset", fill = "dataset")  + 
  scale_y_continuous(trans='log10',
                     breaks = function(x) {
                       brks <- extended_breaks(Q = c(1, 5))(log10(x))
                       10^(brks[brks %% 1 == 0])
                     },
                     labels = math_format(format = log10)) +
  theme_default +
  labs(y = "Proportion reporting at\nleast x partners (%)") +
       # caption = "**Figure S6. Comparison of RDS-IPC weighted complementary cumulative degree distribution of sexual partners in P6M <br> using full, restricted, age-standardized dataset, as well as using anal sex partner as the study outcome.** The y-axis is displayed <br> in log10 scale.") +
  guides(col = guide_legend(title = "Dataset"),
         fill = guide_legend(title = "Dataset"))

dev.off()

# plot select timepoints
png("./figures-3cities/negbin-ipcw/figS_mega_comp_sel.png",
    width = 30, height = 18, units = "cm", res = 600)

plot_cdf_sel(filter(data_fit_cdf, y_pred %in% c(10, 25, 50, 100, 150)),
             col_var = "dataset", file_suff = "all") +
  facet_grid(time_pt ~ city) +
  scale_colour_viridis_d(option = "F", end = 0.9) +
  scale_fill_viridis_d(option = "F", end = 0.9) +
  theme_default +
  labs(y = "Proportion reporting at\nleast x partners (%)") +
       # caption = "**Figure S7. Comparison of RDS-IPC weighted complementary cumulative degree distribution of sexual partners in P6M using full, restricted, age-standardized dataset, as well as using anal sex partner as the study outcome.**") +
  guides(col = guide_legend(title = "Dataset"),
         fill = guide_legend(title = "Dataset"))

dev.off()
