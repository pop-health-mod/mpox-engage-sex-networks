# Libraries and data ----
library(tidyverse)
library(gridExtra)
library(ggtext)

theme_set(theme_bw())

source("./src/utils_helper.R")
source("./src/plot.R")
source("./src/plot_config.R")

# Source R0 estimation ----
# only run if necessary
if(!exists("draw_df")){
  source("./24_distr_r0_ngm_case.R")
}

# Plot R0 estimates -----
# TODO: harmonize sizing with fig1
# only using full dataset
draw_df_main <- filter(draw_df, dataset == "data_full" & time_pt == "Post-Restrictions")
p_r0_main <- ggplot(draw_df_main, aes(x = SAR, y = r0, col = city)) +
  # line to show where R0=1 is
  geom_hline(yintercept = 1, linetype = "dotted") +
  
  # draw the R0 NGM estimate
  geom_ribbon(aes(ymin = r0.l, ymax = r0.u, fill = city), alpha = 0.2, colour = NA) +
  geom_line() +
  # draw the R0 estimates from growth rate
  geom_point(data = filter(point_df, dataset == "data_full" & time_pt == "Post-Restrictions"),
             aes(x = SAR, y = r0, col = city),
             size = 1) +
  geom_segment(data = filter(point_df, dataset == "data_full" & time_pt == "Post-Restrictions"),
               aes(x = SAR, xend = SAR, y = 0, yend = r0, col = city),
               linetype = "dashed",
               size = 0.6) +
  geom_segment(data = filter(point_df, dataset == "data_full" & time_pt == "Post-Restrictions"),
               aes(x = 0, xend = SAR, y = r0, yend = r0, col = city),
               linetype = "dashed",
               size = 0.6) +
  coord_cartesian(ylim = c(0, 4)) +
  
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
  theme_cdf

draw_r0_prepand <- filter(draw_df, dataset == "data_full" & time_pt == "Pre-Pandemic")

# plot only the SAR that was estimated
draw_r0_prepand <- draw_r0_prepand %>% 
  filter(SAR == round(SAR_df$value, 2))

p_r0_rel <- ggplot(draw_r0_prepand, aes(x = city, y = r0, col = city)) +
  # line to show where R0=1 is
  geom_hline(yintercept = 1, linetype = "dotted") +
  
  # draw the R0 NGM estimate
  geom_pointrange(aes(ymin = r0.l, ymax = r0.u), fatten = 0.5) +
  coord_cartesian(ylim = c(0, 4)) +
  
  # aesthetics
  scale_colour_viridis_d(option = "C", end = 0.8) +
  
  # axis labels
  labs(title = "B)", x = "City",
       y = expression(italic("R"[0])*" (projected)")) +
  # y = bquote(atop(R["0"]*" assuming", "pre-pandemic sexual activity"))) +
  theme(legend.position = "none") +
  theme_cdf

# get legend
plt_legend <- cowplot::get_legend(p_r0_main)

png("./fig/fig_2_r0_main.png",
    width = 12, height = 5.5, units = "cm", res = 600)
grid.arrange(
  p_r0_main + theme(legend.position = "none"),
  p_r0_rel,
  plt_legend,
  nrow = 1, widths = unit(c(5, 5, 2), "cm")
)
dev.off()
