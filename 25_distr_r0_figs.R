# Libraries and data ----
library(tidyverse)
library(gridExtra)
library(ggtext)

theme_set(theme_bw())

source("./src/utils_helper.R")
source("./src/plot.R")
source("./src/plot_config.R")

# load R0's (from NGM, projected, from cases)
df_r0_ngm <- read_csv("./out/r0-estim-ngm/R0_estim_NGM.csv")
df_r0_proj <- read_csv("./out/r0-estim-ngm/R0_projected.csv")

r0_case <- read_csv("./out/r0-estim-cases.csv")
# r0_case <- r0_case$x

# load SAR
df_sar <- read_csv("./out/r0-estim-ngm/SAR_estimates.csv")

# add accent to MTL
df_r0_ngm$city[df_r0_ngm$city == "Montreal"] <- "Montréal"
df_r0_proj$city[df_r0_proj$city == "Montreal"] <- "Montréal"
df_sar$city[df_sar$city == "Montreal"] <- "Montréal"

## RM ONCE RE-RAN
# names(df_r0_ngm)[1] <- "city"
# names(df_r0_proj)[1] <- "city"
# names(df_sar)[1] <- "city"
# 
# names(df_r0_ngm)[3] <- "mean_r0"
# names(df_r0_proj)[2] <- "mean_r0"

# df_r0_ngm <- df_r0_ngm %>% 
#   rename(r0.l = cr.i_low, 
#          r0.u = cr.i_upp)
# df_r0_proj <- df_r0_proj %>% 
#   rename(r0.l = cr.i_low, 
#          r0.u = cr.i_upp)

## Process data for plotting ----
# add corresponding R0 from NGM (point estimate) to SAR dataframe
df_sar_pt <- df_sar %>% 
  # select(city, SAR = mean) %>% 
  select(ngrp, city, SAR) %>%
  filter(city != "avg")
df_sar_pt

df_sar_pt$SAR <- round(df_sar_pt$SAR, 2)

df_sar_pt <- left_join(df_sar_pt, 
                       df_r0_ngm, 
                       by = c("ngrp", "city", "SAR")) %>% 
  select(ngrp, city, SAR, r0 = mean_r0)

# Plot R0 estimates -----
# TODO: harmonize sizing with fig1, make function for plotting

p_r0_main <- ggplot(filter(df_r0_ngm, ngrp == "grp100"), 
                    aes(x = SAR, y = mean_r0, col = city)) +
  # line to show where R0=1 is
  geom_hline(yintercept = 1, linetype = "dotted") +
  
  # draw the R0 NGM estimate
  geom_ribbon(aes(ymin = r0.l, ymax = r0.u, fill = city), alpha = 0.2, colour = NA) +
  geom_line() +
  # draw the R0 estimates from growth rate
  geom_point(data = df_sar_pt,
             aes(x = SAR, y = r0, col = city),
             size = 1) +
  geom_segment(data = df_sar_pt,
               aes(x = SAR, xend = SAR, y = 0, yend = r0, col = city),
               linetype = "dashed",
               linewidth = 0.6) +
  geom_segment(data = df_sar_pt,
               aes(x = 0, xend = SAR, y = r0, yend = r0, col = city),
               linetype = "dashed",
               linewidth = 0.6) +
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

p_r0_rel <- ggplot(filter(df_r0_proj, ngrp == "grp100"),
                   aes(x = city, y = mean_r0, col = city)) +
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

# sensitivity analysis comparing sexual activity group categorizations
df_r0_ngm <- df_r0_ngm %>% mutate(ngrp = ifelse(ngrp == "grp100", 
                                                "100 equally partitioned sexual activity groups",
                                                "248 sexual activity groups with finer partitions at the tail"))
df_sar_pt <- df_sar_pt %>% mutate(ngrp = ifelse(ngrp == "grp100", 
                                                "100 equally partitioned sexual activity groups",
                                                "248 sexual activity groups with finer partitions at the tail"))

p_r0_sen <- ggplot(df_r0_ngm,  
                   aes(x = SAR, y = mean_r0, col = ngrp)) +
  # line to show where R0=1 is
  geom_hline(yintercept = 1, linetype = "dotted") +
  
  # draw the R0 NGM estimate
  geom_ribbon(aes(ymin = r0.l, ymax = r0.u, fill = ngrp), alpha = 0.2, colour = NA) +
  geom_line() +
  # draw the R0 estimates from growth rate
  geom_point(data = df_sar_pt,
             aes(x = SAR, y = r0, col = ngrp),
             size = 1) +
  geom_segment(data = df_sar_pt,
               aes(x = SAR, xend = SAR, y = 0, yend = r0, col = ngrp),
               linetype = "dashed",
               linewidth = 0.6) +
  geom_segment(data = df_sar_pt,
               aes(x = 0, xend = SAR, y = r0, yend = r0, col = ngrp),
               linetype = "dashed",
               linewidth = 0.6) +
  coord_cartesian(ylim = c(0, 4)) +
  
  facet_wrap( ~ city) +
  
  # aesthetics
  scale_colour_viridis_d(option = "C", end = 0.8) +
  scale_fill_viridis_d(option = "C", end = 0.8) +
  
  # axis labels
  labs(title = "A)",
       x = "Secondary Attack Rate", 
       y = expression(italic("R"[0])*" (estimated)"),
       # y = bquote(atop(R["0"])),
       col = "Group Categorization", fill = "Group Categorization") +
  theme_cdf
