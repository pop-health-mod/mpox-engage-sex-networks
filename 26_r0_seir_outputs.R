# Library and data----
library(tidyverse)
library(gridExtra)
library(ggtext)

theme_set(theme_bw())

# incidence data
data_incid <- read_csv("./data-public/mpox_case_data_phac.csv")
data_incid$city <- case_when(data_incid$prov == "QC" ~ "mtl",
                             data_incid$prov == "ON" ~ "tor",
                             data_incid$prov == "BC" ~ "van")
data_incid$city_name <- case_when(data_incid$prov == "QC" ~ "Montréal",
                                  data_incid$prov == "ON" ~ "Toronto",
                                  data_incid$prov == "BC" ~ "Vancouver")

# indicate up to when data is fit
data_incid <- data_incid %>% 
  mutate(
    fit_target = (prov == "QC" & date <= "2022-06-15" |
                    prov == "ON" & date <= (as.Date("2022-07-10")) |
                    prov == "BC" & date <= (as.Date("2022-07-10")))
  )

## model fit
data_mod_fit <- read_csv("./out-seir/fit_incidence.csv")
data_mod_fit$date <- as.Date(data_mod_fit$date)

# finalize city label
data_mod_fit$city_name <- case_when(data_mod_fit$city == "mtl" ~ "Montréal",
                                    data_mod_fit$city == "tor" ~ "Toronto",
                                    data_mod_fit$city == "van" ~ "Vancouver")

# make figure ----
## panel A: model fit ----

# TODO this should go inside a function so that we can use it for each city
tapply(data_incid$date, data_incid$prov, range)

# manuscript colour palette
# col_pal <- viridis::viridis(n = 3, option = "C", end = 0.8)

p_mod_a <- ggplot(subset(data_mod_fit, time %% 1 == 0), aes(x = date, y = cases, col = city)) +
  # model fit and CrI's
  # geom_ribbon(aes(ymin = cases_cri.l, ymax = cases_cri.u, fill = city, col = NULL), alpha = 0.4) +
  geom_line(aes(linetype = "Model fit"), linewidth = 1) +
  
  # data to which model is fit [should make it a bit bigger]
  geom_point(
    data = subset(data_incid, fit_target),
    aes(y = incidence, shape = "Observed"),
    size = 1.8, alpha = .8
  ) +
  # geom_step(
  #   data = subset(data_incid, fit_target),
  #   aes(y = incidence)
  # ) +
  # full timseries of incidence
  # geom_point(
  #   data = subset(data_incid, !fit_target),
  #   aes(y = incidence),
  #   alpha = .25,
  #   fill = NA
  # ) +

# plot scale and 
  facet_grid(~ city_name) +
  coord_cartesian(xlim = as.Date(c("2022-04-20", "2022-09-01"))) +
  
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  coord_cartesian(xlim = c(min(data_incid$date), as.Date("2022-07-15"))) +
  # scale_y_continuous(breaks = 0:5 * 2) +
  
  # colour scheme
  scale_colour_viridis_d(option = "C", end = 0.8) +
  scale_fill_viridis_d(option = "C", end = 0.8) +
  
  # legends
  labs(title = "A)",
       x = "Date", y = "Reported mpox cases (day)", linetype = NULL, shape = NULL) +
  guides(
    shape = guide_legend(order = 1, override.aes = list(size = 2.8)),
    colour = "none",
    fill = "none"
  ) +
  theme(
    legend.spacing.y = unit(0, 'cm')
  ); p_mod_a

## panel B: NGM R_eff depleting groups ----
data_r0 <- read_csv("./out-seir/fit_r0.csv")

# finalize city label
data_r0$city_name <- case_when(data_r0$city == "mtl" ~ "Montréal",
                               data_r0$city == "tor" ~ "Toronto",
                               data_r0$city == "van" ~ "Vancouver")

min(data_r0$prop_imm[data_r0$prop_imm > 0] * 100) # the lowest % depleted

# x axis breaks
# reff_x_breaks <- 0:4 * .25
# reff_x_labs <- c("0.0%", "0.25%", "0.5%", "0.75%", "1.0%")

p_mod_b <- ggplot(data_r0, aes(x = prop_imm * 100, y = pt_est)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  
  # R0 point estimate and CrI's
  geom_ribbon(aes(ymin = lci, ymax = uci, fill = city), alpha = 0.4) +
  geom_line(aes(col = city), linewidth = 1) +
  facet_grid(~ city_name) +
  
  # use points instead
  # geom_pointrange(aes(ymin = r0_lci, ymax = r0_uci, col = city), size = 0.3, linewidth = 0.8,
  #                 position = position_dodge(width = 0.045)) +
  # geom_line(aes(col = city), position = position_dodge(width = 0.045)) +
  
  # plot scale and 
  coord_cartesian(xlim = c(0, 3), ylim = c(0, 4)) +
  
  # colour scheme
  scale_colour_viridis_d(option = "C", end = 0.8) +
  scale_fill_viridis_d(option = "C", end = 0.8) +
  
  # legends
  # scale_x_continuous(breaks = reff_x_breaks, labels = reff_x_labs) +
  labs(title = "B)",
       x = "Proportion of the highest-activity population immune (%)",
       y = expression(italic("R"["eff"])*" (estimated)"),
       col = "City") +
  guides(
    shape = "none", colour = "none", fill = "none"
  ) +
  theme(
    legend.position = c(.9, .82)
  ); p_mod_b

# fit_pars_tor <- readRDS("./model_pars_TOR_all_fits.rds")

## panel C: projected R0 ----
## now project the R0
data_r0_pre <- read_csv("./out-seir/fit_r0_prepand.csv")

# finalize city label
data_r0_pre$city_name <- case_when(data_r0_pre$city == "mtl" ~ "Montréal",
                                   data_r0_pre$city == "tor" ~ "Toronto",
                                   data_r0_pre$city == "van" ~ "Vancouver")

data_r0_prepost <- bind_rows(
  data_r0 %>% mutate(type = "Estimated from NGM"),
  data_r0_pre %>% 
    mutate(type = "Projected based\non pre-pandemic\nsex. activity")
) #%>% filter(prop_imm == 0)

data_r0_prepost

data_r0_prepost$type <- factor(data_r0_prepost$type, levels = rev(unique(data_r0_prepost$type)))

# finalize city label
data_r0_prepost$city_name <- case_when(data_r0_prepost$city == "mtl" ~ "Montréal",
                                       data_r0_prepost$city == "tor" ~ "Toronto",
                                       data_r0_prepost$city == "van" ~ "Vancouver")

p_mod_c <- ggplot(subset(data_r0_prepost, type != "Estimated from NGM"), aes(x = prop_imm * 100, y = pt_est)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  
  # R0 point estimate and CrI's
  geom_ribbon(aes(ymin = lci, ymax = uci, fill = city), alpha = 0.4) +
  geom_line(data = data_r0_prepost,
            aes(x = prop_imm * 100, col = city, linetype = type), linewidth = 1) +
  
  # geom_line(data = subset(data_r0_prepost, type == "Estimated from NGM"), aes(col = city),
  #           linetype = "dashed", linewidth = 0.5) +
  
  facet_grid(~ city_name) +
  
  # use points instead
  # geom_pointrange(aes(ymin = r0_lci, ymax = r0_uci, col = city), size = 0.3, linewidth = 0.8,
  #                 position = position_dodge(width = 0.045)) +
  # geom_line(aes(col = city), position = position_dodge(width = 0.045)) +
  
  # plot scale and 
  coord_cartesian(xlim = c(0, 3), ylim = c(0, 4)) +
  
  # colour scheme
  scale_colour_viridis_d(option = "C", end = 0.8) +
  scale_fill_viridis_d(option = "C", end = 0.8) +
  
  # legends
  # scale_x_continuous(breaks = reff_x_breaks, labels = reff_x_labs) +
  labs(title = "C)",
       x = "Proportion of the highest-activity population immune (%)",
       y = expression(italic("R"["eff"])*" (projected)"),
       linetype = NULL) +
  guides(
    shape = "none",
    colour = "none",
    fill = "none"
  ) +
  theme(
    legend.position = c(.2, .78)
  ); p_mod_c

grid.arrange(p_mod_a, p_mod_b, p_mod_c, ncol = 1)

png("./fig/fig_2_seir_results_ngm_r0.png",
    width = 24, height = 24, units = "cm", res = 600)

grid.arrange(p_mod_a, p_mod_b, p_mod_c, ncol = 1)

dev.off()
