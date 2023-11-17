# Library and data----
library(tidyverse)
library(gridExtra)
library(ggtext)

theme_set(theme_bw())

source("./src/plot_config.R")

# dates of vaccination start
vax_start <- tibble(prov = c("QC", "ON", "BC"),
                    city_name = c("Montréal", "Toronto", "Vancouver"),
                    date_scaleup = as.Date(c("2022-06-15", "2022-07-10", "2022-07-10")))

# incidence data
data_incid <- read_csv("./data-public/mpox_case_data_phac.csv")
data_incid$city <- case_when(data_incid$prov == "QC" ~ "mtl",
                             data_incid$prov == "ON" ~ "tor",
                             data_incid$prov == "BC" ~ "van")
data_incid$city_name <- case_when(data_incid$prov == "QC" ~ "Montréal",
                                  data_incid$prov == "ON" ~ "Toronto",
                                  data_incid$prov == "BC" ~ "Vancouver")

## model fit
data_mod_fit <- read_csv("./out-seir/fit_incidence.csv")
data_mod_fit$date <- as.Date(data_mod_fit$date)

# finalize city label
data_mod_fit$city_name <- case_when(data_mod_fit$city == "mtl" ~ "Montréal",
                                    data_mod_fit$city == "tor" ~ "Toronto",
                                    data_mod_fit$city == "van" ~ "Vancouver")

## indicate up to when data is fit in observed and modelled
data_incid <- data_incid %>% 
  mutate(
    fit_target = (prov == "QC" & date <= vax_start$date_scaleup[1] |
                    prov == "ON" & date <= vax_start$date_scaleup[2] |
                    prov == "BC" & date <= vax_start$date_scaleup[3])
  )

data_mod_fit <- data_mod_fit %>% 
  mutate(
    fit_target = (city == "mtl" & date <= vax_start$date_scaleup[1] |
                    city == "tor" & date <= vax_start$date_scaleup[2] |
                    city == "van" & date <= vax_start$date_scaleup[3])
  )

# Make figures ----
## Figure 2: model fit ----

tapply(data_incid$date, data_incid$prov, range)

# manuscript colour palette for cities/provinces
# col_pal <- viridis::viridis(n = 3, option = "C", end = 0.8)

# plot fit (for modelled, need to subset to daily data to plot properly)
data_incid <- subset(data_incid, fit_target)
data_mod_fit <- subset(data_mod_fit, fit_target & time %% 1 == 0)

p_mod_fit <- ggplot(data_mod_fit, aes(x = date, y = cases, col = city)) +
  # date of vaccination start
  geom_vline(
    data = vax_start, aes(xintercept = date_scaleup), linetype = "dashed"
  ) +
  
  # model fit and CrI's
  geom_ribbon(aes(ymin = cases_lci, ymax = cases_uci, fill = city, col = NULL), alpha = 0.4) +
  geom_line(aes(linetype = "Model fit"), linewidth = 0.7) +
  
  # observed incidence data
  geom_point(
    data = data_incid,
    aes(y = incidence, shape = "Observed"),
    size = 1.4, alpha = .4
  ) +

  # plot scale and facetting
  facet_wrap(~city_name, ncol = 1) +
  coord_cartesian(xlim = range(data_incid$date)) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  
  # colour scheme
  scale_colour_viridis_d(option = "C", end = 0.8) +
  scale_fill_viridis_d(option = "C", end = 0.8) +
  
  # legends
  labs(x = "Date", y = "Reported mpox cases (day)", linetype = NULL, shape = NULL) +
  guides(
    shape = guide_legend(order = 1, override.aes = list(size = 2.8)),
    colour = "none",
    fill = "none"
  ) +
  theme_seir +
  theme(
    # legend.position = "bottom",
    legend.box.background = element_rect(fill = "transparent", colour = "transparent"),
    legend.background = element_rect(colour = "transparent", fill = alpha("white", 0.4)),
    legend.position = c(.82, .91),
    legend.spacing.y = unit(-2, "pt"),
    legend.direction = "vertical"
  ); p_mod_fit

png("./fig/fig_2_seir_model_fit.png",
    width = 8, height = 12, units = "cm", res = 600)
p_mod_fit
dev.off()
pdf("./fig/fig_2_seir_model_fit.pdf",
    width = 8 / 2.54, height = 12 / 2.54)
p_mod_fit
dev.off()

### Panel A: NGM R_e depleting groups ----
data_r0 <- read_csv("./out-seir/fit_r0.csv")

# finalize city label
data_r0$city_name <- case_when(data_r0$city == "mtl" ~ "Montréal",
                               data_r0$city == "tor" ~ "Toronto",
                               data_r0$city == "van" ~ "Vancouver")

min(data_r0$prop_imm[data_r0$prop_imm > 0] * 100) # the lowest % depleted

# x axis breaks
# re_x_breaks <- 0:4 * .25
# re_x_labs <- c("0.0%", "0.25%", "0.5%", "0.75%", "1.0%")

p_r0_a <- ggplot(data_r0, aes(x = prop_imm * 100, y = pt_est)) +
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
  # scale_x_continuous(breaks = re_x_breaks, labels = re_x_labs) +
  theme_seir +
  labs(title = "A)",
       x = "Proportion of the highest-activity population immune (%)",
       y = expression(italic("R"["e"])*" (estimated)"),
       col = "City") +
  guides(
    shape = "none", colour = "none", fill = "none"
  ); #p_r0_a

### Panel B: projected R0 ----
data_r0_pre <- read_csv("./out-seir/fit_r0_prepand.csv")

# finalize city label
data_r0_pre$city_name <- case_when(data_r0_pre$city == "mtl" ~ "Montréal",
                                   data_r0_pre$city == "tor" ~ "Toronto",
                                   data_r0_pre$city == "van" ~ "Vancouver")

data_r0_prepost <- bind_rows(
  data_r0 %>% mutate(type = "Estimated"),
  data_r0_pre %>% 
    mutate(type = "Projected based\non pre-pandemic\nsexual activity")
) #%>% filter(prop_imm == 0)

data_r0_prepost

data_r0_prepost$type <- factor(data_r0_prepost$type, levels = rev(unique(data_r0_prepost$type)))

# finalize city label
data_r0_prepost$city_name <- case_when(data_r0_prepost$city == "mtl" ~ "Montréal",
                                       data_r0_prepost$city == "tor" ~ "Toronto",
                                       data_r0_prepost$city == "van" ~ "Vancouver")

p_r0_b <- ggplot(subset(data_r0_prepost, type != "Estimated"), aes(x = prop_imm * 100, y = pt_est)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  
  # R0 point estimate and CrI's
  geom_ribbon(aes(ymin = lci, ymax = uci, fill = city), alpha = 0.4) +
  geom_line(data = data_r0_prepost,
            aes(x = prop_imm * 100, col = city, linetype = type), linewidth = 0.7) +
  
  facet_grid(~ city_name) +
  
  # plot scale and 
  coord_cartesian(xlim = c(0, 3), ylim = c(0, 4)) +
  
  # colour scheme
  scale_colour_viridis_d(option = "C", end = 0.8) +
  scale_fill_viridis_d(option = "C", end = 0.8) +
  
  # legends
  # scale_x_continuous(breaks = re_x_breaks, labels = re_x_labs) +
  labs(title = "B)",
       x = "Proportion of the highest-activity population immune (%)",
       y = expression(italic("R"["e"])*" (projected)"),
       linetype = NULL) +
  guides(
    shape = "none",
    colour = "none",
    fill = "none"
  ) +
  theme_seir +
  theme(
    legend.position = c(.235, .76),
    legend.background = element_rect(colour = "transparent", fill = alpha("white", 0.4)),
    legend.spacing.y = unit(-2, "pt"),
  ); #p_mod_c

### Output Figure 3 ----
grid.arrange(p_r0_a, p_r0_b, ncol = 1)

png("./fig/fig_3_seir_ngm_r0.png",
    width = 18, height = 12, units = "cm", res = 600)
grid.arrange(p_r0_a, p_r0_b, ncol = 1)
dev.off()
pdf("./fig/fig_3_seir_ngm_r0.pdf",
    width = 18 / 2.54, height = 12 / 2.54)
grid.arrange(p_r0_a, p_r0_b, ncol = 1)
dev.off()


