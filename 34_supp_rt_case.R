# Libraries and data -----
library(tidyverse)
library(lubridate)
library(EpiEstim)
library(gridExtra)

theme_set(theme_bw())

## data for Rt estimation (only need the 3 concerned provinces)
PROVS <- c("Québec", "Ontario", "British Columbia")
data_mpox <- read_csv("./data-public/monkeypox-detailed-2023june13.csv")
data_mpox$reporting_pt_en[data_mpox$reporting_pt_en == "Quebec"] <- "Québec"

data_mpox <- subset(data_mpox, reporting_pt_en %in% PROVS)

# Rt estimation to assess outbreak trajectories ----
## Data processing ----
# rename variables
data_mpox <- data_mpox %>% 
  select(pruid, prov_en = reporting_pt_en, prov_fr = reporting_pt_fr,
         date, case_delta = num_confirmedcases_delta, case_cumul = num_confirmedcases_cumulative)

data_mpox$prov_en <- factor(data_mpox$prov_en, levels = PROVS)

# pad with 0's for days without cases
df_pad <- expand.grid(prov_en = PROVS, date = seq(min(data_mpox$date), max(data_mpox$date), by = "1 day"))
df_pad <- left_join(df_pad,
                    unique(data_mpox[, c("pruid", "prov_en", "prov_fr")]),
                    by = "prov_en")

data_mpox <- full_join(data_mpox, df_pad, by = c("pruid", "prov_en", "prov_fr", "date")) %>% 
  arrange(prov_en, date)

## keep only from day since first case up to 200 days after
data_mpox %>% filter(!is.na(case_delta)) %>% group_by(prov_en) %>% summarize(min(date), max(date))
data_mpox <- data_mpox %>% 
  filter(prov_en == "Québec"           & date >= "2022-04-28" & date <= as.Date("2022-04-28") + 200 |
         prov_en == "Ontario"          & date >= "2022-05-13" & date <= as.Date("2022-05-13") + 200 |
         prov_en == "British Columbia" & date >= "2022-05-25" & date <= as.Date("2022-05-25") + 200)

data_mpox <- data_mpox %>% 
  # complete incidence
  mutate(case_delta = ifelse(is.na(case_delta), 0, case_delta)) %>% 
  # complete cumulative
  group_by(pruid) %>% 
  fill(case_cumul, .direction = "down") %>% 
  ungroup()

# create dummy date
data_mpox <- data_mpox %>% 
  group_by(pruid) %>% 
  mutate(t_since_1 = as.integer(date - min(date)), .after = date) %>% 
  ungroup()

# Estimate Rt and plot ----
window <- 14 # how much to smooth the Rt estimate

data_rt <- vector("list", 3)
for(i in 1:length(PROVS)){
  incid <- filter(data_mpox, prov_en == PROVS[i])
  
  # define time windows over which to compute the Rt
  t_start <- seq(2, as.numeric(max(incid$date) - min(incid$date)) - window + 1)
  t_end <- t_start + window 
  
  # note: default when t_start and t_end are not used is to use weekly sliding windows
  #       this should be OK for now, could maybe consider a shorter window
  r_estim <- estimate_R(
    incid = incid$case_delta,
    method = "parametric_si",
    config = make_config(
      t_start = t_start,
      t_end = t_end,
      mean_si = 8.8,
      std_si = 8.7,
      mean_prior = 3,
      std_prior = 2
    )
  )
  
  data_rt[[i]] <- mutate(r_estim$R, prov_en = PROVS[i], .before = 1)
  data_rt[[i]] <- rename(data_rt[[i]], rt_mean = `Mean(R)`,
                         lci = `Quantile.0.025(R)`, uci = `Quantile.0.975(R)`)
}
data_rt <- bind_rows(data_rt) %>% as_tibble()

# add dates (plot each estimate at the end of the window)
df_time <- data_mpox %>% 
  mutate(t_end = t_since_1 + 1) %>% 
  select(prov_en, date, t_end)

data_rt <- full_join(df_time, data_rt, by = c("prov_en", "t_end"))

# provinces in order
data_rt$prov_en <- factor(data_rt$prov_en, levels = PROVS)

# see when Rt stablizes
png("./fig/fig_S7_rt.png",
    width = 16, height = 12, units = "cm", res = 600)

ggplot(data_rt, aes(x = t_end)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_ribbon(aes(ymin = lci, ymax = uci, fill = prov_en), alpha = 0.2) +
  geom_line(aes(y = rt_mean, col = prov_en)) +
  facet_wrap(~prov_en, ncol = 1) +
  coord_cartesian(ylim = c(0, 3)) +
  labs(x = "Days since first case", y = expression(italic('R'["e"])), col = "Province", fill = "Province") +
  theme(legend.position = "top") +
  scale_colour_viridis_d(option = "C", end = 0.8) +
  scale_fill_viridis_d(option = "C", end = 0.8) +
  theme(
    legend.position = "none"
  )

dev.off()

