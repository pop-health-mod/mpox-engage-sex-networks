# Library and data----
library(tidyverse)

## main analyses
outcome_var <- "nb_part_ttl"

## define paths & prefixes based on the analysis being done
fig_path <- "./misc-data-proc-JK/outputs"
out_distr_path <- "./misc-data-proc-JK/outputs"

## load data
data_3cities_pre_ipcw <- read_csv("../mpx-engage-params/data-processed/pre_ipcw_3cities.csv")
# data_3cities_pand_ipcw <- read_csv("../mpx-engage-params/data-processed/pand_ipcw_3cities.csv")
# data_3cities_post_ipcw <- read_csv("../mpx-engage-params/data-processed/post_ipcw_3cities.csv")

# create single dataset with all time periods & cities
# data_3cities <- bind_rows(
#   data_3cities_pre_ipcw,
#   data_3cities_pand_ipcw,
#   data_3cities_post_ipcw
# ) %>%
#   mutate(time_pt = factor(time_pt, 
#                           levels = c("Pre-Pandemic", "Pandemic", "Post-Restrictions"))) %>%
#   mutate(city = recode_factor(city, mtl = "Montreal", trt = "Toronto", van = "Vancouver"))
data_3cities <- data_3cities_pre_ipcw

# create city marker
CITIES <- c("Montreal", "Toronto", "Vancouver")
# TIMEPTS <- c("Pre-Pandemic", "Pandemic", "Post-Restrictions")

# table(data_3cities$time_pt, 
#       data_3cities$city, 
#       useNA = "ifany")
# 
# data_3cities <- data_3cities %>% 
#   mutate(
#     data_pt = factor(paste(city, time_pt, sep = "-"),
#                      levels = paste(rep(CITIES, each = 3), TIMEPTS, sep = "-"))
#   )

# CITIES_DATAPTS <- paste(
#   rep(CITIES, each = 3), 
#   rep(c("Pre-Pandemic", "Pandemic", "Post-Restrictions"), 
#       times = 3), 
#   sep = "-"
# )


# Compute empiric distribution ----
# compute number of participants in each relationship status
dat_nb_rel <- data_3cities %>% 
  group_by(city, rel_status) %>% 
  
  summarize(
    n = n(), n.rds = sum(wt_rds_norm),
    .groups = "drop_last"
  ) %>% 
  mutate(
    n_city = sum(n), n_city.rds = sum(n.rds),
    .after = city
  ) %>% 
  mutate(p_rel = n / sum(n), p_rel.rds = n.rds / sum(n.rds)) %>% 
  ungroup()
  
# verify that n_city.rds matches unweighted, and drop
dat_nb_rel

dat_nb_rel <- select(dat_nb_rel, -n_city.rds)

# verify numbers and proportions
tapply(dat_nb_rel$n, dat_nb_rel$city, sum)
tapply(dat_nb_rel$n.rds, dat_nb_rel$city, sum)

tapply(dat_nb_rel$p_rel, dat_nb_rel$city, sum)
tapply(dat_nb_rel$p_rel.rds, dat_nb_rel$city, sum)

# compute the empiric distribution
table(cty = data_3cities$city, nb_part_ttl = data_3cities$nb_part_ttl)
dat_partn_distr <- data_3cities %>% 
  group_by(city, rel_status, x = nb_part_ttl) %>% 
  
  summarize(
    n = n(), n.rds = sum(wt_rds_norm),
    .groups = "drop_last"
  ) %>% 
  mutate(p = n / sum(n), p.rds = n.rds / sum(n.rds)) %>% 
  ungroup()


## verify the proportions sum to 1 and the n's sum to the city
unique(dat_nb_rel[, c("city", "n_city")])

tapply(dat_partn_distr$n, dat_partn_distr$city, sum)
tapply(dat_partn_distr$n.rds, dat_partn_distr$city, sum)

tapply(dat_partn_distr$p, paste(dat_partn_distr$city, dat_partn_distr$rel_status, sep = "-"), sum)
tapply(dat_partn_distr$p.rds, paste(dat_partn_distr$city, dat_partn_distr$rel_status, sep = "-"), sum)

## save data ----
write.csv(x = dat_nb_rel, file = "./misc-data-proc-JK/outputs/engage_relat_distribution.csv", row.names = F)
write.csv(x = dat_partn_distr, file = "./misc-data-proc-JK/outputs/engage_nb_partn_p6m_by_relat_stat.csv", row.names = F)

