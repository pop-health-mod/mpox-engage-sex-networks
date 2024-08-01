# Library and data----
library(tidyverse)

## load data
data_3cities <- read_csv("../mpx-engage-params/data-processed/engage_baseline_3cities.csv")

# create city marker
CITIES <- c("Montreal", "Toronto", "Vancouver")
# TIMEPTS <- c("Pre-Pandemic", "Pandemic", "Post-Restrictions")

# recode names of partnership status
unique(data_3cities$rel_status)
data_3cities <- data_3cities %>% 
  mutate(
    rel_status_orig = rel_status,
    rel_status = factor(rel_status,
                        levels = c("open", "exclusive", "unclear", "no relationship"),
                        labels = c("main-open", "main-exclusive", "main-unclear", "no-main"))
  )

table(og = data_3cities$rel_status_orig, recode = data_3cities$rel_status)

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

## pad dataframe
dat_partn_distr <- dat_partn_distr %>% 
  complete(
    city, rel_status, x = 0:300,
    fill = list(n = 0, n.rds = 0, p = 0, p.rds = 0)
  )


# verify padding
nrow(dat_partn_distr)
length(0:300) * length(unique(data_3cities$city)) * length(unique(data_3cities$rel_status))

## verify the proportions sum to 1 and the n's sum to the city
unique(dat_nb_rel[, c("city", "n_city")])

tapply(dat_partn_distr$n, dat_partn_distr$city, sum)
tapply(dat_partn_distr$n.rds, dat_partn_distr$city, sum)

tapply(dat_partn_distr$p, paste(dat_partn_distr$city, dat_partn_distr$rel_status, sep = "-"), sum)
tapply(dat_partn_distr$p.rds, paste(dat_partn_distr$city, dat_partn_distr$rel_status, sep = "-"), sum)

## save data ----
write.csv(x = dat_nb_rel, file = "./misc-data-proc/outputs-jk/engage_relat_distribution.csv", row.names = F)
write.csv(x = dat_partn_distr, file = "./misc-data-proc/outputs-jk/engage_nb_partn_p6m_by_relat_stat.csv", row.names = F)

