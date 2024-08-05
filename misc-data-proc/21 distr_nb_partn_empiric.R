# Library and data----
library(tidyverse)

## load data
data_3cities_pre_ipcw <- read_csv("../mpx-engage-params/data-processed/pre_ipcw_3cities.csv")
data_3cities_pand_ipcw <- read_csv("../mpx-engage-params/data-processed/pand_ipcw_3cities.csv")
data_3cities_post_ipcw <- read_csv("../mpx-engage-params/data-processed/post_ipcw_3cities.csv")

# create single dataset with all time periods & cities
data_3cities <- bind_rows(
  data_3cities_pre_ipcw,
  data_3cities_pand_ipcw,
  data_3cities_post_ipcw
) %>%
  mutate(time_pt = factor(time_pt, 
                          levels = c("Pre-Pandemic", "Pandemic", "Post-Restrictions")))

data_van <- filter(data_3cities, city == "van")

# recode age groups
unique(data_van$age_grp_cat)
data_van <- data_van %>% 
  mutate(
    age_grp_cat = case_when(age <  30 ~ "<30",
                            age <  50 ~ "30-49",
                            age >= 50 ~ "50+"),
    age_grp_cat = factor(age_grp_cat,
                        levels = c("<30", "30-49", "50+"))
  )

table(og = data_van$age, recode = data_van$age_grp_cat)

# Compute empiric distribution ----
# compute number of participants in each age x HIVx group
dat_nb_age_hiv <- data_van %>% 
  filter(time_pt == "Pre-Pandemic") %>% 
  group_by(city, age_grp_cat, hiv_stat) %>% 
  
  summarize(
    n = n(), n.rds = sum(wt_rds_norm),
    .groups = "drop"
  ) %>% 
  mutate(
    n_city = sum(n), n_city.rds = sum(n.rds),
    .after = city
  ) %>% 
  mutate(p_age = n / sum(n), p_age.rds = n.rds / sum(n.rds))
  
# verify that n_city.rds matches unweighted, and drop
dat_nb_age_hiv

dat_nb_age_hiv <- select(dat_nb_age_hiv, -n_city.rds)

# verify numbers and proportions
colSums(dat_nb_age_hiv[, c("n", "n.rds", "p_age", "p_age.rds")])

## compute the empiric distribution
table(data_van$nb_part_ttl)

# overall
dat_partn_distr_cty <- data_van %>% 
  group_by(city, time_pt, x = nb_part_ttl) %>% 
  
  summarize(
    n = n(), n.rds = sum(wt_rds_norm),
    .groups = "drop_last"
  ) %>% 
  mutate(p = n / sum(n), p.rds = n.rds / sum(n.rds)) %>% 
  ungroup()

# by age group
dat_partn_distr_by_age_hiv <- data_van %>% 
  group_by(city, time_pt, age_grp_cat, hiv = hiv_stat, x = nb_part_ttl) %>% 
  
  summarize(
    n = n(), n.rds = sum(wt_rds_norm),
    .groups = "drop_last"
  ) %>% 
  mutate(p = n / sum(n), p.rds = n.rds / sum(n.rds)) %>% 
  ungroup()

## pad dataframes
dat_partn_distr_cty <- dat_partn_distr_cty %>% 
  complete(
    city, time_pt, x = 0:300,
    fill = list(n = 0, n.rds = 0, p = 0, p.rds = 0)
  )

dat_partn_distr_by_age_hiv <- dat_partn_distr_by_age_hiv %>% 
  complete(
    city, time_pt, age_grp_cat, hiv, x = 0:300,
    fill = list(n = 0, n.rds = 0, p = 0, p.rds = 0)
  )

# verify padding (length(x) * [nb age groups * 2 for HIV groups] * nb time periods)
nrow(dat_partn_distr_cty)
length(0:300) * 3

nrow(dat_partn_distr_by_age_hiv)
length(0:300) * length(unique(data_van$age_grp_cat)) * 2 * 3

## verify the proportions sum to 1 and the n's sum to the city
# proportion in each age group
unique(dat_nb_age_hiv[, c("city", "n_city")])
colSums(dat_nb_age_hiv[, c("n", "n.rds", "p_age", "p_age.rds")])

# city-wide
tapply(dat_partn_distr_cty[, c("n", "n.rds", "p", "p.rds")], dat_partn_distr_cty$time_pt, colSums)

# by age & HIV group
tapply(dat_partn_distr_by_age_hiv[, c("n", "n.rds")], dat_partn_distr_by_age_hiv$time_pt, colSums)
tapply(dat_partn_distr_by_age_hiv$p,
       paste0(dat_partn_distr_by_age_hiv$time_pt, "_",
              dat_partn_distr_by_age_hiv$age_grp_cat,
              "_hiv", dat_partn_distr_by_age_hiv$hiv),
       sum)
tapply(dat_partn_distr_by_age_hiv$p.rds,
       paste0(dat_partn_distr_by_age_hiv$time_pt, "_",
              dat_partn_distr_by_age_hiv$age_grp_cat,
              "_hiv", dat_partn_distr_by_age_hiv$hiv),
       sum)

## save data ----
write.csv(x = dat_nb_age_hiv, file = "./misc-data-proc/outputs-cm/engage_van_age_hiv_distribution.csv", row.names = F)
write.csv(x = dat_partn_distr_cty, file = "./misc-data-proc/outputs-cm/engage_van_empiric_nb_partn_p6m.csv", row.names = F)
write.csv(x = dat_partn_distr_by_age_hiv, file = "./misc-data-proc/outputs-cm/engage_van_empiric_nb_partn_p6m_age_hiv.csv", row.names = F)

