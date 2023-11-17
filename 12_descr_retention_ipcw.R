# Libraries, functions and data ----
# libraries and functions
library(tidyverse)
library(lubridate)
library(survey)
theme_set(theme_bw())

source("./src/utils_helper.R")
source("./src/utils_ipcw_engage.R")

select <- dplyr::select
### key dates
## WHO declared a pandemic on March 11
# DATE_PAND_START <- as.Date("2020-03-11") %m+% months(3)

## earliest Engage visits after pandemic declared was June 17th 2020
DATE_PAND_ENGAGE <- as.Date("2020-06-01")

## new measures for fully vaccinated international travellers to Canada came into force 
# https://www.canada.ca/en/border-services-agency/news/2021/09/travel-advisory-reminder--on-september-7-new-measures-for-fully-vaccinated-international-travellers-to-canada-will-come-into-force.html
# (shifted by 3 months to account for P6M recall period)
DATE_RES_END <- as.Date("2021-09-07") %m+% months(3)

# load data
data_fu <- read_csv("../mpx-engage-params/data-3cities-feb-2023/engage_visits_3cities.csv")

# check visits between March-June 2020
# sort(unique(data_fu$date_intv[data_fu$date_intv >= "2020-03-01" & data_fu$date_intv <= "2020-06-30"]))

# create city marker
CITIES <- c("mtl", "trt", "van")
table(data_fu$city[data_fu$visit_num == 1], useNA = "ifany")

## date ranges for visits
tapply(data_fu$date_intv, data_fu$visit_num, summary)

# IPCW analysis (main analysis) ----
## Create datasets for every time period ----
## pre-pandemic time period
data_fu_pre <- data_fu %>% 
  group_by(part_id) %>%
  subset(visit_num == 1) %>% 
  mutate(time_pt = "Pre-Pandemic")
sum(data_fu_pre$wt_rds_norm)
length(unique(data_fu_pre$part_id))

## pandemic time period
data_fu_pand <- data_fu %>% 
  # keep only visits after pandemic start and before restrictions end
  filter(date_intv >= DATE_PAND_ENGAGE & date_intv < DATE_RES_END) %>%
  # keep earliest visit available
  group_by(part_id) %>% 
  filter(date_intv == min(date_intv)) %>% 
  mutate(time_pt = "Pandemic")

# post-restrictions time period
data_fu_post <- data_fu %>% 
  # keep only visits after restrictions end
  filter(date_intv >= DATE_RES_END) %>%
  # keep oldest visit available
  group_by(part_id) %>% 
  filter(date_intv == max(date_intv)) %>% 
  mutate(time_pt = "Post-Restrictions")

data_visit_key_date <- bind_rows(data_fu_pre, data_fu_pand, data_fu_post)

## Retention rate ----
# compute proportion retained at every time point
tbl_retain <- data_visit_key_date %>% 
  mutate(time_pt = factor(time_pt, levels = c("Pre-Pandemic", "Pandemic", "Post-Restrictions"))) %>% 
  group_by(city, time_pt) %>% 
  summarize(nb = n(), .groups = "drop_last") %>% 
  mutate(prop = round(nb / max(nb) * 100)) %>% 
  ungroup()

# format proportion (apostrophe used to prevent excel from reading as negatives)
tbl_retain <- tbl_retain %>% 
  mutate(prop = sprintf("'(%s%%)", prop))

# format to one column per city
tbl_retain <- tbl_retain %>% 
  pivot_wider(names_from = "city", values_from = c("nb", "prop")) %>% 
  select(time_pt, ends_with("_mtl"), ends_with("_trt"), ends_with("_van"))

write.csv(tbl_retain, "./out/manuscript-tables/table_S3_retention.csv", row.names = FALSE)

## Compute IPCWs ----
# lists to store each city's dataset with IPCWs
ipcw_pand <- create_city_list()
ipcw_post <- create_city_list()

covariates <- c("apps_partn_m",
                "apps_partn_d",
                "rel_status",
                "age", 
                "hiv_stat",
                "education_level_cat",
                "income_level_cat",
                "ethnicity_cat",
                "sex_work_m",
                "sex_work_d",
                "nb_part_ttl") # variables associated with ltfu

## check for LTFU and create indicator variable
# check if participant has a pandemic period visit
data_ipcw_pand <- data_fu_pre %>% 
  mutate(
    ltfu = (part_id %in% unique(data_fu_pand$part_id))
  )

# check if participant has a post-restrictions period visit
data_ipcw_post <- data_fu_pre %>% 
  mutate(
    ltfu = (part_id %in% unique(data_fu_post$part_id))
  )

# for both datasets, transform variables
data_ipcw_pand <- data_ipcw_pand %>% 
  mutate(across(all_of(covariates), as.factor)) %>%
  mutate(across(all_of(covariates), as.integer))
data_ipcw_post <- data_ipcw_post %>% 
  mutate(across(all_of(covariates), as.factor)) %>%
  mutate(across(all_of(covariates), as.integer))

## generate IPCW weights for all cities
for (cur_city in CITIES){
  ipcw_pand <- compute_ipcw(cur_city, covariates,
                            data_timept = data_ipcw_pand,
                            data_list = ipcw_pand)
}

for (cur_city in CITIES){
  ipcw_post <- compute_ipcw(cur_city, covariates,
                            data_timept = data_ipcw_post,
                            data_list = ipcw_pand)
}

# turn into a single dataset
ipcw_pand <- do.call(bind_rows, ipcw_pand)
ipcw_post <- do.call(bind_rows, ipcw_post)

ipcw_pre <- data_visit_key_date %>%
  subset(time_pt == "Pre-Pandemic") %>%
  mutate(ipw_rds = wt_rds_norm)

ipcw_pand <- merge.data.frame(data_visit_key_date, ipcw_pand, by = "part_id") %>%
  subset(time_pt == "Pandemic")

ipcw_post <- merge.data.frame(data_visit_key_date, ipcw_post, by = "part_id") %>%
  subset(time_pt == "Post-Restrictions")

# save datasets
write.csv(ipcw_pre,"../mpx-engage-params/data-3cities-feb-2023/pre_ipcw_3cities.csv", row.names = F)
write.csv(ipcw_pand,"../mpx-engage-params/data-3cities-feb-2023/pand_ipcw_3cities.csv", row.names = F)
write.csv(ipcw_post,"../mpx-engage-params/data-3cities-feb-2023/post_ipcw_3cities.csv", row.names = F)

## effective sample size at baseline ----
df_ess <- vector("list", 3)
names(df_ess) <- CITIES

for (cty in CITIES){
  print(sprintf("%s ===================", cty))
  data_tmp <- ipcw_pre %>% filter(city == cty & time_pt == "Pre-Pandemic")
  
  print( nrow(data_tmp) )
  print( sum(data_tmp$wt_rds) )
  print( sum(data_tmp$wt_rds_norm) )
  
  # get deff
  df_svy <- svydesign(ids = ~0, data = data_tmp, weights = ~wt_rds)
  avg <- svymean(~ nb_part_ttl, design = df_svy, deff = TRUE)
  
  df_ess[[cty]] <- data.frame(city = cty,
                              n = nrow(data_tmp),
                              deff = data.frame(avg)$deff,
                              ess = nrow(data_tmp) / data.frame(avg)$deff)
  
  # print(df_ess[[cty]]$ess)
  print( round(df_ess[[cty]]$ess, 2) )
}
df_ess <- bind_rows(df_ess)

write.csv(df_ess, "./out/manuscript-tables/table_1_ess.csv")

# Restriction analysis (sensitivity analysis) ----
## Get number of participants with complete data ----
# turn dataset into one row per participant, with columns for all visits
data_visits <- data_fu %>% 
  select(city, part_id, visit_num, date_intv) %>% 
  pivot_wider(names_from = visit_num, values_from = date_intv,
              names_prefix = "visit_")

# how many participants with at least {2,3} visits
data_visits %>% 
  group_by(city,
           has_2 = as.integer(!is.na(visit_2)), 
           has_3 = as.integer(!is.na(visit_3))) %>% 
  count() %>% 
  group_by(city) %>% 
  mutate(prop = n / sum(n)) %>% 
  arrange(city, -has_2, -has_3)

### create indicators for each subgroups
# (1) those with at least one data point in the pandemic (post March 2020)
# (2) those in (1) AND with at least one data point in 2022
names(data_visits)
data_visits <- data_visits %>% 
  mutate(
    pand_vis = case_when(
      visit_2 >= DATE_PAND_ENGAGE & visit_2 < DATE_RES_END | visit_3 >= DATE_PAND_ENGAGE & visit_3 < DATE_RES_END |
      visit_4 >= DATE_PAND_ENGAGE & visit_4 < DATE_RES_END | visit_5 >= DATE_PAND_ENGAGE & visit_5 < DATE_RES_END |
      visit_6 >= DATE_PAND_ENGAGE & visit_6 < DATE_RES_END | visit_7 >= DATE_PAND_ENGAGE & visit_7 < DATE_RES_END |
      visit_8 >= DATE_PAND_ENGAGE & visit_8 < DATE_RES_END | visit_9 >= DATE_PAND_ENGAGE & visit_9 < DATE_RES_END ~ 1,
      T ~ 0
    ),
    post_vis = case_when(
      visit_2 >= DATE_RES_END | visit_3 >= DATE_RES_END |
      visit_4 >= DATE_RES_END | visit_5 >= DATE_RES_END |
      visit_6 >= DATE_RES_END | visit_7 >= DATE_RES_END |
      visit_8 >= DATE_RES_END | visit_9 >= DATE_RES_END  ~ 1,
      T ~ 0
    )
  )

data_visits %>%
  filter(post_vis == 1 & pand_vis == 0) # note: 84 participants have post restriction visit but 
# don't have a pandemic visit

# see how many participants meet the criteria
data_visits %>% 
  group_by(city, pand_vis) %>% 
  summarize(n = n(), .groups = "drop_last") %>% 
  mutate(prop = n / sum(n)) %>% 
  filter(pand_vis == 1)

data_visits %>% 
  group_by(city, post_vis) %>% 
  summarize(n = n(), .groups = "drop_last") %>% 
  mutate(prop = n / sum(n)) %>% 
  filter(post_vis == 1)

data_fu <- full_join(
  select(data_visits, city, part_id, pand_vis, post_vis),
  data_fu,
  by = c("city", "part_id")
)

# participants with *both* pandemic and post restriction visit
## note: baseline visit, as indicated by the visit timeline, was before pandemic
## see code below "date ranges for visits"
data_res <- filter(data_fu, post_vis == 1 & pand_vis == 1) 
length(unique(data_res$part_id))

## create data.frame for baseline data 
data_res_pre <- data_res %>% 
  group_by(part_id) %>%
  subset(visit_num == 1) %>%
  mutate(time_pt = "Pre-Pandemic", .before = 1)

length(unique(data_res_pre$part_id))

## data of the closest point after the pandemic
data_res_pand <- data_res %>% 
  filter(date_intv >= DATE_PAND_ENGAGE & date_intv < DATE_RES_END) %>% 
  group_by(part_id) %>% 
  filter(date_intv == min(date_intv)) %>% 
  mutate(time_pt = "Pandemic", .before = 1)

# verify that the same visits as in IPCW are used
length(unique(data_res_pand$part_id))
ipcw_pand %>% 
  select(part_id, city, date_ipcw = date_intv) %>% 
  right_join(select(data_res_pand, part_id, date_restr = date_intv), by = "part_id") %>% 
  # check for discrepancies
  group_by(city) %>% 
  summarize(nb = n(), nb_dates_same = sum(date_ipcw == date_restr), nb_dates_diff = sum(date_ipcw != date_restr))

## data of the most recent point after the post restriction
# note: used most recent visit to capture post restriction behaviour
data_res_post <- data_res %>% 
  filter(date_intv >= DATE_RES_END)  %>%
  group_by(part_id) %>%
  filter(date_intv == max(date_intv)) %>%
  mutate(time_pt = "Post-Restrictions", .before = 1)

# verify that the same visits as in IPCW are used
length(unique(data_res_post$part_id))
ipcw_post %>% 
  select(part_id, city, date_ipcw = date_intv) %>% 
  right_join(select(data_res_post, part_id, date_restr = date_intv), by = "part_id") %>% 
  # check for discrepancies
  group_by(city) %>% 
  summarize(nb = n(), nb_dates_same = sum(date_ipcw == date_restr), nb_dates_diff = sum(date_ipcw != date_restr))

write.csv(data_res_pre,  "../mpx-engage-params/data-3cities-feb-2023/restriction/res_pre.csv", row.names = F)
write.csv(data_res_pand,  "../mpx-engage-params/data-3cities-feb-2023/restriction/res_pand.csv", row.names = F)
write.csv(data_res_post,  "../mpx-engage-params/data-3cities-feb-2023/restriction/res_post.csv", row.names = F)

# verify that visits don't overlap
range(data_res_pre$date_intv)
range(data_res_pand$date_intv)
range(data_res_post$date_intv)

# save dataset
data_res_fu <- rbind(data_res_pre,data_res_pand, data_res_post)
table(visit = data_res_fu$time_pt, data_res_fu$city)

# save dataset
write_csv(data_res_fu,  "../mpx-engage-params/data-3cities-feb-2023/restriction/res_fu.csv")
