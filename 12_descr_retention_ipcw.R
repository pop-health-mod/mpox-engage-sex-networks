# Libraries, functions and data ----
# libraries and functions
library(tidyverse)
library(lubridate)
library(survival)
library(survminer)
library(gridExtra)

theme_set(theme_bw())

source("./src/utils.R")
source("./src/utils_helper.R")
source("./src/plot.R")
source("./src/ipcw_engage.R")

# key dates
# urgence sanitaire declared on March 13
# DATE_PAND_START <- as.Date("2020-03-13")%m+% months(3)
# WHO declared a pandemic on March 11
DATE_PAND_START <- as.Date("2020-03-11") %m+% months(3) #to account for recall period "P6M"
DATE_RES_END <- as.Date("2021-09-07") %m+% months(3) #new measures for fully vaccinated international travellers to Canada came into force 
#https://www.canada.ca/en/border-services-agency/news/2021/09/travel-advisory-reminder--on-september-7-new-measures-for-fully-vaccinated-international-travellers-to-canada-will-come-into-force.html
# DATE_MPOX_START <- as.Date("2022-05-19") %m+% months(3)
# DATE_MPOX_STABILIZED <- as.Date("2022-09-01")%m+% months(3)

## which variable to use as outcome
outcome_var <- "nb_part_ttl" 
file_suff <- case_when(outcome_var == "nb_part_ttl" ~ "p6m_all",
                       outcome_var == "nb_part_new" ~ "p6m_new",
                       outcome_var == "nb_part_anal" ~ "p6m_anal")

data_fu <- read_csv("./data-3cities-feb-2023/engage_visits_3cities.csv")

# create city marker
CITIES <- c("mtl", "trt", "van")
table(data_fu$city[data_fu$visit_num == 1], useNA = "ifany")

## date ranges for visits
tapply(data_fu$date_intv, data_fu$visit_num, summary)

# Restriction analysis ----
## Get nb of participants with complete data ----
data_visits <- data_fu %>% 
  dplyr::select(city, part_id, visit_num, date_intv) %>% 
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
data_visits <- data_visits %>% 
  mutate(
    pand_vis = case_when(visit_2 >= DATE_PAND_START&visit_2 < DATE_RES_END | visit_3 >= DATE_PAND_START&visit_3 < DATE_RES_END |
                           visit_4 >= DATE_PAND_START&visit_4 < DATE_RES_END | visit_5 >= DATE_PAND_START&visit_5 < DATE_RES_END |
                           visit_6 >= DATE_PAND_START&visit_6 < DATE_RES_END | visit_7 >= DATE_PAND_START&visit_7 < DATE_RES_END |
                           visit_8 >= DATE_PAND_START&visit_8 < DATE_RES_END | visit_9 >= DATE_PAND_START&visit_9 < DATE_RES_END ~ 1,
                         T ~ 0),
    post_vis = case_when(visit_2 >= DATE_RES_END | visit_3 >= DATE_RES_END |
                          visit_4 >= DATE_RES_END | visit_5 >= DATE_RES_END |
                          visit_6 >= DATE_RES_END | visit_7 >= DATE_RES_END |
                          visit_8 >= DATE_RES_END | visit_9 >= DATE_RES_END  ~ 1,
                        T ~ 0))
    # pre_mpox_vis = case_when(visit_2 >= DATE_RES_END&visit_2 < DATE_MPOX_START | visit_3 >= DATE_RES_END&visit_3 < DATE_MPOX_START |
    #                            visit_4 >= DATE_RES_END&visit_4 < DATE_MPOX_START | visit_5 >= DATE_RES_END&visit_5 < DATE_MPOX_START |
    #                            visit_6 >= DATE_RES_END&visit_6 < DATE_MPOX_START | visit_7 >= DATE_RES_END&visit_7 < DATE_MPOX_START |
    #                            visit_8 >= DATE_RES_END&visit_8 < DATE_MPOX_START | visit_9 >= DATE_RES_END&visit_9 < DATE_MPOX_START ~ 1,
    #                          T ~ 0),
    # during_mpox_vis = case_when(visit_2 >= DATE_MPOX_START&visit_2 < DATE_MPOX_STABILIZED | visit_3 >= DATE_MPOX_START&visit_3 < DATE_MPOX_STABILIZED |
    #                        visit_4 >= DATE_MPOX_START&visit_4 < DATE_MPOX_STABILIZED | visit_5 >= DATE_MPOX_START&visit_5 < DATE_MPOX_STABILIZED |
    #                        visit_6 >= DATE_MPOX_START&visit_6 < DATE_MPOX_STABILIZED | visit_7 >= DATE_MPOX_START&visit_7 < DATE_MPOX_STABILIZED |
    #                        visit_8 >= DATE_MPOX_START&visit_8 < DATE_MPOX_STABILIZED | visit_9 >= DATE_MPOX_START&visit_9 < DATE_MPOX_STABILIZED ~ 1,
    #                      T ~ 0),
    # post_mpox_vis = case_when(visit_2 >= DATE_MPOX_STABILIZED | visit_3 >= DATE_MPOX_STABILIZED |
    #                        visit_4 >= DATE_MPOX_STABILIZED | visit_5 >= DATE_MPOX_STABILIZED |
    #                        visit_6 >= DATE_MPOX_STABILIZED | visit_7 >= DATE_MPOX_STABILIZED |
    #                        visit_8 >= DATE_MPOX_STABILIZED | visit_9 >= DATE_MPOX_STABILIZED  ~ 1,
    #                      T ~ 0)) 

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
# data_visits %>% 
#   group_by(city, pre_mpox_vis) %>% 
#   summarize(n = n(), .groups = "drop_last") %>% 
#   mutate(prop = n / sum(n)) %>% 
#   filter(pre_mpox_vis == 1)
# data_visits %>% 
#   group_by(city, during_mpox_vis) %>% 
#   summarize(n = n(), .groups = "drop_last") %>% 
#   mutate(prop = n / sum(n)) %>% 
#   filter(during_mpox_vis == 1)
# data_visits %>% 
#   group_by(city, post_mpox_vis) %>% 
#   summarize(n = n(), .groups = "drop_last") %>% 
#   mutate(prop = n / sum(n)) %>% 
#   filter(post_mpox_vis == 1)

data_fu <- full_join(
  dplyr::select(data_visits, city, part_id, pand_vis, post_vis),
  data_fu,
  by = c("city", "part_id")
)

# participants with *both* pandemic and post restriction visit
## note: baseline visit, as indicated by the visit timeline, was before pandemic
## see code below "date ranges for visits"
data_res <- filter(data_fu, post_vis == 1 & pand_vis == 1) 
length(unique(data_res$part_id))

## baseline data 
data_res_pre <- data_res %>% 
  group_by(part_id) %>%
  subset(visit_num == 1) %>%
  mutate(time_pt = "Pre-Pandemic", .before = 1)
length(unique(data_res_pre$part_id))

## data of the most recent point after the post restriction
## note: used most recent visit to capture post restriction behaviour
data_res_post <- data_res %>% 
  filter(date_intv >= DATE_RES_END)  %>%
  group_by(part_id) %>%
  mutate(post_date = max(date_intv)) %>%
  subset(date_intv == post_date)%>%
  dplyr::select(-post_date) %>%
  mutate(time_pt = "Post-Restriction", .before = 1)
length(unique(data_res_post$part_id))

## data of the closest point after the pandemic
data_res_pand <- data_res %>% 
  filter(date_intv >= DATE_PAND_START & date_intv < DATE_RES_END) %>%
  group_by(part_id) %>%
  mutate(pand_date = min(date_intv)) %>%
  subset(date_intv == pand_date) %>%
  dplyr::select(-pand_date) %>%
  mutate(time_pt = "Pandemic", .before = 1)
length(unique(data_res_pand$part_id))
write_csv(data_res_pre,  "./data-3cities-feb-2023/restriction/res_pre.csv")
write_csv(data_res_pand,  "./data-3cities-feb-2023/restriction/res_pand.csv")
write_csv(data_res_post,  "./data-3cities-feb-2023/restriction/res_post.csv")

# verify that visits don't overlap
range(data_res_pre$date_intv)
range(data_res_pand$date_intv)
range(data_res_post$date_intv)

# save dataset
data_res_fu <- rbind(data_res_pre,data_res_pand, data_res_post)
table(visit = data_res_fu$time_pt, data_res_fu$city)

# save dataset
write_csv(data_res_fu,  "./data-3cities-feb-2023/restriction/res_fu.csv")

# data with LTFU analysis  -----------------------------------------------------
data_fu_pre <- data_fu %>% 
  group_by(part_id) %>%
  subset(visit_num == 1) %>% 
  mutate(time_pt = "Pre-Pandemic")
sum(data_fu_pre$wt_rds_norm)
length(unique(data_fu_pre$part_id))

data_fu_pand <- data_fu %>% 
  group_by(part_id) %>%
  subset(pand_vis == 1) %>% 
  filter(date_intv >= DATE_PAND_START & date_intv < DATE_RES_END) %>%
  group_by(part_id) %>%
  mutate(pand_date = min(date_intv)) %>%
  subset(date_intv == pand_date) %>%
  dplyr::select(-pand_date) %>% 
  mutate(time_pt = "Pandemic")

data_fu_post <- data_fu %>% 
  group_by(part_id) %>%
  subset(post_vis == 1) %>% 
  filter(date_intv >= DATE_RES_END) %>%
  group_by(part_id) %>%
  mutate(post_date = max(date_intv)) %>%
  subset(date_intv == post_date)%>%
  dplyr::select(-post_date)%>% 
  mutate(time_pt = "Post-Restriction")

data_visit_key_date <- rbind.data.frame(data_fu_pre,data_fu_pand,data_fu_post) %>%
  group_by(city,part_id) %>%
  mutate_at(c("bath_m", "bath_d", 
              "education_level_cat","income_level_cat", "ethnicity_cat", 
              "groupsex_m", "groupsex_d", 
              "sex_work_m", "sex_work_d"),
            as.factor)

# retention rate
table(data_visit_key_date$time_pt,
      data_visit_key_date$city)

data_visit_key_date %>%
  group_by(time_pt,city) %>%
  summarise(count = n(),
            normalized_rds_sum = sum(wt_rds_norm))

# full dataset analysis (add ipcw) --------------------------------------------------------------------
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
               "nb_part_ttl") #variables associated with ltfu

data_ipcw_pand <- data_fu_pre %>% #use baseline values of covariate
  mutate_at(covariates, as.factor) %>%
  mutate_at(covariates, as.integer) %>%
  mutate(ltfu = as.logical(ifelse(pand_vis == 1,0,1)))

data_ipcw_post<- data_fu_pre %>% #use baseline values of covariate
  mutate_at(covariates, as.factor) %>%
  mutate_at(covariates, as.integer) %>%
  mutate(ltfu = as.logical(ifelse(post_vis == 1,0,1)))

for (cur_city in CITIES){
  data<-data_ipcw_pand[data_ipcw_pand$city == cur_city,]
  covariates_include <- find.imbalance(data, covariates, weights = "wt_rds_norm")
  data$ipw_rds <- make.ipcw(data, covariates_include, 'wt_rds_norm', 1) # determine capping weight
  print(summary(data$ipw_rds))
  print(calculate.smd(data, covariates_include,'ipw_rds'))
  
  #check that rds x ipcw sum to rds adj. numbers of participants (both for ltfu and non-ltfu separately), i.e. make sure that weighing does not inflate the sample size
  data_ltfu <- data %>%
    subset(ltfu == 1)
  print(sum(data_ltfu$ipw_rds))
  print(sum(data_ltfu$wt_rds_norm))
  
  data_fu <- data %>%
    subset(ltfu == 0)
  print(sum(data_fu$ipw_rds))
  print(sum(data_fu$wt_rds_norm))
  
  ipcw_pand[[cur_city]]$part_id <- data$part_id
  ipcw_pand[[cur_city]]$ipw_rds <- data$ipw_rds
  
}

for (cur_city in CITIES){
  data <- data_ipcw_post[data_ipcw_post$city == cur_city, ]
  covariates_include <- find.imbalance(data, covariates, weights = "wt_rds_norm")
  data$ipw_rds <- make.ipcw(data, covariates_include, 'wt_rds_norm', 1) #determine capping weight
  
  print(summary(data$ipw_rds))
  print(calculate.smd(data, covariates_include, 'ipw_rds'))
  
  #check that rds x ipcw sum to rds adj. numbers of participants (both for ltfu and non-ltfu separately), i.e. make sure that weighing does not inflate the sample size
  data_ltfu<-data %>%
    subset(ltfu == 1)
  print(sum(data_ltfu$ipw_rds))
  print(sum(data_ltfu$wt_rds_norm))
  
  data_fu<-data %>%
    subset(ltfu == 0)
  print(sum(data_fu$ipw_rds))
  print(sum(data_fu$wt_rds_norm))
  
  ipcw_post[[cur_city]]$part_id <- data$part_id
  ipcw_post[[cur_city]]$ipw_rds <- data$ipw_rds
  
}
ipcw_pand <- do.call(rbind.data.frame, ipcw_pand)
ipcw_post <- do.call(rbind.data.frame, ipcw_post)
ipcw_pre <- data_visit_key_date %>%
  subset(time_pt == "Pre-Pandemic") %>%
  mutate(ipw_rds = wt_rds_norm)
ipcw_pand <- merge.data.frame(data_visit_key_date, ipcw_pand, by = "part_id") %>%
  subset(time_pt == "Pandemic")
ipcw_post <- merge.data.frame(data_visit_key_date, ipcw_post, by = "part_id") %>%
  subset(time_pt == "Post-Restriction")

# write.csv(ipcw_pre,"./data-3cities-feb-2023/pre_ipcw_3cities.csv")
# write.csv(ipcw_pand,"./data-3cities-feb-2023/pand_ipcw_3cities.csv")
# write.csv(ipcw_post,"./data-3cities-feb-2023/post_ipcw_3cities.csv")
