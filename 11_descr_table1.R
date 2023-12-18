# Library and data -------------
library(tidyverse)
source("./src/utils_helper.R")
data_3cities <- read_csv("../mpx-engage-params/data-processed/engage_baseline_3cities.csv")

# Data cleaning and formating ------------------
## create dummy variables
# age; reference is 18-29
data_3cities <- make_ind_age(data_3cities)
table(data_3cities$city, data_3cities$age_grp, useNA = "ifany")

# partnership status
data_3cities <- make_ind_rel(data_3cities)
table(data_3cities$city, data_3cities$rel_status, useNA = "ifany")

# check SPVs and groupsex
table(data_3cities$city, data_3cities$bath, useNA = "ifany")
table(data_3cities$city, data_3cities$groupsex, useNA = "ifany")

data_3cities$bath[is.na(data_3cities$bath)] <- "Missing"
data_3cities$groupsex[is.na(data_3cities$groupsex)] <- "Missing"

# format apps variable
table(data_3cities$city, data_3cities$apps_partn, useNA = "ifany")

# format sex work variable (assume NAs are no)
table(data_3cities$city, data_3cities$sex_work, useNA = "ifany")

# hiv status 
table(data_3cities$city,data_3cities$hiv_stat, useNA = "ifany")

# turn into factors to ensure correct ordering and coding
data_3cities <- data_3cities %>%
  mutate(
    rel_status = factor(rel_status,
                        levels = c("no relationship", "open", "exclusive", "unclear"),
                        labels = c("Single", "Open", "Exclusive", "Unclear")),
    sex_work = factor(sex_work, levels = c("yes", "no", "PNA"),
                      labels = c("Yes", "No", "Missing")),
    hiv_stat = factor(hiv_stat,
                      levels = c(1, 0),
                      labels = c("Seropositive", "Seronegative")),
    groupsex = factor(groupsex,
                      levels = c(1, 0, "Missing"),
                      labels = c("Yes", "No", "Missing")),
    bath = factor(bath,
                  levels = c(1, 0, "Missing"),
                  labels = c("Yes", "No", "Missing")),
    apps_partn = factor(ifelse(is.na(apps_partn), "Missing", apps_partn),
                        levels = c(1, 0, "NA"),
                        labels = c("Yes", "No", "Missing"))
  )

# Table 1, unadjusted & RDS (pre-pandemic only) ----
## Table for proportions ----
# duplicate data to have outputs for 'overall' (the 3 cities combined)
data_dbl <- bind_rows(
  data_3cities,
  mutate(data_3cities, city = "all")
)

# generate variable of total nb of participants
data_dbl <- data_dbl %>% 
  group_by(city) %>% 
  mutate(city_n = n()) %>% 
  ungroup()

# city
tbl_1 <- data_dbl %>% 
  group_by(city) %>% 
  summarize(nb = n(), .groups = "drop")

ls_var <- c("age_grp", "rel_status", "hiv_stat", "bath",
            "groupsex", "apps_partn", "sex_work")

# for each categorical variable compute the proportions
for(var in ls_var){
  tbl_1 <- tbl_1 %>% 
    bind_rows(
      data_dbl %>% 
        group_by(city_n, city, response = get(var)) %>% 
        summarize(nb = n(), nb_rds = sum(wt_rds_norm),
                  .groups = "drop_last") %>% 
        mutate(
          # unadjusted proportion
          prop = nb / city_n,
          # RDS-adjusted
          prop_rds = nb_rds / city_n,
          # RDS 95% CI
          me = sqrt(prop_rds*(1 - prop_rds))/sqrt(city_n) * 1.96,
          ci.lb = prop_rds - me,                                   
          ci.ub = prop_rds + me,
          # variable label
          char = var
        )
    )
}

## Table for mean number of partners ----
ls_vars_partn <- c("nb_part_ttl", "nb_part_anal")
tbl_1_partn <- tibble()

for(var in ls_vars_partn){
  tbl_1_partn <- tbl_1_partn %>% 
    bind_rows(
      data_dbl %>%
        mutate(nb_part_rds = get(var) * wt_rds_norm) %>% 
        group_by(city) %>% 
        summarize(
          # unadjusted mean
          mean_nb = mean(get(var)),
          sd_nb = sprintf("(SD=%s)", round(sd(get(var)), 1)),
          # RDS-adjusted
          city_n = n(),
          mean_rds = sum(nb_part_rds) / city_n,
          # RDS 95% CI
          me = sqrt( sum((nb_part_rds - mean_rds)^2) ) / city_n,
          ci.lb = mean_rds - me * 1.96,
          ci.ub = mean_rds + me * 1.96,
          .groups = "drop"
        ) %>%
        mutate(char = var)
    )
}

## Pivot and join tables ----
## proportions
# round proportions (need to turn into character to join with SD later)
tbl_1 <- tbl_1 %>% 
  mutate(across(c(prop, prop_rds, ci.lb, ci.ub), \(x) round_prop(x)))

# format for table
tbl_1 <- tbl_1 %>% 
  mutate(rds_ci = sprintf("(%s\u2013%s)", ci.lb, ci.ub))

# switch to one column per city
tbl_1 <- tbl_1 %>% 
  select(city, nb, response, prop, prop_rds, char, rds_ci) %>% 
  pivot_wider(names_from = "city", values_from = c("nb", "prop", "prop_rds", "rds_ci"))


## number of patners
# round
tbl_1_partn <- tbl_1_partn %>% 
  mutate(across(c(mean_nb, mean_rds, ci.lb, ci.ub), \(x) round(x, 1)))

# turn RDS mean into character to join later
tbl_1_partn <- tbl_1_partn %>% 
  mutate(mean_rds = as.character(mean_rds))

# format for table
tbl_1_partn <- tbl_1_partn %>% 
  mutate(rds_ci = sprintf("(%s\u2013%s)", ci.lb, ci.ub))

# switch to one column per city
tbl_1_partn <- tbl_1_partn %>% 
  select(city, mean_nb, sd_nb, mean_rds, rds_ci, char) %>% 
  pivot_wider(names_from = "city", values_from = c("mean_nb", "sd_nb", "mean_rds", "rds_ci"))

# rename columns (match mean with counts and SD with proportions)
names(tbl_1_partn) <- gsub("mean_nb", "nb", names(tbl_1_partn))
names(tbl_1_partn) <- gsub("mean_rds", "prop_rds", names(tbl_1_partn))
names(tbl_1_partn) <- gsub("sd_nb_", "prop_", names(tbl_1_partn))

## join together and reorder columns
tbl_1 <- bind_rows(tbl_1, tbl_1_partn)

## Output table ----
tbl_1 <- tbl_1 %>% 
  select(char, response,
         ends_with("_mtl"), ends_with("_trt"), ends_with("_van"), ends_with("all"))

write.csv(tbl_1, "./out/manuscript-tables/table_1_unadj_rds.csv", row.names = FALSE)
