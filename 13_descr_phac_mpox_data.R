# Libraries and data ----
library(tidyverse)

source("./src/plot.R")

theme_set(theme_bw())

## case data for model fitting
#     from PHAC, https://health-infobase.canada.ca/mpox/
#     (use monkeypox-detailed.csv file)
case_data <- read.csv("./data-public/monkeypox-detailed-2023june13.csv")
case_data$date <- as.Date(case_data$date)

# keep only 2022 data
case_data <- case_data %>% filter(date <= "2022-12-15")

# only need the 3 concerned provinces
PROVS <- c("Quebec", "Ontario", "British Columbia")
case_data <- case_data %>%
  subset(reporting_pt_en %in% PROVS)

# rename provinces
case_data <- case_data %>% 
  mutate(prov = case_when(reporting_pt_en == "Quebec" ~ "QC",
                          reporting_pt_en == "Ontario" ~ "ON",
                          reporting_pt_en == "British Columbia" ~ "BC"))

# format data, to get time from first case in each province and rename variables
case_data <- case_data %>% 
  group_by(reporting_pt_en) %>%
  mutate(
    earliest_date = as.Date(min(date)),
    time_conti = as.integer(as.Date(date) - earliest_date),
    incidence = num_confirmedcases_delta,
    cumul_cases = num_confirmedcases_cumulative,
  ) %>% 
  ungroup()

# drop unnecessary variables
case_data <- case_data %>% 
  select(prov, date, time_conti, incidence, cumul_cases)

## visualize data
# incidence
plot_phac_case(case_data) +
  geom_point()

# cumulative cases
plot_phac_case(case_data, yvar = "cumul_cases") +
  geom_line()

## pad with 0's for days without cases
make_pad_df <- function(df, province = "QC"){
  df <- df[df$prov == province, ]
  
  # create padding for dates
  df_pad <- expand.grid(prov = province,
                        date = seq(min(df$date), max(df$date), by = "1 day"))
  
  # add continuous time
  df_pad$time_conti <- 0:(nrow(df_pad) - 1)
  
  return(df_pad)
}

## add 0-incidence rows to each province
# quebec
df_pad_qc <- make_pad_df(case_data, "QC")
case_data <- case_data %>% full_join(df_pad_qc, by = c("prov", "date", "time_conti"))

# ontario
df_pad_on <- make_pad_df(case_data, "ON")
case_data <- case_data %>% full_join(df_pad_on, by = c("prov", "date", "time_conti"))

# british columbia
df_pad_bc <- make_pad_df(case_data, "BC")
case_data <- case_data %>% full_join(df_pad_bc, by = c("prov", "date", "time_conti"))

case_data <- case_data %>% arrange(prov, time_conti)

# pad incidence
case_data$incidence[is.na(case_data$incidence)] <- 0

# pad cumulative
case_data <- case_data %>% 
  group_by(prov) %>% 
  fill(cumul_cases, .direction = "down") %>% 
  ungroup()

case_data <- arrange(case_data, prov, date)

write_csv(case_data, "./data-public/mpox_case_data_phac.csv")

## visualize data
# incidence
plot_phac_case(case_data) +
  geom_point()

# cumulative cases
plot_phac_case(case_data, yvar = "cumul_cases") +
  geom_line()

