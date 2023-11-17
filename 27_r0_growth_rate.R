# Libraries and data ----
library(tidyverse)
library(ggtext)

theme_set(theme_bw())

## mpox parameters
# get infectious duration from fit model
df_pars <- read.csv("./out-seir/fit_pars.csv")
# infectious duration
D = subset(df_pars, city == "all" & parameter == "duration infectiousness (1/gamma)")$value
# latency duration
D_apostr = 5.1

## case data for R0 based on growth rate
#     from PHAC, https://health-infobase.canada.ca/mpox/
#     (use monkeypox-detailed.csv file)
case_data <- read.csv("./data-public/monkeypox-detailed-2023june13.csv")
case_data$reporting_pt_en[case_data$reporting_pt_en == "Quebec"] <- "Québec"

# only need the 3 concerned provinces
PROVS <- c("Québec", "Ontario", "British Columbia")
case_data <- case_data %>%
  subset(reporting_pt_en %in% PROVS)

CITIES <- c("Montreal", "Toronto", "Vancouver")

# R0 from cases ----
## format data, to get time from first case in each province and ln(cumulative cases)
case_data <- case_data %>% 
  group_by(reporting_pt_en) %>%
  mutate(
    earliest_date = as.Date(min(date)),
    time_conti = as.integer(as.Date(date) - earliest_date),
    ln_cumu_cases = log(num_confirmedcases_cumulative),
    ln_incid = log(num_confirmedcases_delta),
    reporting_pt_en = factor(reporting_pt_en, level = PROVS)
  )

## Main results (first 50 days) ----
# cutoff data to use in growth rate estimation for R0
cutoff <- 50
case_data_cutoff <- case_data %>% 
  filter(time_conti <= cutoff & reporting_pt_en == "British Columbia" |
           time_conti <= cutoff & reporting_pt_en == "Ontario" |
           time_conti <= cutoff & reporting_pt_en == "Québec")

## plot cumulative incidence and fitted log-growth rate
png("./fig/fig_S2_cumul_incidence.png",
    width = 15, height = 7, units = "cm", res = 600)
ggplot(case_data, aes(x = time_conti,
                      y = ln_cumu_cases,
                      group = reporting_pt_en,
                      col = reporting_pt_en)) +
  geom_point(size = 0.6, alpha = 0.8) + 
  # geom_line(aes(y = ln_incid), linetype = "dashed") + # can use to compare growth in cumulative cases vs incidence
  # plot a straight line for the data used in the regression
  # this is equivalent to doing lm(ln_cumu_cases ~ time_conti) for each province
  geom_smooth(data = case_data_cutoff, method = "lm", se = FALSE, linewidth = 0.5) +
  # change axis titles
  
  coord_cartesian(xlim = c(0, 200)) +
  scale_colour_viridis_d(option = "C", end = 0.8)  +
  labs(x = "Days since first case",
       y = "ln(Cumulative confirmed mpox cases)",
       col = "Province") +
  
  theme(axis.title = element_text(size = 9),
        axis.text = element_text(size = 7),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 7))
dev.off()

lambda <- vector(mode = "numeric", 
                 length = 3)
names(lambda) <- PROVS
for (province in PROVS){
  case_data_province <- case_data_cutoff %>% 
    subset(reporting_pt_en == province)
  
  # regression
  linear_model <- lm(formula = ln_cumu_cases ~ time_conti, case_data_province)
  # linear_model <- lm(formula = incid ~ time_conti, case_data_province)
  
  # growth rate
  lambda[[province]] <- linear_model$coefficients[2] %>%
    unname()
}

# formula is (1 + lambda*D)(1 + lambda*D')
r0_case <- (1 + lambda * D) * (1 + lambda * D_apostr)
round(r0_case, 2)

r0_main <- data.frame(time_cutoff = rep(cutoff, 3), prov = names(r0_case), r0 = r0_case)

## Sensitivity (first 40 or 60 days) ----
r0_sens <- vector("list", 2)
sens_cutoffs <- c(40, 60)

for(cur_cutoff in sens_cutoffs) {
  # subset case data
  case_tmp <- case_data %>% 
    filter(time_conti <= cur_cutoff & reporting_pt_en == "British Columbia" |
             time_conti <= cur_cutoff & reporting_pt_en == "Ontario" |
             time_conti <= cur_cutoff & reporting_pt_en == "Québec")
  
  # compute growth rate and R0 in each province
  lambda <- vector(mode = "numeric", 
                   length = 3)
  names(lambda) <- PROVS
  
  for (province in PROVS){
    case_tmp_province <- case_tmp %>% 
      subset(reporting_pt_en == province)
    
    linear_model <- lm(formula = ln_cumu_cases ~ time_conti, case_tmp_province)
    
    lambda[[province]] <- linear_model$coefficients[2] %>%
      unname()
  }
  
  r0_tmp <- (1 + lambda * D) * (1 + lambda * D_apostr)
  r0_sens[[which(cur_cutoff == sens_cutoffs)]] <- data.frame(
    time_cutoff = rep(cur_cutoff, 3), prov = names(r0_tmp), r0 = r0_tmp
  )
}

r0_sens

## output results ----
# put outputs together
r0_out <- bind_rows(r0_main, r0_sens)

# format
r0_out <- r0_out %>% pivot_wider(names_from = prov, values_from = r0)
r0_out <- r0_out %>% arrange(time_cutoff)

write.csv(r0_out, "./out/r0-estim-cases.csv", row.names = FALSE)

