# Libraries and data ----
library(tidyverse)
library(ggtext)

theme_set(theme_bw())

## case data for R0 based on growth rate
#     from PHAC, https://health-infobase.canada.ca/mpox/
#     (use monkeypox-detailed.csv file)
case_data <- read.csv("../mpx-engage-params/misc-grant-app/data-case/monkeypox-detailed-2023june13.csv")
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
  linear_model <- lm(formula = ln_cumu_cases ~ time_conti, case_data_province)
  # linear_model <- lm(formula = incid ~ time_conti, case_data_province)
  lambda[[province]] <- linear_model$coefficients[2] %>%
    unname()
}

# formula is (1 + lambda*D)(1 + lambda*D') where D: infectious period & D': 
SI = 8.99
D_apostr = (7.12 - 2.5)
D = SI - D_apostr

r0_case <- (1 + lambda * D) * (1 + lambda * D_apostr)
round(r0_case, 2)

# simplified formula, should roughly correspond but use above in main paper
# round( (1 + lambda * SI) , 2)

write.csv(r0_case, "./out/r0-estim-cases.csv", row.names = FALSE)
