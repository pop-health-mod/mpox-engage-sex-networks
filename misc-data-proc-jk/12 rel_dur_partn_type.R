# Library and data----
library(tidyverse)

## main analyses
outcome_var <- "nb_part_ttl"

## define paths & prefixes based on the analysis being done
fig_path <- "./misc-data-proc-JK/outputs"
out_distr_path <- "./misc-data-proc-JK/outputs"

## load data
data_3cities_pre_ipcw <- read_csv("../mpx-engage-params/data-processed/pre_ipcw_3cities.csv")

## proceed with grouping by relationship status
# reorder by turning into a factor
data_base <- data_base %>% 
  mutate(reltn_type = factor(
    reltn_type,
    c("one-time", "casual", "main", "sex work", "PNA"))
  )

# exclude 3 individuals who didn't report any relations
data_base <- filter(data_base, !is.na(partner_id))

# exclude 257 relationships that happened >6 months ago (6% of entries)
table(data_base$reltn_too_long_ago, useNA = "ifany")

# exclude entries for implausible relationship dates (5%)
table(data_base$implausible_reltn, useNA = "ifany")

data_base <- filter(data_base, !reltn_too_long_ago | is.na(reltn_too_long_ago))
data_base <- filter(data_base, is.na(implausible_reltn))

tapply(data_base$par_duration, data_base$reltn_type, summary)
tapply(data_base$par_duration/30, data_base$reltn_type, summary)

# Descriptive analyses ----
## Number of partners reported ----
# tabulate number of respondents in each relationship status
data_nb_partn <- data_demo %>% 
  mutate() %>% 
  group_by(rel_status) %>% 
  summarize(N_total = n(),
            wt_total = sum(wt))

# tabulate total number of partners reported by relationship status and partner type
tmp_nb_partn <- data_base %>% 
  count(rel_status, reltn_type, name = "nb_partn_rep")

data_nb_partn <- full_join(data_nb_partn, tmp_nb_partn, by = "rel_status")

# compute mean number of partners and prepare output table
data_nb_partn <- data_nb_partn %>% 
  mutate(mean_partn_rep = nb_partn_rep / N_total)

nb_partn_out <- data_nb_partn %>% 
  select(rel_status, reltn_type, mean_partn_rep) %>% 
  pivot_wider(names_from = rel_status, values_from = mean_partn_rep)

# add overall (not rel_status-stratified) mean number of partners per category
tmp_nb_partn <- tmp_nb_partn %>% 
  group_by(reltn_type) %>% 
  summarize(nb_reported = sum(nb_partn_rep))
tmp_nb_partn$overall <- tmp_nb_partn$nb_reported / n_respondent

nb_partn_out <- full_join(tmp_nb_partn, nb_partn_out, by = c("reltn_type"))
rm(tmp_nb_partn)

## verify that column mean matches with mean reported number 
# overall
nb_partn_out[, 3] %>% colSums(na.rm = T)
nrow(data_base) / nrow(data_demo)

# by relationship status
nb_partn_out[, 4:7] %>% colSums(na.rm = T)
data_base %>% 
  count(rel_status, name = "nb_partn_rep") %>% 
  left_join(count(data_demo, rel_status, name = "n_respondents"),
            by = c("rel_status")) %>% 
  mutate(mean = nb_partn_rep / n_respondents)

write.csv(nb_partn_out, "./output/partn-duration-freq/mean_nb_partn.csv",
          row.names = FALSE)

## Duration of partnerships ----
### sheer number
ggplot(data_base, aes(x = par_duration/30, y = reltn_type, fill = reltn_type)) +
  # 6-month limit
  geom_vline(xintercept = 180/30) +
  geom_boxplot() +
  coord_cartesian(xlim = c(0, 5000)/30) +
  scale_x_continuous("Partnership duration in months") +
  scale_y_discrete(limits = rev(levels(data_base$reltn_type))) +
  ggtitle("Partnership duration (months), by Type of Partnership") +
  theme(legend.position = "none")

ggsave("./figures/partn-sex-freq/distr_partn_length_by_type.png", device = "png",
       width = 30, height = 14, units = "cm", dpi = 300)

ggplot(data_base, aes(x = par_duration, group = reltn_type, fill = reltn_type)) +
  geom_histogram(bins = 100) +
  coord_cartesian(xlim = c(0, 5000)) +
  scale_x_continuous("Partnership duration in days") +
  ggtitle("Partnership duration in days, by Type of Partnership") +
  scale_fill_discrete(name = "type of partnership")

### faceted by reltn_type
duration_by_type <- data_base %>% 
  count(reltn_type, par_duration) %>%
  group_by(reltn_type) 

# ggplot(data_base, aes(x = par_duration)) +
#   geom_histogram(y = n, bins = 100) + facet_wrap(~reltn_type) +
#   scale_fill_viridis_d() +
#   coord_cartesian(xlim = c(0, 5000)) +
#   scale_x_continuous("Partnership duration in days") +
#   ylab("Count") +
#   ggtitle("Partnership duration in days for Partnership Types")

## faceted by rel_status
ggplot(data_base,
       aes(x = par_duration/30, y = reltn_type, fill = reltn_type)) +
  # 6-month limit
  geom_vline(xintercept = 180/30) +
  geom_boxplot() +
  
  facet_wrap(~rel_status) +
  
  coord_cartesian(xlim = c(0, 5000)/30) +
  scale_x_continuous("Partnership duration in months") +
  scale_y_discrete(limits = rev(levels(data_base$reltn_type))) +
  ggtitle("Partnership duration (months), by Type of Partnership") +
  theme(legend.position = "none")

ggsave("./figures/partn-sex-freq/distr_partn_length_by_type_relstatus.png", device = "png",
       width = 30, height = 14, units = "cm", dpi = 300)

### Output partnership duration ----
# compute raw duration
data_duration <- data_base %>% 
  mutate(par_duration = ifelse(par_duration == 0, 1, par_duration)) %>% 
  group_by(reltn_type) %>% 
  summarize(
    across(par_duration, list(mean=~mean(., na.rm = TRUE), sd = ~sd(., na.rm = TRUE),
                              Q0 = ~min(., na.rm = TRUE), Q1 = ~quantile(.,0.25,na.rm = TRUE), Q2 = ~median(., na.rm = TRUE),
                              Q3 = ~quantile(.,0.75,na.rm = TRUE), Q4 = ~max(., na.rm = TRUE)),
           .names = "dur_{.fn}"))

data_duration_100quant <- data_base %>% 
  mutate(par_duration = ifelse(par_duration == 0, 1, par_duration)) %>% 
  group_by(reltn_type) %>% 
  summarize(q=list(quantile(par_duration, probs = seq(0,1,0.01), na.rm=T)))%>%
  unnest_wider(q, names_repair = ~paste0('dur_Q', sub('%', '', .)))

# compute duration setting ceiling to 90% of data
# 90% of entries last <1642 days (4.5 years), 99% <8802 days (24 yrs)
quantile(data_base$par_duration, c(.9, .99), na.rm = T)
tapply(data_base$par_duration, data_base$reltn_type, na.rm = T, quantile, c(.9, .99))

data_duration_adj <- data_base %>% 
  mutate(par_duration = case_when(
    par_duration == 0 ~ 1,
    par_duration >= 365*5 ~ 365*5,
    TRUE ~ par_duration
    )) %>% 
  group_by(reltn_type) %>% 
  summarize(
    across(par_duration, list(mean=~mean(., na.rm = TRUE), sd = ~sd(., na.rm = TRUE),
                              Q0 = ~min(., na.rm = TRUE), Q1 = ~quantile(.,0.25,na.rm = TRUE), Q2 = ~median(., na.rm = TRUE),
                              Q3 = ~quantile(.,0.75,na.rm = TRUE), Q4 = ~max(., na.rm = TRUE)),
           .names = "dur_{.fn}")
  )

data_duration_adj_100quant <- data_base %>% 
  mutate(par_duration = case_when(
    par_duration == 0 ~ 1,
    par_duration >= 365*5 ~ 365*5,
    TRUE ~ par_duration
  )) %>% 
  group_by(reltn_type) %>% 
  summarize(q=list(quantile(par_duration, probs = seq(0,1,0.01), na.rm=T)))%>%
  unnest_wider(q)


write.csv(data_duration, "./output/partn-duration-freq/partnership_duration_obs.csv",
          row.names = FALSE)
write.csv(data_duration_adj, "./output/partn-duration-freq/partnership_duration_adjusted.csv",
          row.names = FALSE)

write.csv(data_duration_100quant, "./output/partn-duration-freq/partnership_duration_obs_quantiles.csv",
          row.names = FALSE)
write.csv(data_duration_adj_100quant, "./output/partn-duration-freq/partnership_duration_adjusted_quantiles.csv",
          row.names = FALSE)
