# Library and data----
library(tidyverse)
library(lubridate)

## define paths & prefixes based on the analysis being done
fig_path <- "./misc-data-proc-jk/outputs"
out_distr_path <- "./misc-data-proc-jk/outputs"

## load data
data_3cities <- read_csv("../mpx-engage-params/data-processed/engage_baseline_3cities.csv")
data_partn <- read_csv("../mpx-engage-params/data-processed/engage_partners_3cities.csv")
data_partn <- data_partn[data_partn$visit_num == 1, ] # N=9,043

# create city marker
CITIES <- sort( unique(data_3cities$city) )

## proceed with grouping by relationship status
# reorder by turning into a factor
unique(data_partn$reltn_type)
data_partn$reltn_type[data_partn$reltn_type == "sex work"] <- "sex-work"
data_partn <- data_partn %>% 
  mutate(reltn_type = factor(
    reltn_type,
    c("one-time", "casual", "main", "sex-work"))
  )

RELTN_TYPES <- levels(data_partn$reltn_type)

data_partn <- left_join(data_partn,
                        data_3cities[, c("city", "part_id", "date_intv",
                                         "wt_rds", "wt_rds_norm")],
                        by = c("part_id")) %>% 
  select(city, part_id, visit_num, date_intv, wt_rds, wt_rds_norm,
         partner_id:p6m_sex_num_pna)

## Data cleaning - partnership duration ----
data_partn <- mutate(data_partn,
                     rel_length = as.numeric(recent_sex_date - first_sex_date),
                     .after = recent_sex_date)

# fix NAs (n=44)
data_partn[is.na(data_partn$rel_length), ]
table(data_partn$reltn_type[is.na(data_partn$rel_length)])
data_partn$rel_length[is.na(data_partn$rel_length)] <- 1

# fix relationship lengths caused by not having day for first_sex (n=2,156)
table(data_partn$rel_length <= 0)
table(data_partn$rel_length <= 0 & 
        year(data_partn$first_sex_date) == year(data_partn$recent_sex_date) & 
        month(data_partn$first_sex_date) == month(data_partn$recent_sex_date))
data_partn$rel_length[data_partn$rel_length <= 0 & 
                        year(data_partn$first_sex_date) == year(data_partn$recent_sex_date) & 
                        month(data_partn$first_sex_date) == month(data_partn$recent_sex_date)] <- 1

# assume if month is off by just one same situation as above (n=126)
table(data_partn$rel_length <= 0)
table(data_partn$rel_length <= 0 & 
        # both within same year
        (year(data_partn$first_sex_date) == year(data_partn$recent_sex_date) & 
        month(data_partn$first_sex_date) == (month(data_partn$recent_sex_date) + 1)) |
        # first recorded as year Jan and recent happened in (year-1) Dec
        ((year(data_partn$first_sex_date) - 1) == year(data_partn$recent_sex_date) & 
           month(data_partn$first_sex_date) == (month(data_partn$recent_sex_date) == 12)))
data_partn$rel_length[data_partn$rel_length <= 0 & 
                        (year(data_partn$first_sex_date) == year(data_partn$recent_sex_date) & 
                           month(data_partn$first_sex_date) == (month(data_partn$recent_sex_date) + 1)) |
                        ((year(data_partn$first_sex_date) - 1) == year(data_partn$recent_sex_date) & 
                           month(data_partn$first_sex_date) == (month(data_partn$recent_sex_date) == 12))] <- 1

# assume typo if month is same but year off by 1, if month matches then assume recent=first (n=34)
table(data_partn$rel_length <= 0)
table(data_partn$rel_length <= 0 & 
        (year(data_partn$first_sex_date) - 1) == year(data_partn$recent_sex_date) & 
        month(data_partn$first_sex_date) == month(data_partn$recent_sex_date))
data_partn$rel_length[data_partn$rel_length <= 0 & 
                        (year(data_partn$first_sex_date) - 1) == year(data_partn$recent_sex_date) & 
                        month(data_partn$first_sex_date) == month(data_partn$recent_sex_date)] <- 1



# assume typo if month is same but year off by 1, if recent month is after then use that (n=23)
table(data_partn$rel_length <= 0)
table(data_partn$rel_length <= 0 & 
        (year(data_partn$first_sex_date) - 1) == year(data_partn$recent_sex_date) & 
        month(data_partn$first_sex_date) <= month(data_partn$recent_sex_date))
data_partn <- mutate(
  data_partn, rel_length = ifelse(rel_length <= 0 & 
                                    (year(first_sex_date) - 1) == year(recent_sex_date) & 
                                    month(first_sex_date) <= month(recent_sex_date),
                                  as.numeric(recent_sex_date - (first_sex_date - years(1))),
                                  rel_length
  )
)

# assume if one-time or casual then duration of 1 (n = 59)
table(data_partn$reltn_type[data_partn$rel_length <= 0])
data_partn$rel_length[data_partn$rel_length <= 0 & 
                        (data_partn$reltn_type == "one-time" | data_partn$reltn_type == "casual")] <- 1

# remaining 12 impute 1 as well
table(data_partn$rel_length <= 0)
data_partn[data_partn$rel_length <= 0, ]

data_partn$rel_length[data_partn$rel_length <= 0] <- 1

# simple descriptive tables
tapply(data_partn$rel_length, paste(data_partn$city, data_partn$reltn_type, sep = "-"), summary)
tapply(data_partn$rel_length/30, paste(data_partn$city, data_partn$reltn_type, sep = "-"), summary)

## Data cleaning - sex frequency & rate ----
# exclude PNA's (n=228)
table(data_partn$p6m_sex_num_pna, useNA = "ifany")
data_partn_freq <- filter(data_partn, p6m_sex_num_pna == 0)
sum(is.na(data_partn_freq$p6m_sex_num))

# simple descriptive tables
tapply(data_partn_freq$p6m_sex_num, paste(data_partn_freq$city, data_partn_freq$reltn_type, sep = "-"), summary)

# assume cap at 1 per day (n=25)
sum(data_partn_freq$p6m_sex_num > 180)
data_partn_freq$p6m_sex_num <- ifelse(data_partn_freq$p6m_sex_num > 180,
                                      180, data_partn_freq$p6m_sex_num)

# check reports of nb sex acts > min(length of relationship, 6 months) (n = 757)
data_partn_freq <- data_partn_freq %>% 
  group_by(part_id, visit_num, partner_id) %>% 
  mutate(rel_length_cap = min(rel_length, 180), .after = rel_length) %>% 
  ungroup()

sum(data_partn_freq$p6m_sex_num > data_partn_freq$rel_length_cap)

# impute minimum as 1 given that person was reported as partner (n = 200)
sum(data_partn_freq$p6m_sex_num == 0)
data_partn_freq$p6m_sex_num <- ifelse(data_partn_freq$p6m_sex_num == 0,
                                      1, data_partn_freq$p6m_sex_num)

# Descriptive -----
## Partnership duration ----
# visualize
ggplot(data_partn, aes(x = city, y = rel_length, fill = city)) +
  geom_violin(alpha = 0.5) +
  
  coord_cartesian(ylim = c(0, 2000)) +
  facet_wrap(~reltn_type)

## summarize
# by key statistics
data_partn_num_tbl <- data_partn %>% 
  group_by(city, reltn_type) %>% 
  summarize(
    n_indv = length(unique(part_id)), n_rel = n(),
    mean = mean(rel_length), sd = sd(rel_length),
    q0 = min(rel_length), q1 = quantile(rel_length, .25),
    q2 = median(rel_length), q3 = quantile(rel_length, .75),
    q4 = max(rel_length), .groups = "drop"
  )

# by 1-percentile
quants <- seq(0, 1, .01)
data_partn_num_perc <- data.frame()

for(cur_cty in CITIES){
  for(cur_type in RELTN_TYPES){
    data_tmp <- filter(data_partn, city == cur_cty & cur_type == reltn_type)
    
    # compute summary statistics
    data_summ_tmp <- data.frame(
      city = cur_cty, reltn_type = cur_type,
      n_indv = length(unique(data_tmp$part_id)), n_rel = nrow(data_tmp),
      mean = mean(data_tmp$rel_length), sd = sd(data_tmp$rel_length)
    )
    
    # compute percentiles
    data_perc_tmp <- t(quantile(data_tmp$rel_length, probs = quants))
    colnames(data_perc_tmp) <- paste0("q_", format(quants, digits = 3))
    
    # join
    data_partn_num_perc <- rbind(data_partn_num_perc, cbind(data_summ_tmp, data_perc_tmp))
  }
}

rm(data_tmp, cur_cty, cur_type)

# verify results
data_partn_num_tbl
data_partn_num_perc[1:12, c(1:6, which(colnames(data_partn_num_perc) %in% c("q_0.00", "q_0.25", "q_0.50", "q_0.75", "q_1.00")))]

## Sex frequency & rate ----
data_partn_freq <- mutate(data_partn_freq, p6m_sex_rate = p6m_sex_num / rel_length_cap, .after = p6m_sex_num)

data_partn_freq <- data_partn_freq %>% 
  group_by(part_id, visit_num, partner_id) %>% 
  mutate(p6m_sex_rate_cap = min(p6m_sex_rate, 1), .after = p6m_sex_rate) %>% 
  ungroup()

# visualize (frequency)
ggplot(data_partn_freq, aes(x = city, y = p6m_sex_num, fill = city)) +
  geom_violin(alpha = 0.5) +
  
  coord_cartesian(ylim = c(0, 200)) +
  facet_wrap(~reltn_type)

# visualize (rate, both uncapped and capped)
ggplot(data_partn_freq, aes(x = city, y = p6m_sex_rate, fill = city)) +
  geom_jitter(alpha = 0.2) +
  geom_violin(alpha = 0.5) +
  
  coord_cartesian(ylim = c(0, 2)) +
  facet_wrap(~reltn_type)

## summarize
# by key statistics
data_partn_rate_tbl <- data_partn_freq %>% 
  group_by(city, reltn_type) %>% 
  summarize(
    n_indv = length(unique(part_id)), n_rel = n(),
    mean = mean(p6m_sex_rate_cap), sd = sd(p6m_sex_rate_cap),
    q0 = min(p6m_sex_rate_cap), q1 = quantile(p6m_sex_rate_cap, .25),
    q2 = median(p6m_sex_rate_cap), q3 = quantile(p6m_sex_rate_cap, .75),
    q4 = max(p6m_sex_rate_cap), .groups = "drop"
  )

# by 1-percentile
quants <- seq(0, 1, .01)
data_partn_rate_perc <- data.frame()

for(cur_cty in CITIES){
  for(cur_type in RELTN_TYPES){
    data_tmp <- filter(data_partn_freq, city == cur_cty & cur_type == reltn_type)
    
    # compute summary statistics
    data_summ_tmp <- data.frame(
      city = cur_cty, reltn_type = cur_type,
      n_indv = length(unique(data_tmp$part_id)), n_rel = nrow(data_tmp),
      mean = mean(data_tmp$p6m_sex_rate_cap), sd = sd(data_tmp$p6m_sex_rate_cap)
    )
    
    # compute percentiles
    data_perc_tmp <- t(quantile(data_tmp$p6m_sex_rate_cap, probs = quants))
    colnames(data_perc_tmp) <- paste0("q_", format(quants, digits = 3))
    
    # join
    data_partn_rate_perc <- rbind(data_partn_rate_perc, cbind(data_summ_tmp, data_perc_tmp))
  }
}

rm(data_tmp, cur_cty, cur_type)

# verify results
data_partn_rate_tbl
data_partn_rate_perc[1:12, c(1:6, which(colnames(data_partn_rate_perc) %in% c("q_0.00", "q_0.25", "q_0.50", "q_0.75", "q_1.00")))]

# Output ----

write.csv(data_partn_rate_perc, "./misc-data-proc/outputs-jk/engage_partnership_sexrate.csv",
          row.names = FALSE)
write.csv(data_partn_num_perc, "./misc-data-proc/outputs-jk/engage_partnership_duration.csv",
          row.names = FALSE)



