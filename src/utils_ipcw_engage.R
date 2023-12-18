################################################################################
# This code uses inverse probability of censoring weights to correct for 
# dependent censoring in RDS studies. 
# Author: Alain Fourmigue
# Created On: 2023-11-23
# © 2023 Alain Fourmigue <alain.fourmigue@free.fr>
###############################################################################

# calculate means and std. diff. mean (smd) of binary variables between (non)ltfu participants
# data: data frame containing binary variables
# ltfu: a binary variable indicating the participants who were ltfu
# wgt: an optional numerical vector of weights (e.g. ipcw weights)
calc.smd <- function(data, ltfu, wgt = rep(1., nrow(data))) {
  stopifnot(all(ltfu == 0 | ltfu == 1))
  n.0 <- sum(1 - ltfu)
  n.1 <- sum(ltfu)
  prop.0 <- sapply(subset(data, ltfu == 0), weighted.mean, wgt[ltfu == 0], na.rm = T)
  prop.1 <- sapply(subset(data, ltfu == 1), weighted.mean, wgt[ltfu == 1], na.rm = T)
  ss.0 <- n.0 * prop.0 * (1 - prop.0)
  ss.1 <- n.1 * prop.1 * (1 - prop.1)
  sd.pooled <- sqrt((ss.0 + ss.1) / (n.0 + n.1 - 2))
  smd <- abs(prop.0 - prop.1) / sd.pooled
  round(cbind(non.ltfu = prop.0, ltfu = prop.1, smd = smd), 2)
}

# calculate inverse probability of censoring weights (ipcw)
# data: data frame containing the predictors used to model ltfu
# ltfu: a binary variable indicating the participants who were ltfu
# prds: vector of strings indicating the predictors used to model ltfu
# wgt: an optional numerical vector of weights (e.g. rds-ii weights)
calc.ipcw <- function(data, ltfu, prds, wgt = rep(1., nrow(data))) {
  
  stopifnot(all(ltfu == 0 | ltfu == 1))
  
  wgt <- wgt * nrow(data) / sum(wgt) # weights normalization
  
  mdl <- formula(paste('ltfu', paste(prds, collapse='+'), sep='~'))
  fit <- glm(mdl, quasibinomial, data, wgt)
  ps <- predict(fit, data, 'response')
  ps[is.na(ps)] <- mean(ps, na.rm=T)
  p <- weighted.mean(ltfu, wgt)
  ifelse(ltfu, p/ps, (1-p) / (1-ps)) * wgt
} 

# function author -- Vanessa Xiu
# compute IPCWs and perform checks
compute_ipcw <- function(cur_city, covariates, data_timept, data_list){
  data <- data_timept[data_timept$city == cur_city, ]
  # covariates_include <- find.imbalance(data, covariates, weights = "wt_rds_norm")
  data$ipw_rds <- calc.ipcw(data = data,
                            ltfu = data$ltfu,
                            prds = covariates,
                            wgt = data$wt_rds_norm)
  
  print("Summary and SMD of dataset")
  print(summary(data$ipw_rds))
  print(calc.smd(data = data[, covariates], ltfu = data$ltfu, wgt = data$ipw_rds))
  
  # check that rds x ipcw sum to rds adj. numbers of participants (both for ltfu and non-ltfu separately),
  # i.e. make sure that weighing does not inflate the sample size
  print("LTFU participants ====")
  data_ltfu <- data %>%
    subset(ltfu == 1)
  print(sum(data_ltfu$ipw_rds))
  print(sum(data_ltfu$wt_rds_norm))
  
  print("retained participants ====")
  data_fu <- data %>%
    subset(ltfu == 0)
  print(sum(data_fu$ipw_rds))
  print(sum(data_fu$wt_rds_norm))
  
  # add to list
  data_list[[cur_city]]$part_id <- data$part_id
  data_list[[cur_city]]$ipw_rds <- data$ipw_rds
  
  return(data_list)
}
