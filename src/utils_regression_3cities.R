library(Rcpp)
# sourceCpp("./code-3cities/negbin.cpp")
# Process observed data ----
#'#' compute observed distribution
compute_density_obs <- function(data, var_city = "city",wt="wt_rds_norm"){
  # weighted and unweighted frequencies
  density_observed <- data %>% 
    group_by(get(var_city), nb_part_ttl) %>% 
    summarize(nb_ppl = n(), wt_sum = sum(get(wt)), .groups = "drop_last")
  
  names(density_observed)[names(density_observed) == "get(var_city)"] <- var_city
  
  # weighted and unweighted probabilities/densities
  density_observed <- density_observed %>% 
    mutate(pop = sum(nb_ppl), dens_un = nb_ppl / pop,
           pop_wt = sum(wt_sum), dens_wt = wt_sum / pop_wt) %>% 
    ungroup()
  
  return(density_observed[, c(var_city, "nb_part_ttl", "nb_ppl", "pop", "dens_un",
                              "wt_sum", "pop_wt", "dens_wt")])
}

compute_density_obs_new <- function(data, var_city = "city",wt="wt_rds_norm"){
  # weighted and unweighted frequencies
  density_observed <- data %>% 
    group_by(get(var_city), nb_part_new) %>% 
    summarize(nb_ppl = n(), wt_sum = sum(get(wt)), .groups = "drop_last")
  
  names(density_observed)[names(density_observed) == "get(var_city)"] <- var_city
  
  # weighted and unweighted probabilities/densities
  density_observed <- density_observed %>% 
    mutate(pop = sum(nb_ppl), dens_un = nb_ppl / pop,
           pop_wt = sum(wt_sum), dens_wt = wt_sum / pop_wt) %>% 
    ungroup()
  
  return(density_observed[, c(var_city, "nb_part_new", "nb_ppl", "pop", "dens_un",
                              "wt_sum", "pop_wt", "dens_wt")])
}
# Regression result posterior checks ----
#' turn posterior predictive distributions into a data.frame
#' with a row per iteration and a column for each individual's prediction
extract_yhat_long <- function(stan_fit, row_id){
  # rows: iterations, columns: observations
  stan_pred <- extract(stan_fit, pars = "y_hat")$y_hat
  
  # assign names to each observation
  colnames(stan_pred) <- row_id
  
  # turn into a data.frame in long format
  stan_pred <- as_tibble(stan_pred)
  stan_pred <- stan_pred %>% 
    pivot_longer(cols = everything(), names_to = "part_id", values_to = "y_hat")
  
  # add iteration numbering
  stan_pred <- stan_pred %>% 
    group_by(part_id) %>% 
    mutate(iter = 1:n()) %>% 
    ungroup()
  
  # reorder and return
  return(stan_pred[, c("iter", "part_id", "y_hat")])
}

#' provide a summary of the posterior predictive distributions
summarize_yhat <- function(data, cred_int = .95){
  data %>% 
    group_by(part_id) %>% 
    summarize(mean = mean(y_hat),
              mdn = median(y_hat),
              cr.i_low = quantile(y_hat, .5 - (cred_int/2),na.rm=TRUE),
              cr.i_upp = quantile(y_hat, .5 + (cred_int/2),na.rm=TRUE),
              min = min(y_hat),
              max = max(y_hat))
}

# Process regression results ----
#' compute the PMF of a negative binomial distribution with mean mu
#' and overdispersion parameter phi, over range_x
#' returns the list of mass values over range_x
dnbinom_single <- function(mu, range_x, phi){
 dnbinom(x = range_x, size = phi, mu = mu, log=FALSE)
}

#' computes the PMF over a list of means
# and returns an array where row: values of x, and a column for each individual
# layer is each iteration
compute_dens_negbin <- function(mu_list, range_x = 0:500, phi)
  {
  pmf_array = array(NA, dim = c(length(range_x),ncol(mu_list),nrow(mu_list)))
  for (i in 1:nrow(mu_list)) #nrow(mu_list) is the # of iterations
    {
    for (j in 1:ncol(mu_list)) #ncol(mu_list) is the # of individuals
    {
      pmf_array[,j,i]<- dnbinom_single(range_x = range_x, mu = mu_list[i,j], phi = phi)
      # pmf_array[,j,i]<- dnbinom_single_cpp(range_x = range_x, mu = mu_list[i,j], phi = phi)
    }
    }
  return(pmf_array)
  }

compute_dens_negbin.single <- function(mu_list, range_x = 0:500, phi, run_par = FALSE, cl = NULL){
  if(!run_par){
    apply(cbind(mu_list), 1, dnbinom_single, range_x = range_x, phi = phi, simplify = T)
  # using this parallelization seems to actually slow down the code....
  } else {
    # start a cluster (seems to take time, better to initalize outside function)
    if(is.null(cl)) stop("Please provide an initialized cluster.")
    
    # send objects into cluster
    clusterExport(cl, varlist = c("mu_list", "dnbinom_single", "phi"), envir = environment())
    
    # perform operation in parallel
    out <- parApply(
      cl,
      X = cbind(mu_list), MARGIN = 1, FUN = dnbinom_single, # apply arguments
      range_x = range_x, phi = phi#, simplify = T # function arguments
    )
    
    # return object
    return(out)
  }
}

#' computes the expected density over all individuals in density_list
#' subsetting is possible via subset_pos, which should be a logical vector
#' indicating which individuals to keep in the computation
compute_pmf <- function(density_list, wt_values, subset_pos = NULL,
                        nb_aggregate = NULL)
  {
  
  # if no filtering is specified, use the entire sample
  if(is.null(subset_pos))
    {
    subset_pos <- rep(TRUE, length(wt_values))
    }
  
  # subset weights
  wt_values <- wt_values[subset_pos]
  
  # compute population size and sum of w_i based on arguments
  pop_ttl <- sum(subset_pos, na.rm=TRUE)
  wt_ttl <- sum(wt_values, na.rm=TRUE)
  
  # when using the aggregation method, population size is not pop_ttl but
  # the sum of how many individuals there are in each row
  # and weights need to be constructed for each row/datapoint
  if(!is.null(nb_aggregate)){
    pop_ttl <- sum(nb_aggregate)
    wt_rows <- nb_aggregate
  }
  
  ## subset data to individuals in specified relationship
  for(i in 1:length(density_list)){
    density_list[[i]] <- density_list[[i]][, subset_pos]
  }
  
  ## compute the unweighted densities over range_x **in each iteration**
  # density_unwt <- vector("list", length(density_list))
  density_all <- vector("list", length(density_list))
  
  # using foreach and dopar is slower in this instance
  for(i in 1:length(density_list)){
  # density_all <- foreach(i = 1:length(density_list), .packages = "data.table") %dopar% {
    # collapse the densities of the i-th iteration into a single data.frame, where
    # row: individual, column: value of pmf
    # with current edits it is the **opposite**
    # and subset to the chosen individuals
    # cur_dens <- as.data.frame(t(density_list[,,i])) #transpose so col = pmf, row = individual
    cur_dens <- density_list[[i]]
    
    ## compute the unweighted mean of all densities
    # by collpasing each column pmf into a simple mean
    density_unwt_tmp <- rowSums(cur_dens, na.rm=TRUE) / pop_ttl #x-range*1 
    
    ## compute the weighted mean of all densities
    # first make the matrix of weights, such that each row (individual) has
    # length(range_x) repetitions of the person's weight
    matrx_wt <- matrix(rep(wt_values, each = nrow(cur_dens)), ncol = ncol(cur_dens), byrow = F)
    
    # use a similar method for the unweighted population if using the aggregation method
    if(!is.null(nb_aggregate)){
      matrx_wt_rows <- matrix(rep(wt_rows, each = nrow(cur_dens)), ncol = ncol(cur_dens), byrow = F)
    }
    
    # verify that weights were repeated properly
    # sum(wt_values)
    # rowSums(matrx_wt)
    
    # apply weights to density
    dens_cur_iter_wted <- cur_dens * matrx_wt # element-wise multiplication
    if(!is.null(nb_aggregate)){
      dens_cur_iter_unwted <- cur_dens * matrx_wt_rows
    }
    
    # collapse each column (y value) into a weighted mean
    density_wt_tmp <- rowSums(dens_cur_iter_wted, na.rm=TRUE) / wt_ttl
    
    # X aggregation method: replace unweighted density matrix
    if(!is.null(nb_aggregate)){
      density_unwt_tmp <- rowSums(dens_cur_iter_unwted, na.rm=TRUE) / pop_ttl
    }
    
    ## format row names (y values) and columns (y_pred and corresponding density)
    density_all[[i]] <- data.table(y_pred = seq(0, nrow(cur_dens) - 1, 1),
                                   dens_un = density_unwt_tmp,
                                   dens_wt = density_wt_tmp,
                                   iter = i)
  }
  density_all <- rbindlist(density_all)
  
  # density_unwt <- density_all[type == "unwt"]
  # density_wt <- density_all[type == "wt"]
  
  # format into tibbles with a column for each density
  # tbl_density <- left_join(
  #   density_all[type == "unwt"],
  #   density_all[type == "wt"],
  #   by = c("iter", "y_pred")
  # )
  
  return(density_all)
  }


summarize_density <- function(
    data, which_dens = "dens_un", cred_int = .95
){
  # choose which density to summarize
  names(data)[names(data) == which_dens] <- "target_density"
  
  # summarize
  data[, .(mean = mean(target_density),
           mdn = median(target_density),
           cr.i_low = quantile(target_density, .5 - (cred_int/2)),
           cr.i_upp = quantile(target_density, .5 + (cred_int/2))),
       by = .(y_pred)]
}

#' computes the expected density over all individuals in density_list
#' subsetting is possible via subset_pos, which should be a logical vector
#' indicating which individuals to keep in the computation
compute_cdf <- function(density_list, wt_values, subset_pos = NULL,
                        nb_aggregate = NULL){
  
  # if no filtering is specified, use the entire sample
  # if(is.null(subset_pos)){
  #   subset_pos <- rep(TRUE, length(wt_values))
  # }
  
  # subset weights
  # wt_values <- wt_values[subset_pos]
  
  # compute population size and sum of w_i based on arguments
  # pop_ttl <- sum(subset_pos, na.rm=TRUE)
  pop_ttl <- length(wt_values)
  wt_ttl <- sum(wt_values, na.rm=TRUE)
  
  # when using the aggregation method, population size is not pop_ttl but
  # the sum of how many individuals there are in each row
  # and weights need to be constructed for each row/datapoint
  if(!is.null(nb_aggregate)){
    pop_ttl <- sum(nb_aggregate)
    wt_rows <- nb_aggregate
  }
  
  ## subset data to individuals in specified relationship
  # no longer necessary
  # for(i in 1:length(density_list)){
  #   density_list[[i]] <- density_list[[i]][, subset_pos]
  # }
  
  ## compute the unweighted densities over range_x **in each iteration**
  density_all <- vector("list", length(density_list))
  
  for(i in 1:length(density_list)){
    # subset to i-th iteration
    # row: each x-value of pmf, column: individuals
    cur_dens <- density_list[[i]]
    
    ## compute the unweighted cumulative sum of all densities by row
    # first compute CDF based on PMF
    # cur_cdf <- apply(cur_dens, 2, cumsum)
    cur_cdf <- apply(cur_dens, 2, spatstat.utils::revcumsum)
    
    # then collapse each column cdf into a simple mean
    cdf_unwt_tmp <- rowSums(cur_cdf, na.rm=TRUE) / pop_ttl #x-range*1
    
    ## compute the weighted mean of all densities
    # first make the matrix of weights, such that each row (individual) has
    # length(range_x) repetitions of the person's weight
    matrx_wt <- matrix(rep(wt_values, each = nrow(cur_dens)), ncol = ncol(cur_dens), byrow = F)
    
    # use a similar method for the unweighted population if using the aggregation method
    if(!is.null(nb_aggregate)){
      matrx_wt_rows <- matrix(rep(wt_rows, each = nrow(cur_dens)), ncol = ncol(cur_dens), byrow = F)
    }
    
    # verify that weights were repeated properly
    # sum(wt_values)
    # rowSums(matrx_wt)
    
    # apply weights to density
    cdf_cur_iter_wted <- cur_cdf * matrx_wt # element-wise multiplication
    if(!is.null(nb_aggregate)){
      cdf_cur_iter_unwted <- cur_cdf * matrx_wt_rows
    }
    
    # collapse each column (y value) into a weighted mean
    cdf_wt_tmp <- rowSums(cdf_cur_iter_wted, na.rm=TRUE) / wt_ttl
    
    # X aggregation method: replace unweighted density matrix
    if(!is.null(nb_aggregate)){
      cdf_unwt_tmp <- rowSums(cdf_cur_iter_unwted, na.rm=TRUE) / pop_ttl
    }
    
    ## format row names (y values) and columns (y_pred and corresponding cdf)
    density_all[[i]] <- data.table(y_pred = seq(0, nrow(cur_cdf) - 1, 1),
                                   cdf_un = cdf_unwt_tmp,
                                   cdf_wt = cdf_wt_tmp,
                                   iter = i)
  }
  density_all <- rbindlist(density_all)
  
  return(density_all)
}

# TODO: add confidence intervals here?
#' compute the inverse CDF, i.e. Pr(y >= k)
compute_inv_cdf <- function(data){
  # reorder
  data <- data %>% arrange(relationship)
  
  # compute CDF by summing over PMF
  data <- data %>% 
    group_by(relationship) %>%
    mutate(inv_cdf = cumsum(mean))
  
  # reorder and return
  data <- data %>% arrange(relationship)
  
  return(data)
}
