compute_grp_size <- function(pmf_estim_m, order_grps = FALSE,
                             finer_categorization = FALSE){
  # create data.frame with a y_pred and a density column (PMF at Y=k)
  data_pmf_summ <- data.frame(y_pred = 0:(length(pmf_estim_m) - 1),
                              pmf_m = pmf_estim_m)
  if(!finer_categorization){
  # divide into 100 sexual activity groups
  for (j in 1:99) {
    i = 99 - j + 2
    # print(c("j", j, "i", i, 3 * i))
    data_pmf_summ <- data_pmf_summ %>% 
      mutate(quant_grp = ifelse(y_pred <= 3 * i, paste0("grp_", i), quant_grp))
  }

  data_pmf_summ <- data_pmf_summ %>% 
    mutate(quant_grp = case_when(
      y_pred == 0         ~ "grp_1",
      TRUE ~ quant_grp))
    # data_pmf_summ <- data_pmf_summ %>% 
    #   mutate(quant_grp = case_when(y_pred == 0 ~ "grp_1",
    #                                y_pred <= 5 ~ "grp_2",
    #                                y_pred <= 10 ~ "grp_3",
    #                                y_pred <= 25 ~ "grp_4",
    #                                y_pred >  25 ~ "grp_5"))
  # group_num = 5
  group_num = 100
  
  }
  else{
  #  248 sexual activity groups, finer group (2 sexual partner number per group) for >= 25 partners
  # (1 sexual partner number per group) for >= 100 partners
  for (i in 300:100) {
    index = i - 52
    data_pmf_summ <- data_pmf_summ %>%
      mutate(quant_grp = ifelse(y_pred <=  i, paste0("grp_", index), quant_grp))
  }
  for (i in 50:13) {
    index = i - 3
    data_pmf_summ <- data_pmf_summ %>%
      mutate(quant_grp = ifelse(y_pred <  2 * i, paste0("grp_", index), quant_grp))
  }
  for (i in 8:0) {
    index = i + 1
    data_pmf_summ <- data_pmf_summ %>%
      mutate(quant_grp = ifelse(y_pred <= 3 * i, paste0("grp_", index), quant_grp))
  }
    group_num = 248
  }
  
  # order groups for data inspection
  if(order_grps){
    data_pmf_summ$quant_grp <- factor(data_pmf_summ$quant_grp,
                                      levels = paste0("grp_", 1:group_num))
  }
  
  # add group density
  data_pmf_summ <- data_pmf_summ %>% 
    # group_by(data_pt, quant_grp) %>%
    group_by(quant_grp) %>% 
    summarize(y_mean = mean(y_pred),
    # summarize(y_mean = sum(y_pred * pmf_m) / sum(pmf_m),
              qt = sum(pmf_m),
              min = min(y_pred), max = max(y_pred),
              .groups = "drop")
  
  return(data_pmf_summ)
}

#' @param pfm_estim_m    PMF of the distribution of sexual partners in the P6M
#' @param beta_range     range of per-partnership infection risks to use
#' @param epsilon        assortativity coefficient
compute_r0_ngm <- function(pmf_estim_m, beta_range, epsilon = 0.8,
                           run_par = FALSE,
                           cap = -99,
                           # k = k){
                           finer = TRUE){
  if(run_par){
    library(dplyr)
  }
  
  # R0 estimation code ----
  ## Compute group densities ----
  # get contact rates and group density (dij = Nij/Nj where i is a sexual activity group, j is a city-timept)
  # pmf_estim_m <- data_pmf_city_timept[cur_iter, ] # delete??
  # data_pmf_summ <- compute_grp_size(pmf_estim_m, T, k = k)
  if(cap == -99){
    data_pmf_summ <- compute_grp_size(pmf_estim_m, T, finer_categorization = FALSE)
  } else {
    data_pmf_summ <- compute_grp_size(pmf_estim_m, T, cap = cap)
  }
  
  # data_pmf_summ$y_mean <- c(rep(0, 99), 299)
  # data_pmf_summ$y_mean <- c(rep(50, 100))
  # data_pmf_summ$qt <- rev(data_pmf_summ$qt)
  # data_pmf_summ$qt <- data_pmf_summ$qt * 1000
  #' adding an **empty group** should NOT change results whatsoever
  #' however it DOES change the results when using assortativity
  # data_pmf_summ <- data_pmf_summ %>% add_row(y_mean = 50000, qt = min(data_pmf_summ$qt) / 10000000000)
  
  # removes 0-partnership group from NGM
  # data_pmf_summ <- data_pmf_summ %>% filter(y_mean > 0)
  
  n = length(unique(data_pmf_summ$quant_grp)) # number of sexual activity groups
  # browser()
  # compute the contact rates
  # c: number of sexual contacts per day (estimate from P6M data and adjusted to daily timescale)
  # d: proportion of the population in that sexual activity group
  contact <- data_pmf_summ %>% 
    transmute(c = y_mean / (6*30), d = qt)
  
  ## Disease parameters ----
  si = 9.0               # estimate from literature
  incub = 7.1            # estimate from literature
  latent = incub - 2     # assume 2 days of pre-symptomatic transmission
  duration = si - latent # latent and infectious periods add up to serial interval
  
  # incubation period
  v <- 1 / latent
  
  # disease duration
  gamma <- 1 / duration
  # gamma <- 1 / 17
  
  ## Constituents of NGM ----
  R0 <- rep(-1, length(beta_range))
  
  # set-up
  m = 2 * n # need E and I compartments for each group
  Sigma = matrix(0, m, m)
  E = matrix(0, m, n)
  
  # Transition matrix sigma
  for(i in 1:n){
    Sigma[2*i - 1, 2*i - 1] <- -v
    Sigma[2*i, 2*i - 1] <- v
    Sigma[2*i, 2*i] <- -gamma
  }
  
  # Auxiliary matrix E
  for(i in 1:n){
    E[2*i - 1, i] <- 1
  }
  
  ## Computation of NGM & R0 ----
  {
    # Transmission matrix T
    T_pt = matrix(0, m, m)
    
    # c: city and time-specific sexual activity rate
    # d: proportion in sexual activity group (instead of using population size N)
    c_pt = contact$c
    d_pt = contact$d
    
    # c_ttl:    total sexual activity rate
    #           c_1 + c_2 + c_3 + ... + c_n
    # denom_pt: total number of partnerships generated in the population
    #           c_1*d_1 + c_2*d_2 + ... c_n*d_n
    c_ttl = sum(c_pt)
    
    denom_pt = sum(c_pt * d_pt)
    
    ## create assortative mixing matrix
    # TODO coarsen the activity groups
    # assortative
    # mix_assort <- diag(n)
    
    # sexual activity grouped as 0, 1-6, 7-9, 10-12, ....
    # TODO ADAPT FOR FINER SCALE
    mix_assort = matrix(0, n, n)
    
    # the upper bound of the index of each group
    # super_grps <- c(1, 2, 5, 10, 17, 33, 100)
    # super_grps <- c(1, 2, 5, 10, 17, n)
    # super_grps <- c(1, 2, 3, n)
    super_grps <- c(1, 2, 5, 8, n)
    
    for(indx in 1:length(super_grps)){
      i <- super_grps[indx]
      j <- super_grps[indx - 1]
      if(i == 1 | i == 2){
        # fully assortative for first two groups
        mix_assort[i, i] = 1
      } else {
        # assign a square matrix from the start of the last group
        # up to the upper boundary of the current group
        mix_assort[(j + 1):i, (j + 1):i] = 1 / (i - j)
      }
    }
    rm(super_grps, i, j)
    # mix_assort[i, j] = case_when(i == 1 | i == 2 ~ 1,
    #                              i <= 5          ~ 1/(5   - 2),
    #                              i <= 10         ~ 1/(10  - 5),
    #                              i <= 17         ~ 1/(17  - 10),
    #                              i <= 33         ~ 1/(33  - 17),
    #                              i <= 67         ~ 1/(67  - 33),
    #                              i <= 100        ~ 1/(100 - 67))
    
    # check that the assortative mixing matrix has been properly constructed
    if(epsilon == 0){
      range_mix_assort <- range(rowSums(mix_assort))
      if(range_mix_assort[1] != 1 | range_mix_assort[2] != 1){
        stop("Verify assortative mixing matrix")
      }
    }
    # rowSums(mix_assort)
    
    # FOR SIMPLER R0
    R0_simple <- matrix(0, n, n)
    
    #' the columns are in the order of E_i -> I_i -> E_i+1 ....
    #' where i is the sexual activity group
    for(beta in beta_range){
      for(i in 1:n){
        for(j in 1:n){
          # transmission matrix using fully proportionate mixing
          # T_pt[2 * i - 1, 2 * j] <- c_pt[i] * d_pt[i] * c_pt[j] * beta / denom_pt

          # transmission matrix using fully proportionate mixing
          T_pt[2 * i - 1, 2 * j] <- beta * c_pt[i] * d_pt[i] * (1 / d_pt[j]) *
            (
              ( (c_pt[j] * d_pt[j] / denom_pt) * (1 - epsilon) ) + (mix_assort[i, j] * epsilon)
            )
          
          # R0_simple <- beta * c_pt[i] * duration * g_ij
          # R0_simple[i, j] <- beta * c_pt[i] * duration * (
          #   (c_pt[j] * d_pt[j] / denom_pt) * (1 - epsilon) + (mix_assort[i, j] * epsilon)
          # )
        }
      }
      
      # compute NGM using -E' * T * Sigma^-1 * E
      K_pt = -t(E) %*% T_pt %*% solve(Sigma) %*% E
      
      # verified that using NGM with large domain K_l = -T *sigma^-1
      # gives same result as below
      # K_alt <- -T_pt %*% solve(Sigma)
      # Re(eigen(K_alt)$values[1])
      
      R0_pt <- Re(eigen(K_pt)$values[1])
      # R0_pt <- Re(eigen(R0_simple)$values[1])
      
      R0[match(beta, beta_range)] <- R0_pt
      # R0 <- R0_pt
    }
  }
  
  return(R0)
}


