compute_grp_size <- function(pmf_estim_m, order_grps = FALSE){
  # create data.frame with a y_pred and a density column (PMF at Y=k)
  data_pmf_summ <- data.frame(y_pred = 0:(length(pmf_estim_m) - 1),
                              pmf_m = pmf_estim_m)
  
  # divide into 100 sexual activity groups
  for (j in 1:99) {
    i = 99 - j + 2
    data_pmf_summ <- data_pmf_summ %>%
      mutate(quant_grp = ifelse(y_pred <= 3 * i, paste0("grp_", i), quant_grp))
  }
  
  data_pmf_summ <- data_pmf_summ %>%
    mutate(quant_grp = case_when(
      y_pred == 0         ~ "grp_1",
      TRUE ~ quant_grp))
  
  # order groups for data inspection
  if(order_grps){
    data_pmf_summ$quant_grp <- factor(data_pmf_summ$quant_grp,
                                      levels = paste0("grp_", 1:100))
  }
  
  # add group density
  data_pmf_summ <- data_pmf_summ %>% 
    # group_by(data_pt, quant_grp) %>%
    group_by(quant_grp) %>% 
    summarize(y_mean = mean(y_pred), 
              qt = sum(pmf_m),
              .groups = "drop")
  
  return(data_pmf_summ)
}

compute_r0_ngm <- function(pmf_estim_m, beta_range, run_par = FALSE){
  if(run_par){
    library(dplyr)
  }
  
  # R0 estimation code ----
  ## Compute group densities ----
  # get contact rates and group density (dij = Nij/Nj where i is a sexual activity group, j is a city-timept)
  # pmf_estim_m <- data_pmf_city_timept[cur_iter, ] # delete??
  data_pmf_summ <- compute_grp_size(pmf_estim_m, T)
  
  n = length(unique(data_pmf_summ$quant_grp)) # number of sexual activity groups
  
  # compute the contact rates
  # c: number of sexual contacts per day (estimate from P6M data)
  # d: proportion of the population in that sexual activity group
  contact <- data_pmf_summ %>% 
    transmute(c = y_mean / (6*30), d = qt)
  
  ## Disease parameters ----
  # incubation period
  v <- 1/(7.9)
  
  # disease duration
  gamma <- 1/(17.3)
  
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
  # for(pt in data_pt_name){
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
    
    #' the columns are in the order of E_i -> I_i -> E_i+1 ....
    #' where i is the sexual activity group
    for(beta in beta_range){
      for(i in 1:n){ 
        for(j in 1:n){
          T_pt[2 * i - 1, 2 * j] <- c_pt[i] * d_pt[i] * c_pt[j] * beta / denom_pt
        }
      }
      K_pt = -t(E) %*% T_pt %*% solve(Sigma) %*% E
      
      # verified that using NGM with large domain K_l = -T *sigma^-1
      # gives same result as below
      # K_alt <- -T_pt %*% solve(Sigma)
      # Re(eigen(K_alt)$values[1])
      
      R0_pt <- Re(eigen(K_pt)$values[1])
      
      R0[match(beta, beta_range)] <- R0_pt
    }
  }
  
  return(R0)
}

