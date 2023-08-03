# Descriptive ----
# basic functions to plot density colouring by specified characteristic
plot_density <- function(
    data, fill_by,
    log_scale = FALSE, alpha = 0.2
  ){
  p <- ggplot(data, aes_string(x = "nb_part_ttl", fill = fill_by)) +
    geom_density(alpha = alpha) +
     
    labs(x = "Number of partners in P6M") +
    theme(legend.position = "bottom")
  
  if(log_scale){
    p <- p +
      scale_x_continuous(trans = scales::pseudo_log_trans(),
                         breaks = c(0, 10, 20, 40, 100, 200, 300))
  }
  p
}

plot_violin <- function(
    data, x_name = "nb_part_ttl", strat_by = NULL,
    log_scale = FALSE, show_pts = TRUE, show_central_tendency = TRUE,
    alpha = 0.2
  ){
  p <- ggplot(data, aes_string(x = x_name, y = strat_by)) +
    
    labs(x = "Number of partners in P6M") +
    theme(legend.position = "none")
  
  # order y axis correctly
  if(strat_by != 1){
    # get variable name
    y_name <- gsub("as\\.character|\\(|\\)", "", strat_by)
    
    y_levels <- unique(data[[y_name]])
    y_levels <- sort(y_levels, na.last = T)
    y_levels <- as.character(y_levels)
    y_levels[y_levels == "NA"] <- NA
    
    p <- p +
      scale_y_discrete(limits = rev(y_levels))
  }
  
  # additional layers/options
  if(log_scale){
    p <- p +
      scale_x_continuous(trans = scales::pseudo_log_trans(),
                         breaks = c(0, 10, 20, 40, 100, 200, 300))
  }
  if(show_pts){
    p <- p +
      geom_jitter(alpha = 0.05, size = 0.2, height = 0.20)
  }
  # add main layer (done after to superimpose over individual data points)
  p <- p +
    geom_violin(aes_string(fill = strat_by), alpha = 0.5)
  # add central tendency over violin distribution
  if(show_central_tendency){
    if(strat_by != 1){
      data_central <- data %>% 
        group_by(across(all_of(y_name))) %>% 
        summarize(mean = mean(get(x_name)), median = median(get(x_name)),
                  quant.90 = quantile(get(x_name), .90))
    } else {
      data_central <- data %>% 
        summarize(mean = mean(get(x_name)), median = median(get(x_name)),
                  quant.90 = quantile(get(x_name), .90))
    }
    
    p <- p +
      # mean
      geom_point(data = data_central, aes_string(x = "mean", y = strat_by),
                 shape = 4, size = 1) +
      # median
      geom_point(data = data_central, aes_string(x = "median", y = strat_by),
                 shape = 7, size = 1) +
      # 90th quantile
      geom_point(data = data_central, aes_string(x = "quant.90", y = strat_by),
                 shape = 17, size = 1)
  }
  
  p
}

# Data verification ----
plot_hist_indiv <- function(data_pred, data_obs, id_row){
  # set layout
  par(mfrow = c(2, 1))
  
  # regular scale
  hist(data_pred[, id_row], breaks = 50); abline(v = data_obs$nb_part_ttl[id_row], col = "red")
  # log scale
  hist(log(data_pred[, id_row]), breaks = 50); abline(v = log(data_obs$nb_part_ttl[id_row]), col = "red")
}

# Outputs -----
plot_dens_fit_obs <- function(
    density_fit, density_observed, name_density = "dens_un",
    xmax = 150, log_scale = FALSE,
    title = NULL
    ){
  ## expand so that every number is represented in the data frame
  density_observed <- left_join(
    expand.grid(nb_part_ttl = min(density_observed$nb_part_ttl):max(density_observed$nb_part_ttl)),
    density_observed, by = "nb_part_ttl"
  )
  density_observed <- density_observed %>% 
    mutate(dens_un = ifelse(is.na(dens_un), 0, dens_un),
           dens_wt = ifelse(is.na(dens_wt), 0, dens_wt))
  
  ## plot
  # base layer with fit from regression
  p <- ggplot(density_fit, aes(x = y_pred)) +
    # empirical distribution
    geom_col(data = density_observed, aes(x = nb_part_ttl, y = get(name_density)),
             fill = "#c3f4f6", alpha = 0.9) +
    geom_line(data = density_observed, aes(x = nb_part_ttl, y = get(name_density), col = "Observed"),
              alpha = 0) +
    # fitted distribution
    geom_line(aes(y = mean, col = "Fitted")) +
    geom_ribbon(aes(ymin = cr.i_low, ymax = cr.i_upp), alpha = 0.2) +
    
    # aesthetic settings
    scale_x_continuous(breaks = 20 * 0:8) +
    coord_cartesian(xlim = c(0, xmax), ylim = c(0, .20)) +
    # colours from met.brewer("Homer1")[1:8]
    scale_colour_manual(values = c("Fitted" = "#a62f00", "Observed" = "#32b2da")) +
    
    labs(title = title,
         x = "Number of partners in P6M", y = "Density", col = NULL) +
    theme(panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank())
  
  if(log_scale){
    p <- p + 
      scale_x_continuous(trans = scales::pseudo_log_trans(),
                         breaks = c(0, 1, 5, 25, 50, 100, 150))
  }
  
  return(p)
}


plot_dens_fit_obs_new <- function(
    density_fit, density_observed, name_density = "dens_un",
    xmax = 150, log_scale = FALSE,
    title = NULL
){
  ## expand so that every number is represented in the data frame
  density_observed <- left_join(
    expand.grid(nb_part_new = min(density_observed$nb_part_new):max(density_observed$nb_part_new)),
    density_observed, by = "nb_part_new"
  )
  density_observed <- density_observed %>% 
    mutate(dens_un = ifelse(is.na(dens_un), 0, dens_un),
           dens_wt = ifelse(is.na(dens_wt), 0, dens_wt))
  
  ## plot
  # base layer with fit from regression
  p <- ggplot(density_fit, aes(x = y_pred)) +
    # empirical distribution
    geom_col(data = density_observed, aes(x = nb_part_new, y = get(name_density)),
             fill = "#c3f4f6", alpha = 0.9) +
    geom_line(data = density_observed, aes(x = nb_part_new, y = get(name_density), col = "Observed"),
              alpha = 0) +
    # fitted distribution
    geom_line(aes(y = mean, col = "Fitted")) +
    geom_ribbon(aes(ymin = cr.i_low, ymax = cr.i_upp), alpha = 0.2) +
    
    # aesthetic settings
    scale_x_continuous(breaks = 20 * 0:8) +
    coord_cartesian(xlim = c(0, xmax), ylim = c(0, .20)) +
    # colours from met.brewer("Homer1")[1:8]
    scale_colour_manual(values = c("Fitted" = "#a62f00", "Observed" = "#32b2da")) +
    
    labs(title = title,
         x = "Number of partners in P6M", y = "Density", col = NULL) +
    theme(panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank())
  
  if(log_scale){
    p <- p + 
      scale_x_continuous(trans = scales::pseudo_log_trans(),
                         breaks = c(0, 1, 5, 25, 50, 100, 150))
  }
  
  return(p)
}







plot_dens_fit_compare <- function(
    density_fits, col_var,
    xmax = 150, log_scale = FALSE, cred_int = TRUE,
    title = NULL
    ){
  ## plot
  p <- ggplot(density_fits, aes(y_pred)) +
    geom_line(aes(y = mean, col = get(col_var))) +
    
    scale_x_continuous(breaks = 20 * 0:8) +
    coord_cartesian(xlim = c(0, xmax)) +
    labs(x = "Number of partners in P6M", y = "Density", col = "", fill = "") +
    theme(panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank())
  
  ## optional plot layers
  if(cred_int){
    p <- p +
      geom_ribbon(aes(ymin = cr.i_low, ymax = cr.i_upp, fill = get(col_var)), alpha = 0.2)
  }
  
  if(log_scale){
    p <- p + 
      scale_x_continuous(trans = scales::pseudo_log_trans(),
                         breaks = c(0, 1, 5, 25, 50, 100, 150))
  }
  
  return(p)
}

plot_cdf_fit <- function(
    data, col_var = NULL, cred_int = TRUE,
    p_title = NULL,
    outcome_type,
    linewidth = 0.25){
  ## plot
  p <- ggplot(data, aes(x = y_pred, y = mean)) +
    labs(title = p_title,
         x = "Number of reported sexual partners in the P6M",
         y = "Proportion reporting at least x partners")
  
  # credible interval
  if(cred_int){
    p <- p +
      geom_ribbon(aes(ymin = cr.i_low, ymax = cr.i_upp, fill = get(col_var)), alpha = 0.6)
  }
  
  # base layer
  p <- p +
    geom_line(aes(col = get(col_var)), linewidth = linewidth)
  
  # colour setting for 3-city comparison
  if(col_var %in% c("city", "city_name", "interval")){
    p <- p +
      scale_fill_viridis_d(option = "C", end = 0.8) +
      scale_colour_viridis_d(option = "C", end = 0.8)
  }
  
  return(p)
}

# plot select CDF values as point estimates
plot_cdf_sel <- function(
    data, col_var = NULL, file_suff = NULL,
    p_title = NULL, size_pt = .8, size_line = 0.4){
  ## plot
  p <- ggplot(data, aes(x = factor(y_pred), col = get(col_var))) +
    # point estimate
    geom_point(
      aes(y = mean),
      size = size_pt,
      position = position_dodge(0.4)
    ) +
    # credible interval
    geom_linerange(
      aes(ymin = cr.i_low, ymax = cr.i_upp),
      linewidth = size_line,
      position = position_dodge(0.4)
    ) +
    labs(title = p_title,
         # x = sprintf("Number of reported sexual partners in P6M", gsub("p6m_", "", file_suff)),
         x = "Number of reported sexual partners in P6M",
         y = "Proportion reporting at least x partners",
         col = col_var)
  
  # colour setting for 3-city comparison
  if(col_var %in% c("city", "city_name", "interval")){
    p <- p +
      scale_fill_viridis_d(option = "C", end = 0.8) +
      scale_colour_viridis_d(option = "C", end = 0.8)
  }
  
  return(p)
}

plot_inverse_cdf <- function(
    cdf_fit, col_var,
    xmax = 100, log_scale = FALSE, cred_int = TRUE,
    title = NULL
    ){
  ## plot
  # base layer with fit from regression
  p <- ggplot(cdf_fit, aes(x = y_pred)) +
    # fitted distribution
    geom_line(aes(y = inv_cdf, col = get(col_var))) +
    
    # aesthetic settings
    scale_x_continuous(breaks = 20 * 0:8) +
    coord_cartesian(xlim = c(0, xmax), ylim = c(0, 1)) +
    
    labs(title = title,
         x = "Number of partners in P6M",
         y = "Pr(having x or more partners)", col = NULL, fill = NULL) +
    theme(panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank())
  
  ## optional plot layers
  if(cred_int){
    p <- p +
      geom_ribbon(aes(ymin = cr.i_low, ymax = cr.i_upp, fill = get(col_var)), alpha = 0.2)
  }
  
  if(log_scale){
    p <- p + 
      scale_x_continuous(trans = scales::pseudo_log_trans(),
                         breaks = c(0, 1, 5, 25, 50, 100, 150))
  }
  
  return(p)
}
