# Descriptive ----
# incidence from PHAC data
plot_phac_case <- function(data, yvar = "incidence", facet = FALSE){
  p <- ggplot(mutate(data, prov = factor(prov, levels = c("QC", "ON", "BC"))),
              aes(x = date, y = get(yvar), col = prov)) +
    
    coord_cartesian(ylim = c(0, max(data[, yvar]))) +
    
    scale_colour_viridis_d(option = "C", end = 0.8)
  
  # separate into facets
  if(facet){
    p <- p +
      facet_wrap(~prov) +
      theme(legend.position = "none")
  }
  
  return(p)
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

### add the verification plot

# Outputs -----
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
