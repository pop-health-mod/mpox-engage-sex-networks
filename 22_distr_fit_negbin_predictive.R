# source code ---------------------------
source("./21_distr_fit_negbin_coef.R")
# get posterior predictive distribution ----
# get distribution of expected predictions
# rows: iterations, columns: each individual's prediction (mean of the NB)
data_ypred <- create_city_list(CITIES_DATAPTS)
phi_all <- create_city_list(CITIES_DATAPTS)

for(cur_city in CITIES_DATAPTS){
  data_ypred[[cur_city]] <- extract(fit_bayes_ls[[cur_city]], pars = "y_pred")$y_pred
  phi_all[[cur_city]] <- extract(fit_bayes_ls[[cur_city]], pars = "phi")$phi
}

density_list_by_iter <- create_city_list(CITIES_DATAPTS)

# start a cluster
numCores <- parallel::detectCores()
numCores
numCores <- numCores - 1
cl <- parallel::makeCluster(numCores, type = "PSOCK")

# register
doParallel::registerDoParallel(cl = cl) 
# foreach::getDoParRegistered()

t0 <- Sys.time()
for(cur_city in CITIES_DATAPTS){
  
  print(cur_city)
  
  density_list_by_iter[[cur_city]] <- foreach(i = 1:nrow(data_ypred[[cur_city]])) %dopar% {
    compute_dens_negbin.single(
      data_ypred[[cur_city]][i, ],
      range_x = 0:300,
      phi = phi_all[[cur_city]][i]
    )
  }
}
t1 <- Sys.time()
t1 - t0 

# close cluster
stopCluster(cl)

#   list                          list      matrix
# density_list_by_iter[[city]][[iteration]][0:300, individuals]

# make lists to hold data for each city
tmp_density <- create_city_list(CITIES_DATAPTS) # 4000 iterations per city-data_pt 
dens_wt_by_city <- create_city_list(CITIES_DATAPTS)

t0 <- Sys.time()

for(cur_city in CITIES_DATAPTS){
  cat(sprintf("%s\n", cur_city))
  # compute density in the city and visit
  tmp_density[[cur_city]] <- compute_pmf(
    density_list_by_iter[[cur_city]],
    filter(data_x_aggrt, data_pt == cur_city)$ipw_rds,
    nb_aggregate = filter(data_x_aggrt, data_pt == cur_city)$nb
  )
  
  # # all densities should sum up to 1 (VERIFIED)
  # tmp_density[[cur_city]] %>% group_by(iter) %>% summarize(dens_ttl = sum(dens_wt)) %>%
  #   pull(dens_ttl) %>% unique() %>% print()
  
  dens_wt_by_city[[cur_city]] <- summarize_density(tmp_density[[cur_city]], "dens_wt")
}

t1 <- Sys.time()
t1 - t0 


# tail comparison ---------------------------------------------------------

degree <- 100  # comparing density of participants reporting at least x degree of partners

# across timepoint comparison

prop_by_city <- create_city_list(CITIES)

for(city_name in CITIES){

  prop_by_city[[city_name]] <- data.frame(pre_pand = 0, post_pand = 0, pre_post = 0)

  for(iter_num in 1:4000){
    cdf_city_pre <- sum(filter(tmp_density[[paste0(city_name, "-Pre-Pandemic")]],
                               iter == iter_num & y_pred >= degree)$dens_wt)
    cdf_city_pand <- sum(filter(tmp_density[[paste0(city_name, "-Pandemic")]],
                                iter == iter_num & y_pred >= degree)$dens_wt)
    cdf_city_post <- sum(filter(tmp_density[[paste0(city_name, "-Post-Restriction")]],
                                iter == iter_num & y_pred >= degree)$dens_wt)

    prop_by_city[[city_name]]$pre_pand <- prop_by_city[[city_name]]$pre_pand + as.numeric(cdf_city_pre >= cdf_city_pand)
    prop_by_city[[city_name]]$post_pand <- prop_by_city[[city_name]]$post_pand + as.numeric(cdf_city_post >= cdf_city_pand)
    prop_by_city[[city_name]]$pre_post <- prop_by_city[[city_name]]$pre_post + as.numeric(cdf_city_pre >= cdf_city_post)
  }

  prop_by_city[[city_name]]$pre_pand <- prop_by_city[[city_name]]$pre_pand/4000
  prop_by_city[[city_name]]$post_pand <- prop_by_city[[city_name]]$post_pand/4000
  prop_by_city[[city_name]]$pre_post <- prop_by_city[[city_name]]$pre_post/4000
}

# across city comparison

prop_by_timept <- create_city_list(TIMEPTS)

for(time_pt_name in TIMEPTS){

  prop_by_timept[[time_pt_name]] <- data.frame(trt_mtl = 0, trt_van = 0, van_mtl = 0)

  for(iter_num in 1:4000){
    cdf_mtl_timept <- sum(filter(tmp_density[[paste0("Montreal-", time_pt_name)]],
                                 iter == iter_num & y_pred >= degree)$dens_wt)
    cdf_trt_timept <- sum(filter(tmp_density[[paste0("Toronto-", time_pt_name)]],
                                 iter == iter_num & y_pred >= degree)$dens_wt)
    cdf_van_timept <- sum(filter(tmp_density[[paste0("Vancouver-", time_pt_name)]],
                                 iter == iter_num & y_pred >= degree)$dens_wt)

    prop_by_timept[[time_pt_name]]$trt_mtl <- prop_by_timept[[time_pt_name]]$trt_mtl + as.numeric(cdf_trt_timept >= cdf_mtl_timept)
    prop_by_timept[[time_pt_name]]$trt_van <- prop_by_timept[[time_pt_name]]$trt_van + as.numeric(cdf_trt_timept >= cdf_van_timept)
    prop_by_timept[[time_pt_name]]$van_mtl <- prop_by_timept[[time_pt_name]]$van_mtl + + as.numeric(cdf_van_timept >= cdf_mtl_timept)
  }
  prop_by_timept[[time_pt_name]]$trt_mtl <- prop_by_timept[[time_pt_name]]$trt_mtl/4000
  prop_by_timept[[time_pt_name]]$trt_van <- prop_by_timept[[time_pt_name]]$trt_van/4000
  prop_by_timept[[time_pt_name]]$van_mtl <- prop_by_timept[[time_pt_name]]$van_mtl/4000
}

# collapse into single dataset --------------------------------------------
# collapse each city

for(cur_city in CITIES_DATAPTS){
  dens_wt_by_city[[cur_city]]$data_pt <- cur_city
}


dens_wt_by_city <- bind_rows(dens_wt_by_city)

# reorder to put city and relationship first
dens_wt_by_city <- dens_wt_by_city[, .(data_pt, y_pred, mean, mdn, cr.i_low, cr.i_upp)]

# verification of posterior distributions and densities ----
# verify results by looking at the mean number of partners
data_mean_nb_partn <- dens_wt_by_city %>%
    group_by(data_pt) %>%
    summarize(mean_wt = sum(y_pred * mean),
              cr.i_low = sum(y_pred * cr.i_low),
              cr.i_upp = sum(y_pred * cr.i_upp),
              .groups = "drop") %>% 
    mutate(type = "neg. bin.")

data_mean_nb_partn$data_pt <- factor(data_mean_nb_partn$data_pt, levels = CITIES_DATAPTS)

# verify that weights sum up to 1
dens_wt_by_city %>% 
  group_by(data_pt) %>% 
  summarize(dens_ttl = sum(mean), .groups = "drop")

# output tables (pmf and cdf) ----
# full fitted distribution pmf
write.csv(dens_wt_by_city, sprintf("./output-3cities/negbin-ipcw/full_fitted_pmf_wt_%s.csv", file_suff),
          row.names = F)

# compute fitted CDF in each datapoint (similar procedure as for pmf)
# make lists to hold data for each city
cdf_wt_by_city <- create_city_list(CITIES_DATAPTS)

t0 <- Sys.time()
for(cur_city in CITIES_DATAPTS){
  cat(sprintf("%s\n", cur_city))
  # compute CDF in the city/timepoint
  tmp_density <- compute_cdf(
    density_list_by_iter[[cur_city]],
    filter(data_x_aggrt, data_pt == cur_city)$ipw_rds,
    nb_aggregate = filter(data_x_aggrt, data_pt == cur_city)$nb
  )
  
  # compute weighted densities
  cdf_wt_by_city[[cur_city]] <- summarize_density(tmp_density, "cdf_wt")
}
t1 <- Sys.time()
t1 - t0 

# collapse into single dataset
for(cur_city in CITIES_DATAPTS){
  cdf_wt_by_city[[cur_city]]$data_pt <- cur_city
}

# collapse into a single dataframe
cdf_wt_by_city <- bind_rows(cdf_wt_by_city)

# reorder to put city and relationship first
cdf_wt_by_city <- cdf_wt_by_city[, .(data_pt, y_pred, mean, mdn, cr.i_low, cr.i_upp)]


cdf_wt_by_city <- cdf_wt_by_city %>% 
  mutate(city = gsub("-Pre-Pandemic|-Pandemic|-Post-Restriction", "", data_pt),
         time_pt = gsub("Montreal-|Toronto-|Vancouver-", "", data_pt),
         .before = 1)

# factorize time points to order correctly
cdf_wt_by_city$time_pt <- factor(cdf_wt_by_city$time_pt, 
                                 levels = c("Pre-Pandemic", "Pandemic", "Post-Restriction"))

# save full fitted distribution
write.csv(cdf_wt_by_city, sprintf("./output-3cities/negbin-ipcw/full_fitted_cdf_wt_%s.csv", file_suff),
          row.names = F)