# Library and data ----
library(tidyverse)

data_pmf <- read_csv("./misc-data-proc/outputs-cm/engage_van_fitted_pmf_nb_partn_p6m.csv")

## create the sexual activity groups
# see Xiu et al. (2024) for the detailed definitions of the time-periods
# https://doi.org/10.1093/infdis/jiae033
#    Prepandemic             = Feb 2017 to Aug 2019 (baseline visit)
#    (COVID) pandemic period = Jun 2020 to Nov 2021
#    Postrestrictions period = Dec 2021 to Feb 2023

# can choose here which time period(s) to use
TIMEPTS <- c("Pre-Pandemic", "Pandemic", "Post-Restrictions")

df_results_ls <- list(pre_pand = NULL, pand = NULL, post_restr = NULL)

for(cur_time_pt in TIMEPTS){
  print(paste(cur_time_pt, "========================"))
  
  # dataframe to store results
  data_pmf_tmp <- filter(data_pmf, time_pt == cur_time_pt)
  
  {
    # Categorize sex. actv. groups based on pre-specified percentiles ----
    ## Setup ----
    # percentiles at which to create the sex actv. groups
    # these can be modified to make larger / smaller groups, and
    # also to increase / decrease the number of groups
    cdf_cutoff <- c(.20, .30, .40, .50,  .60,  .70,  .80,  .90,  .92, .94,
                    .95, .96, .97, .98, .99, .992, .994, .996, .998, 1)
    
    # fitted distributions of sexual partnerships
    pmf_size <- data_pmf_tmp$mean
    pmf_rate <- data_pmf_tmp$y_pred
    
    # sexual activity groups to create
    k_size <- vector("double", 0)
    k_rate <- vector("double", 0)
    
    # which points of the PMF are included
    k_min <- vector("double", 0)
    k_max <- vector("double", 0)
    
    # correct the PMF by reweighing to 1
    if( sum(pmf_size) < 1 ){
      print(sprintf("pmf was reweighted, previously added up to %s", sum(pmf_size)))
      pmf_size <- pmf_size / sum(pmf_size)
    }
    
    ## Categorization loop ----
    for(i in 1:length(cdf_cutoff)){
      # print(i)
      # verify before running loop that the PMF adds up to 1
      if( i == 1 & round(sum(pmf_size), 10) < 1 ) { stop("error; sum(PMF) < 1") }
      
      # take x% from the current group
      x_perc <- ifelse(i == 1,
                       cdf_cutoff[i],
                       cdf_cutoff[i] - cdf_cutoff[i - 1])
      
      if(pmf_size[1] >= x_perc){
        ## if we have enough people we take the required density, assign rate and go to next group
        k_size[i] <- min(x_perc, pmf_size[1])
        pmf_size[1] <- pmf_size[1] - k_size[i]
        
        k_rate[i] <- pmf_rate[1]
        
        # store min and max
        k_min[i] <- k_max[i] <- pmf_rate[1]
      } else {
        ## otherwise pull from additional groups
        ## and get a rate that is the weighted average of the partnerships in all groups
        vec_rate <- vector("double", 0)
        vec_size <- vector("double", 0)
        
        # tracks how much we have left to assign
        x_perc_left <- x_perc
        j <- 1
        
        while(round(sum(vec_size), 10) < round(x_perc, 10)){
          # first note the rate and the contribution of the group
          vec_size[j] <- min(x_perc_left, pmf_size[1])
          pmf_size[1] <- pmf_size[1] - vec_size[j]
          
          vec_rate[j] <- pmf_rate[1]
          
          # check progress
          # if(pmf_rate[1] %% 50 == 0) { print(pmf_rate[1]) }
          
          # then remove PMF point from vector
          if(pmf_size[1] <= 0 & length(pmf_size) > 1){
            pmf_size <- pmf_size[2:length(pmf_size)]
            pmf_rate <- pmf_rate[2:length(pmf_rate)]
          }
          
          # note how much we have added to vec_size
          x_perc_left <- x_perc_left - vec_size[j]
          
          # increase counter
          j <- j + 1
        }
        
        # now we add the size and weighted rate to the group compartment
        k_size[i] <- sum(vec_size)
        k_rate[i] <- sum(vec_size * vec_rate) / sum(vec_size)
        
        # verify that the weighted average makes sense
        if( !(vec_rate[1] <= k_rate[i] & k_rate[i] <= vec_rate[j-1]) ) { print("error; weighted rate calculation") }
        
        # store min and max
        k_min[i] <- vec_rate[1]
        k_max[i] <- vec_rate[j-1]
      }
      
      # remove PMF point from vectors if we have taken all the density
      if(pmf_size[1] <= 0 & length(pmf_size) > 1){
        pmf_size <- pmf_size[2:length(pmf_size)]
        pmf_rate <- pmf_rate[2:length(pmf_rate)]
      }
    }
    
    ## verify results
    # verify group sizes
    if( sum(k_size) < 1 ) { print("error; final group sizes do not add up to 1") }
    
    # verify mean partnership rate
    sum(k_rate * k_size) / sum(k_size)
    sum(data_pmf_tmp$y_pred * data_pmf_tmp$mean) / sum(data_pmf_tmp$mean)
    
    # format results
    df_out <- data.frame(city = "van",
                         time_pt = cur_time_pt,
                         sex_grp = paste("group_", 1:length(cdf_cutoff), sep = ""),
                         mean_rate = k_rate,
                         prop = k_size,
                         min_grp = k_min,
                         max_grp = k_max)
  }
  
  df_results_ls[[which(cur_time_pt == TIMEPTS)]] <- df_out
}

# Output ----
# verify outputs
for(cur_time_pt in names(df_results_ls)){
  print(paste(cur_time_pt, "============================================"))
  
  df_tmp <- df_results_ls[[cur_time_pt]]
  
  # check mean rates
  print(
    sprintf("Mean rate:      %s partners over P6M",
            round(sum(df_tmp$mean_rate * df_tmp$prop), 1))
  )
  
  # check minimum and max rates 
  print(
    sprintf("Rate range:     [%s-%s] partners over P6M",
            round(min(df_tmp$mean_rate), 1), round(max(df_tmp$mean_rate), 1))
  )
}
rm(df_tmp)

# round for output
df_results <- bind_rows(df_results_ls)
df_results$prop <- round(df_results$prop, 5)

df_results$prop

write_csv(df_results, "./misc-data-proc/outputs-cm/engage_van_actv_grps_nb_partn_p6m.csv")


