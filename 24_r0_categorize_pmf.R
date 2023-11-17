# Library and data ----
library(tidyverse)

data_pmf <- read_csv("./out/fitted-distributions/pmf_weighted_all_partn.csv")

## create the groups for the post-restrictions (main analysis) and
## pre-pandemic ('projection' of R_eff assuming those contact rates) results
CITIES <- c("mtl", "tor", "van")
TIME_PERIODS <- c("-Post-Restrictions", "-Pre-Pandemic")

df_results_ls <- list(post_main = NULL, pre_project = NULL)

for(time_period in TIME_PERIODS){
  print(paste(time_period, "========================"))
  
  # list to store results
  df_results <- list(mtl = NULL,
                     tor = NULL,
                     van = NULL)
  
  data_pmf_tmp <- subset(data_pmf, grepl(time_period, data_pt))
  
  for(cty in CITIES){
    # Categorize sex. actv. groups based on pre-specified percentiles ----
    ## Setup ----
    # percentiles at which to group
    cdf_cutoff <- c(.20, .30, .40, .50,  .60,  .70,  .80,  .90,  .92, .94,
                    .95, .96, .97, .98, .99, .992, .994, .996, .998, 1)
    
    data_city <- data_pmf_tmp %>% 
      filter(
        grepl(case_when(cty == "mtl" ~ "Montreal",
                        cty == "tor" ~ "Toronto",
                        cty == "van" ~ "Vancouver"),
              data_pt)
      )
    
    # fitted distributions of sexual partnerships
    pmf_size <- data_city$mean
    pmf_rate <- data_city$y_pred
    
    # sexual activity groups to create
    k_size <- vector("double", 0)
    k_rate <- vector("double", 0)
    
    # which points of the PMF are included
    k_min <- vector("double", 0)
    k_max <- vector("double", 0)
    
    # correct the PMF by reweighing to 1
    if( sum(pmf_size) < 1 ){
      print(sprintf("%s pmf was reweighted, previously added up to %s", cty, sum(pmf_size)))
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
    sum(data_city$y_pred * data_city$mean) / sum(data_city$mean)
    
    # format results
    df_out <- data.frame(city = cty,
                         sex_grp = paste("group_", 1:length(cdf_cutoff), sep = ""),
                         mean_rate = k_rate,
                         prop = k_size,
                         min_grp = k_min,
                         max_grp = k_max)
    
    df_results[[cty]] <- df_out
  }
  
  df_results_ls[[which(time_period == TIME_PERIODS)]] <- bind_rows(df_results)
}

# Output ----
# verify outputs
for(cty in CITIES){
  df_tmp_post <- subset(df_results_ls$post_main, city == cty)
  df_tmp_pre <- subset(df_results_ls$pre_project, city == cty)
  
  # check mean rates
  print(
    sprintf("Mean rate for %s:      post-restrictions -- %s    pre-pandemic -- %s",
            cty, round(sum(df_tmp_post$mean_rate * df_tmp_post$prop), 1),
            round(sum(df_tmp_pre$mean_rate * df_tmp_pre$prop), 1))
  )
  
  # check minimum and max rates 
  print(
    sprintf("Rate range for %s:     post-restrictions -- [%s-%s]    pre-pandemic -- [%s-%s]",
            cty, round(min(df_tmp_post$mean_rate), 1), round(max(df_tmp_post$mean_rate), 1),
            round(min(df_tmp_pre$mean_rate), 1), round(max(df_tmp_pre$mean_rate), 1))
  )
}
rm(df_tmp_post, df_tmp_pre)

# round for output
df_results_ls$post_main$prop; round(df_results_ls$post_main$prop, 4)
df_results_ls$pre_project$prop; round(df_results_ls$pre_project$prop, 4)

df_results_ls$post_main$prop <- round(df_results_ls$post_main$prop, 5)
df_results_ls$pre_project$prop <- round(df_results_ls$pre_project$prop, 5)

write_csv(df_results_ls$post_main, "./out/fitted-distributions-grouped/pmf_weighted_partn_p6m_20grps.csv")
write_csv(df_results_ls$pre_project, "./out/fitted-distributions-grouped/pmf_weighted_partn_p6m_20grps_prepandemic.csv")

