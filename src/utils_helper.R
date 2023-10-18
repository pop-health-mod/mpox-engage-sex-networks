# Variable grouping and formatting ----
## make indicator variables
# age; reference is 18-29
make_ind_age <- function(data){
  data %>% 
    mutate(age_30_39 = as.integer(age_grp == "30-39"),
           age_40_49 = as.integer(age_grp == "40-49"),
           age_50_59 = as.integer(age_grp == "50-59"),
           age_60_   = as.integer(age_grp == "60+"))
}

# partnership status; reference is not partnered
make_ind_rel <- function(data){
  data %>% mutate(rel_status = ifelse(is.na(rel_status),"unclear",rel_status))%>%
    mutate(rel_y_excl = as.integer(rel_status == "exclusive"),
           rel_y_open = as.integer(rel_status == "open"),
           rel_y_uncl = as.integer(rel_status == "unclear"))
}

# TODO MOVE TO CLEANING
# make city variable
make_city <- function(data){
  data %>% 
    mutate(city = case_when(grepl("^M", part_id) ~ "mtl",
                            grepl("^T", part_id) ~ "trt",
                            grepl("^V", part_id) ~ "van"),
           .before = 1)
}

# Simple data manipulation/computation ----
# for use in dplyr pipes
#' turn proportion to percentage and into character
round_prop <- function(x){
  as.character(round(x * 100, 0))
}

# get absolute minimum
abs_min <- function(x){
  min( abs(x) )
}

# Miscellaneous ----
# creates an empty 3-item list, one entry per city
create_city_list <- function(list_names = CITIES){
  ls_new <- vector("list", length = length(list_names))
  names(ls_new) <- list_names
  return(ls_new)
}

