

#function to get measures closest to X months within a window
#d: data
#agem: age in months to get the measure by (assuming average months are 30.4067 days)
#window: window around the target age in months that the measure is allowed to fall 
#    default: 30.4167 days but that drops more observations than using a larger window like a month in either direction
#   NOTE!: if the window is too large, there is a danger of including a growth observation before the -omics sample collection
#  So if using this data, after merging to exposure data make sure to clean by dropping obs in the wrong time order based on child age
get_age_specific_measures <- function(d, agem=6, window=14){
  
  dsub <- d %>% group_by(studyid, subjid) %>%
    filter(abs(agedays-!!(agem)* 30.4167) == min(abs(agedays-!!(agem)* 30.4167)),
           abs(agedays-!!(agem)* 30.4167) <= window) %>%
    ungroup()
  
  return(dsub)
}


calc_centile_delta <- function(centile, agedays, agemax){
  centile <- centile[agedays < agemax]
  if(length(centile)>1){
    centile_delta <- max(centile, na.rm=T)-min(centile, na.rm=T)
  }else{
    centile_delta <- NA
  }
  return(centile_delta)
}



#function to calculate number of WHO centile spaces crossed:
calculateCentileCrossed <- function(measurement1, measurement2) {
  # Define the centile lines
  centiles <- c(95, 90, 75, 50, 25, 10, 5)
  
  # Find the centile index for each measurement
  findCentileIndex <- function(measurement) {
    sum(measurement <= centiles) + 1
  }
  
  index1 <- findCentileIndex(measurement1)
  index2 <- findCentileIndex(measurement2)
  
  # Calculate the number of centile spaces crossed (negative is a fall)
  return(index2 - index1)
}

#function to calculate different velocity measures
calc_velocity_measures <- function(d1, d2){
  
  d <- bind_rows(d1, d2)
  #drop kids without measures at both time points
  d <- d %>% group_by(studyid, subjid, subjido) %>%
    mutate(n_meas=n()) %>% filter(n_meas>1)
  
  # d <- d %>% arrange(studyid, subjid, subjido, agedays)
  # dim(d2)
  # dim(d %>% distinct(studyid, subjid, subjido))
  
  df <- d %>% group_by(studyid, subjid, subjido) %>%
    arrange(agedays, .by_group = TRUE) %>%
    mutate(
      age_diff = agedays - lag(agedays),
      lencm_delta = (lencm - lag(lencm))/age_diff,
      laz_delta = (haz - lag(haz))/age_diff,
      wtkg_delta = (wtkg - lag(wtkg))/age_diff,
      waz_delta = (waz - lag(waz))/age_diff,
      whz_delta = (whz - lag(whz))/age_diff,
      hcircm_delta = (hcircm - lag(hcircm))/age_diff,
      hcaz_delta = (hcaz - lag(hcaz))/age_diff,
      # Growth faltering at 3, 6, 12, 18, 24
      # a fall across 1 or more weight centile spaces, if birthweight was below the 9th centile
      # a fall across 2 or more weight centile spaces, if birthweight was between the 9th and 91st centiles
      # a fall across 3 or more weight centile spaces, if birthweight was above the 91st centile
      #when current weight is below the 2nd centile for age, whatever the birthweight.
      height_centile_vel = height_centile - lag(height_centile),
      weight_centile_vel = weight_centile - lag(weight_centile),
      height_centiles_crossed=mapply(calculateCentileCrossed,height_centile, first(height_centile)),
      weight_centiles_crossed=mapply(calculateCentileCrossed,weight_centile, first(weight_centile)),
      #Add components of growth faltering for Nolan's exploration
      growth_faltering1_lbw = 1*(birthweight_centile<9 & weight_centiles_crossed<0),
      growth_faltering2_1crossed = 1*(birthweight_centile<91 & birthweight_centile>=9 & weight_centiles_crossed < (-1)),
      growth_faltering3_2crossed = 1*(birthweight_centile>=91 & weight_centiles_crossed < (-2)),
      growth_faltering4_low_weight = 1*(weight_centile < 2),
      #calc growth faltering
      growth_faltering=ifelse(birthweight_centile<9 & weight_centiles_crossed<0 |
                                birthweight_centile<91 & birthweight_centile>=9 & weight_centiles_crossed < (-1) |
                                birthweight_centile>=91 & weight_centiles_crossed < (-2) |
                                weight_centile < 2, 1, 0),
      #"Thriving" Moving 2 centile spaces in a positive direction in the previous period regardless of anthropometric status
      thriving=ifelse(weight_centiles_crossed>=2,1,0),
      #Recovery: moving 2 centile spaces in a positive direction out of anthropometric deficit without relapse
      #Note: the relapse is not checked
      recovery=ifelse(weight_centiles_crossed>=2 & first(uwt==1) |
                        height_centiles_crossed>=2 & first(stunt==1) |
                        weight_centiles_crossed>=2 & first(wast==1), 1, 0)
    ) %>% filter(agedays!=min(agedays)) %>%
    ungroup() %>%
    select(studyid, subjid, subjido, lencm_delta,
           laz_delta, wtkg_delta, waz_delta, whz_delta, hcircm_delta, hcaz_delta,
           growth_faltering1_lbw, growth_faltering2_1crossed, growth_faltering3_2crossed, growth_faltering4_low_weight,
           height_centile_vel, weight_centile_vel,
           growth_faltering, thriving, recovery) 

  df <- df %>% arrange(studyid, subjid, subjido)

  return(df)
}




#Load and compile growth velocity standards
library(readxl)

who_girls_2mo <- read_excel(paste0(here::here(),"/data/metadata/growth standards/ttt_length_girls_2mon_z.xlsx")) %>% mutate(sex="female", agerange=2)
who_girls_3mo <- read_excel(paste0(here::here(),"/data/metadata/growth standards/ttt_length_girls_3mon_z.xlsx")) %>% mutate(sex="female", agerange=3)
who_girls_4mo <- read_excel(paste0(here::here(),"/data/metadata/growth standards/ttt_length_girls_4mon_z.xlsx")) %>% mutate(sex="female", agerange=4)
who_girls_6mo <- read_excel(paste0(here::here(),"/data/metadata/growth standards/ttt_length_girls_6mon_z.xlsx")) %>% mutate(sex="female", agerange=6)
who_boys_2mo <- read_excel(paste0(here::here(),"/data/metadata/growth standards/ttt_length_boys_2mon_z.xlsx")) %>% mutate(sex="male", agerange=2)
who_boys_3mo <- read_excel(paste0(here::here(),"/data/metadata/growth standards/ttt_length_boys_3mon_z.xlsx")) %>% mutate(sex="male", agerange=3)
who_boys_4mo <- read_excel(paste0(here::here(),"/data/metadata/growth standards/ttt_length_boys_4mon_z.xlsx")) %>% mutate(sex="male", agerange=4)
who_boys_6mo <- read_excel(paste0(here::here(),"/data/metadata/growth standards/ttt_length_boys_6mon_z.xlsx")) %>% mutate(sex="male", agerange=6)
who_data_length <- bind_rows( who_girls_2mo, who_girls_3mo, who_girls_4mo, who_girls_6mo,
                              who_boys_2mo, who_boys_3mo, who_boys_4mo, who_boys_6mo)

who_weight_girls_1mo <- read_excel(paste0(here::here(),"/data/metadata/growth standards/ttt-weight-girls-1mon-z.xlsx")) %>% mutate(sex="female", agerange=1)
who_weight_girls_2mo <- read_excel(paste0(here::here(),"/data/metadata/growth standards/ttt-weight-girls-2mon-z.xlsx")) %>% mutate(sex="female", agerange=2)
who_weight_girls_3mo <- read_excel(paste0(here::here(),"/data/metadata/growth standards/ttt-weight-girls-3mon-z.xlsx")) %>% mutate(sex="female", agerange=3)
who_weight_girls_4mo <- read_excel(paste0(here::here(),"/data/metadata/growth standards/ttt_weight_girls_4mon_z.xlsx")) %>% mutate(sex="female", agerange=4)
who_weight_girls_6mo <- read_excel(paste0(here::here(),"/data/metadata/growth standards/ttt_weight_girls_6mon_z.xlsx")) %>% mutate(sex="female", agerange=6)
who_weight_boys_1mo <- read_excel(paste0(here::here(),"/data/metadata/growth standards/ttt-weight-boys-1mon-z.xlsx")) %>% mutate(sex="male", agerange=1)
who_weight_boys_2mo <- read_excel(paste0(here::here(),"/data/metadata/growth standards/ttt-weight-boys-2mon-z.xlsx")) %>% mutate(sex="male", agerange=2)
who_weight_boys_3mo <- read_excel(paste0(here::here(),"/data/metadata/growth standards/ttt-weight-boys-3mon-z.xlsx")) %>% mutate(sex="male", agerange=3)
who_weight_boys_4mo <- read_excel(paste0(here::here(),"/data/metadata/growth standards/ttt_weight_boys_4mon_z.xlsx")) %>% mutate(sex="male", agerange=4)
who_weight_boys_6mo <- read_excel(paste0(here::here(),"/data/metadata/growth standards/ttt_weight_boys_6mon_z.xlsx")) %>% mutate(sex="male", agerange=6)
who_data_weight <- bind_rows(who_weight_girls_1mo, who_weight_girls_2mo, who_weight_girls_3mo, who_weight_girls_4mo, who_weight_girls_6mo,
                             who_weight_boys_1mo, who_weight_boys_2mo, who_weight_boys_3mo, who_weight_boys_4mo, who_weight_boys_6mo)




#clean/standardize charts
clean_who_vel_standards <- function(who_data){
  who_data$sex <- str_to_title(who_data$sex)
  
  who_data$Interval <- gsub(" \\u2013 ","-",who_data$Interval)
  who_data$Interval <- gsub(" \\u0080 ","",who_data$Interval)
  who_data$Interval <- gsub(" \\u0093 ","",who_data$Interval)
  who_data$Interval <- gsub(" ","",who_data$Interval)
  who_data$Interval <- gsub("–","-",who_data$Interval)
  
  #first interval
  who_data$interval1 <- str_split(who_data$Interval, "-", simplify = T)[,1]
  #convert to days
  who_data$interval1 <- as.numeric(gsub("4wks","0.9205469",who_data$interval1)) * 30.4167
  who_data$interval1[who_data$interval1==0] <- 1
  
  #second interval
  who_data$interval2 <- str_split(who_data$Interval, "-", simplify = T)[,2]
  #convert to days
  who_data$interval2 <- gsub("mo","",who_data$interval2)
  who_data$interval2 <- as.numeric(gsub("4wks","0.9205469",who_data$interval2)) * 30.4167
  return(who_data)
}


who_data_length <- clean_who_vel_standards(who_data_length)
who_data_weight <- clean_who_vel_standards(who_data_weight)

#Add Delta as 0 to the length velocity standards
who_data_length$Delta <- 0


write.csv(who_data_length, file=paste0(here::here(),"/data/metadata/growth standards/combined_length_vel_standards.csv"))
write.csv(who_data_weight, file=paste0(here::here(),"/data/metadata/growth standards/combined_weight_vel_standards.csv"))



#function to calculate WHO velocity Z-scores
calculate_WHO_velocity_zscore <- function(sex, age1, age2, length1, length2, who_data){
  if(!is.na(sex) & !is.na(age1) & !is.na(age2) & !is.na(length1) & !is.na(length2)){
    
    # Calculate the time interval between measurements (in months)
    time_interval <- (age2 - age1) / 30.4167
    
    # Calculate the linear velocity (in cm/month)
    velocity <- (length2 - length1)
    
    # Subset the WHO data for the given sex and find closest intervals
    
    relevant_data1 <- who_data %>% 
      filter(sex == !!sex) %>%
      arrange(abs(agerange - time_interval), abs(interval1 - age1)) %>%
      head(2)
    # Interpolate the Z-scores between the closest two intervals
    # Extract intervals and calculate weights
    interval_start <- relevant_data1$interval1 
    weights <- 1 - abs(interval_start - age1) / sum(abs(interval_start - age1))
    
    # Compute Z-scores for both intervals
    zscores <- (((velocity+relevant_data1$Delta) / relevant_data1$M)^relevant_data1$L - 1) / (relevant_data1$L * relevant_data1$S)
    
    # Weighted average of Z-scores
    if(sum(is.na(zscores))==0){
      weighted_zscore1 <- sum(zscores * weights)
    }else{
      weighted_zscore1=zscores[!is.na(zscores)]
    }
    
    relevant_data2 <- who_data %>% 
      filter(sex == !!sex) %>%
      arrange(abs(interval1 - age1),abs(agerange - time_interval)) %>%
      head(2)
    intervals <- relevant_data2$agerange
    weights <- 1 - abs(intervals - time_interval) / sum(abs(intervals - time_interval))
    zscores <- (((velocity+relevant_data2$Delta) / relevant_data2$M)^relevant_data2$L - 1) / (relevant_data2$L * relevant_data2$S)
    if(sum(is.na(zscores))==0){
      weighted_zscore2 <- sum(zscores * weights)
    }else{
      weighted_zscore2=zscores[!is.na(zscores)]
    }
    
    weighted_zscore = (weighted_zscore1+weighted_zscore2)/2
    
    
    if(identical(weighted_zscore, numeric(0))){
      weighted_zscore = NA
    }
    if(abs(weighted_zscore)>6 & !is.na(weighted_zscore)){
      weighted_zscore = NA
    }
    
    return(weighted_zscore)
    
  }else{
    return(NA)
  }
}


# #wrapper function
# 
# calc_who_velocity_measures <- function(d, var="vel_z"){
#   
#   res = d %>% 
#     rowwise() %>% 
#     mutate(vel_z = as.numeric(calculate_velocity_zscore(sex=sex,
#                                                             age1=agedays_birth,
#                                                             age2=agedays_3mo,
#                                                             length1=wtkg_birth*1000,
#                                                             length2=wtkg_3mo*1000,
#                                                             who_data=who_data_weight))) %>% ungroup() %>%
#     rename(!!(vel_z)=var)
#   
#   return(res)
# }


#tabulation function
customSummary <- function(df, group_var) {
  tab_strat <- df %>% group_by(!!sym(group_var)) %>%
    summarise(across(-c('subjid','subjido'), .fns = list(
      summary = ~ if (all(. %in% c(0, 1, NA)) && length(unique(.)) < 4) {
        ones <- sum(., na.rm = TRUE)
        total <- sum(!is.na(.))
        paste0( round(ones/total * 100, 2), "% (", ones, "/", total,")")
      }else{
        mean_sd <- c(mean = mean(., na.rm = TRUE), sd = sd(., na.rm = TRUE), n = sum(!is.na(.)))
        paste0(round(mean_sd[1], 2), " (", round(mean_sd[2], 1),", N=",mean_sd[3],")")
      }
    )), .groups = 'drop')
  tab_overall <- df %>%
    summarise(across(-c('studyid','subjid','subjido'), .fns = list(
      summary = ~ if (all(. %in% c(0, 1, NA)) && length(unique(.)) < 4) {
        ones <- sum(., na.rm = TRUE)
        total <- sum(!is.na(.))
        paste0( round(ones/total * 100, 2), "% (", ones, "/", total,")")
      }else{
        mean_sd <- c(mean = mean(., na.rm = TRUE), sd = sd(., na.rm = TRUE), n = sum(!is.na(.)))
        paste0(round(mean_sd[1], 2), " (", round(mean_sd[2], 1),", N=",mean_sd[3],")")
      }
    )))
  
  tab=bind_rows(tab_strat, tab_overall)
  tab$studyid[5] <- 'Combined'   
  tab[tab=='NaN (NA)'] <-''
  tab[tab=="NaN% (0/0)"] <-''
  
  colnames(tab) = gsub('_summary','',colnames(tab) )
  return(tab)
}
