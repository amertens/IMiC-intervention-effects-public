
rm(list=ls())
library(tidyverse)

#https://github.com/nutriverse/  for intergrowth and zscorer
library(zscorer)
library(growthstandards)

source(paste0(here::here(),"/functions/growth_outcome_functions.R"))


#-------------------------------------------------------------------------------
#load IMIC harmonized datasets
#-------------------------------------------------------------------------------

misame <- read.csv(paste0(here::here(),"/data/imic_harmonized/MISAME_3_IMiC_analysis.csv")) 
vital <- read.csv(paste0(here::here(),"/data/imic_harmonized/VITAL_Lactation_IMiC_analysis.csv")) 
vital_mamta <- read.csv(paste0(here::here(),"/data/imic_harmonized/vital_mamta_IMiC_analysis.csv"))


elicit <- read.csv(paste0(here::here(),"/data/imic_harmonized/ELICIT_IMiC_analysis.csv"))
child <- read.csv(paste0(here::here(),"/data/imic_harmonized/CHILD_IMiC_analysis.csv"))


#growth velocity standards
who_data_length <- read.csv(file=paste0(here::here(),"/data/metadata/growth standards/combined_length_vel_standards.csv"))
who_data_weight <- read.csv(file=paste0(here::here(),"/data/metadata/growth standards/combined_weight_vel_standards.csv"))

summary(misame$HCAZ)
summary(elicit$HCAZ)
summary(child$HCAZ)
summary(vital$HCAZ)


#-------------------------------------------------------------------------------
#save just growth data
#-------------------------------------------------------------------------------

dput( colnames(child ))


dfull <- bind_rows(elicit %>% mutate(SUBJIDO=as.character(SUBJIDO)) %>% 
                     select("STUDYID", "SUBJID",  "SUBJIDO", "VISITNUM", "VISIT", "AGEDAYS","SEX", "WTKG",    "LENCM",     "HCIRCM",  "MUACCM",  "WAZ",     "HAZ",     "WHZ",         "HCAZ",    "MUAZ"), 
                         misame %>% mutate(SUBJIDO=as.character(SUBJIDO))  %>% 
                     select("STUDYID", "SUBJID",  "SUBJIDO", "VISITNUM", "VISIT", "AGEDAYS","SEX","GAGEBRTH",  'BIRTHWT',   "WTKG",    "LENCM", "HCIRCM",   "MUACCM",  "WAZ",     "HAZ",     "WHZ","HCAZ","MUAZ"), 
                         vital  %>% mutate(SUBJIDO=as.character(SUBJIDO)) %>% 
                     select("STUDYID", "SUBJID",  "SUBJIDO", "VISITNUM", "VISIT", "AGEDAYS","SEX","GAGEBRTH",  'BIRTHWT',   "WTKG",    "LENCM","HCIRCM", "MUACCM",  "WAZ",     "HAZ",     "WHZ","HCAZ","MUAZ"),
                         child  %>% mutate(SUBJIDO=as.character(SUBJIDO)) %>% 
                     select("STUDYID", "SUBJID",  "SUBJIDO", "VISITNUM", "VISIT", "AGEDAYS","SEX","GAGEBRTH",  'BIRTHWT',   "WTKG",    "LENCM","HCIRCM", "WAZ",     "HAZ",     "WHZ",         "HCAZ"))

#Make dataset for first 6 months with just antho measures for April
dfull_6mo <- dfull %>% filter(AGEDAYS < 30.4167*7) %>% select(STUDYID, SUBJID, SUBJIDO, VISITNUM, VISIT, AGEDAYS, WAZ, HAZ, WHZ,MUAZ, HCAZ) %>% 
  filter(!is.na(WAZ) | !is.na(HAZ) | !is.na(WHZ)) %>% arrange(STUDYID, SUBJID, AGEDAYS)
head(dfull_6mo)

write.csv(dfull_6mo, file=paste0(here::here(),"/results/imic_all_growth_measures_combined.csv"))


#-------------------------------------------------------------------------------
# calculate static growth measures
#-------------------------------------------------------------------------------
colnames(dfull) <-  tolower(colnames(dfull))
d <- dfull %>% filter(!is.na(wtkg) |
                        !is.na(lencm) |
                        !is.na(hcircm) |
                        !is.na(muaccm)) %>%
              mutate(stunt =  1*(haz < -2),
                     sstunt = 1*(haz < -3),
                     wast =  1*(whz < -2),
                     swast = 1*(whz < -3),
                     uwt =  1*(waz < -2),
                     suwt = 1*(waz < -3),
                     low_hc =  1*(hcaz < -2),
                     s_low_hc = 1*(hcaz < -3),
                     muac_wast =  1*(muaz < -2),
                     muac_swast = 1*(muaz < -3),
                     #calculate centiles:
                     height_centile=pnorm(haz)*100,
                     weight_centile=pnorm(waz)*100)

table(!is.na(d$hcircm),!is.na(d$hcaz ), d$studyid)

#note! need to look up hcz and muac outliers

#count and drop outliers
outlier_tab <- d %>% group_by(studyid) %>%
  summarise(laz_outlier = sum(haz < -6 | haz > 6, na.rm=T),
            wlz_outlier = sum(whz < -5 | whz > 5, na.rm=T),
            waz_outlier = sum(waz < -6 | waz > 5, na.rm=T),
            hcz_outlier = sum(hcaz < -5 | hcaz > 5, na.rm=T),
            muaz_outlier = sum(muaz < -5 | muaz > 5, na.rm=T))
outlier_tab

#save for growth report
saveRDS(outlier_tab, file=paste0(here::here(),"/results/growth_outlier_exclusions_tab.RDS"))

#drop biologically implausible values
d$haz[d$haz < -6 | d$haz > 6] <- NA
d$whz[d$whz < -5 | d$whz > 5] <- NA
d$waz[d$waz < -6 | d$waz > 5] <- NA
d$hcaz[d$hcaz < -5 | d$hcaz > 5] <- NA
d$muaz[d$muaz < -5 | d$muaz > 5] <- NA

#-------------------------------------------------------------------------------
#check for length decreases beyond the Technical Error of Measurement 
#-------------------------------------------------------------------------------
#Calculate prevalence of height decrease between measurements
#Count a height decrease if height decreases more than WHO standard (2.8 x expert TEM of 0.29 =  0.812)
#Count a head circumference decrease if it decreases more than WHO standard (2.8 x expert TEM of 0.16 =  0.448)

#From ReliabilityofanthropometricmeasurementsintheWHOMulticentre GrowthReferenceStudy
#https://cdn.who.int/media/docs/default-source/child-growth/child-growth-standards/acta-paediatrica-supplement-on-the-who-child-growth-standards.pdf?sfvrsn=e8c31ab6_0

impl_loss <- d %>% group_by(studyid, subjid, subjido) %>% 
  arrange(agedays, .by_group = TRUE) %>%
  mutate(length_loss_flag = 1*(lencm + 0.812 < lag(lencm)),
         headcirc_loss_flag = 1*(hcircm + 0.448 < lag(hcircm)),
         length_traj_flag = ifelse(length_loss_flag==1|
                         lag(length_loss_flag)==1|
                         lead(length_loss_flag)==1|
                         lead(length_loss_flag,2)==1,1,0),
         headcirc_traj_flag = ifelse(headcirc_loss_flag==1|
                                     lag(headcirc_loss_flag)==1|
                                     lead(headcirc_loss_flag)==1 |
                                     lead(headcirc_loss_flag,2)==1,1,0)) %>%
  select(studyid, subjid, subjido,  agedays, lencm, hcircm, haz, whz, hcaz, length_loss_flag, headcirc_loss_flag, length_traj_flag,headcirc_traj_flag)
table(impl_loss$length_loss_flag)
table(impl_loss$headcirc_loss_flag)
prop.table(table(impl_loss$length_loss_flag))*100
prop.table(table(impl_loss$headcirc_loss_flag))*100
impl_length_loss_trajectories <- impl_loss[impl_loss$length_traj_flag==1,] %>% filter(!is.na(subjido)) %>% as.data.frame()


impl_length_loss_trajectories <- impl_loss[impl_loss$length_traj_flag==1,] %>% filter(!is.na(subjido)) %>% 
  select(studyid, subjid, subjido,  agedays, haz, whz, length_loss_flag, lencm) %>% 
  group_by(studyid, subjid,  subjido) %>%
  mutate(haz_deviation=abs(mean(haz, na.rm=T)-haz), drop_flag=ifelse((length_loss_flag ==1 & haz_deviation > lag(haz_deviation) |
         lead(length_loss_flag) ==1 & haz_deviation >= lead(haz_deviation)), 1,NA)) %>%
  as.data.frame() 


impl_headcir_loss_trajectories <- impl_loss[impl_loss$headcirc_traj_flag==1,] %>% filter(!is.na(subjido)) %>% 
  select(studyid, subjid, subjido,  agedays, hcaz, headcirc_loss_flag, hcircm) %>% 
  group_by(studyid, subjid,  subjido) %>%
  mutate(hcaz_deviation=abs(mean(hcaz, na.rm=T)-hcaz),
         drop_flag=ifelse((headcirc_loss_flag ==1 & hcaz_deviation > lag(hcaz_deviation) |
                  lead(headcirc_loss_flag) ==1 & hcaz_deviation >= lead(hcaz_deviation)), 1,NA)) %>%
  as.data.frame() 

saveRDS(impl_length_loss_trajectories, file=paste0(here::here(),"/results/growth_impl_length_loss_trajectories.RDS"))
saveRDS(impl_headcir_loss_trajectories, file=paste0(here::here(),"/results/growth_impl_headcir_loss_trajectories.RDS"))


#consider also using this package for data cleaning:
#https://carriedaymont.github.io/growthcleanr/articles/usage.html

#check the data for repeated observations and remove
dim(d)
d <- d %>% distinct()
dim(d)
#need to figure out why these are repeated in the raw data

#save birth weight centile for growth faltering definition
d <- d %>% group_by(studyid, subjid, subjido) %>% 
  arrange(agedays, .by_group = TRUE) %>%
  mutate(birthweight_centile = first(weight_centile),
         birthweight_centile = ifelse(first(agedays)>14,NA,birthweight_centile)) %>% ungroup()

#make a plot of the birthweight centile density by study
d %>% ggplot(aes(x=birthweight_centile)) + geom_density() + facet_wrap(~studyid)
#update plot to add pretty colors
d %>% ggplot(aes(x=birthweight_centile, fill=studyid)) + geom_density(alpha=.5) + theme_minimal() + theme(legend.position="none")

#-------------------------------------------------------------------------------
# calculate centile delta (highest minus lowest in a period)
#-------------------------------------------------------------------------------

d <- d %>% group_by(studyid, subjido) %>%
  arrange(agedays, .by_group = TRUE) %>%
  mutate(height_centile_delta_3mo=calc_centile_delta(height_centile, agedays, 30.4167*3+14),
         height_centile_delta_6mo=calc_centile_delta(height_centile, agedays, 30.4167*6+14),
         height_centile_delta_12mo=calc_centile_delta(height_centile, agedays, 30.4167*13),
         height_centile_delta_18mo=calc_centile_delta(height_centile, agedays, 30.4167*19),
         height_centile_delta_24mo=calc_centile_delta(height_centile, agedays, 30.4167*25),
         weight_centile_delta_3mo=calc_centile_delta(weight_centile, agedays, 30.4167*3+14),
         weight_centile_delta_6mo=calc_centile_delta(weight_centile, agedays, 30.4167*6+14),
         weight_centile_delta_12mo=calc_centile_delta(weight_centile, agedays, 30.4167*13),
         weight_centile_delta_18mo=calc_centile_delta(weight_centile, agedays, 30.4167*19),
         weight_centile_delta_24mo=calc_centile_delta(weight_centile, agedays, 30.4167*25)) %>% filter(!is.na(subjido))

summary(d$height_centile_delta_3mo)
summary(d$height_centile_delta_6mo)
summary(d$height_centile_delta_12mo)
summary(d$height_centile_delta_18mo)
summary(d$height_centile_delta_24mo)


#-------------------------------------------------------------------------------
# age specific data for primary outcomes
#-------------------------------------------------------------------------------

colnames(child ) <- tolower(colnames(child ))
agem=3
df <- child%>% group_by(studyid, subjid) %>%
  filter(abs(agedays-!!(agem)* 30.4167) == min(abs(agedays-!!(agem)* 30.4167)),
         abs(agedays-!!(agem)* 30.4167) <= 30.4167) %>%
  ungroup()
plot(hist(df$agedays))
length(unique(df$subjid))


d_birth = get_age_specific_measures(d, agem=0)
d3 = get_age_specific_measures(d, agem=3, window=30.4167)
d6 = get_age_specific_measures(d, agem=6, window=30.4167)
d12 = get_age_specific_measures(d, agem=12, window=30.4167) %>% filter(studyid!='MISAME-3') #use 10 month measure for misame
d18 = get_age_specific_measures(d, agem=18, window=30.4167)
d24 = get_age_specific_measures(d, agem=24, window=30.4167)

table(d3$studyid)
table(d6$studyid)

#for misame: use 10 month measurements as there aren't 12 month measurements
d10 = get_age_specific_measures(d, agem=10, window=30.4167) %>% filter(studyid=='MISAME-3')
d12 = bind_rows(d10, d12)



#examine total unique children versus age-specific measures
d %>% group_by(studyid) %>% distinct(subjido) %>% summarize(n())
d_birth %>% group_by(studyid) %>% distinct(subjido) %>% summarize(n())
d3 %>% group_by(studyid) %>% distinct(subjido) %>% summarize(n())
d6 %>% group_by(studyid) %>% distinct(subjido) %>% summarize(n())
d12 %>% group_by(studyid) %>% distinct(subjido) %>% summarize(n())
d18 %>% group_by(studyid) %>% distinct(subjido) %>% summarize(n())
d24 %>% group_by(studyid) %>% distinct(subjido) %>% summarize(n())

#-------------------------------------------------------------------------------
# Calculate growth velocity measures (other than WHO Z-scores)
#-------------------------------------------------------------------------------

d_birth_3 <- calc_velocity_measures(d_birth, d3) %>% select(-starts_with("growth_faltering"))
d_3_6 <- calc_velocity_measures(d3, d6) %>% select(-starts_with("growth_faltering"))
d_6_12 <- calc_velocity_measures(d6, d12) %>% select(-starts_with("growth_faltering"))
d_12_18 <- calc_velocity_measures(d12, d18) %>% select(-starts_with("growth_faltering"))

prop.table(table(d_3_6$recovery))*100
prop.table(table(d_6_12$recovery))*100

#do birth to 6 months to get faltering measure and velocity measures
d_birth_6 <- calc_velocity_measures(d_birth, d6)

d1=d_birth
d2=d6

table(d_birth_6$growth_faltering)
prop.table(table(d_birth_6$growth_faltering))
table(d_birth_6$studyid, d_birth_6$growth_faltering)
table(d_birth_6$growth_faltering, d_birth_6$growth_faltering1_lbw==1)
table(d_birth_6$growth_faltering, d_birth_6$growth_faltering2_1crossed==1)
table(d_birth_6$growth_faltering, d_birth_6$growth_faltering3_2crossed==1)
table(d_birth_6$growth_faltering, d_birth_6$growth_faltering4_low_weight==1)

#merge into the later time point datasets
dim(d3)
dim(d_birth_3)
d3 <- left_join(d3, d_birth_3, by=c('studyid', 'subjid', 'subjido'))
d6 <- left_join(d6, d_3_6, by=c('studyid', 'subjid', 'subjido'))
d12 <- left_join(d12, d_6_12, by=c('studyid', 'subjid', 'subjido'))
d18 <- left_join(d18, d_12_18, by=c('studyid', 'subjid', 'subjido'))

#check merge
summary(d_3_6$laz_delta)
summary(d6$laz_delta)



d_birth_6 <- left_join(d3 %>% select(studyid,subjid,subjido,sex,gagebrth,birthwt,height_centile_delta_3mo,  height_centile_delta_6mo,  
                                     height_centile_delta_12mo, height_centile_delta_18mo, height_centile_delta_24mo, weight_centile_delta_3mo,
                                     weight_centile_delta_6mo,  weight_centile_delta_12mo, weight_centile_delta_18mo, weight_centile_delta_24mo),
                       d_birth_6, by=c('studyid', 'subjid', 'subjido'))



#-------------------------------------------------------------------------------
# make single growth outcomes dataset for analysis
#-------------------------------------------------------------------------------

#make a longform dataset of just the target ages
d_age_sub <- bind_rows(d_birth %>% mutate(agegroup="birth"), 
                       d_birth_6 %>% mutate(agegroup="birth_6mo"), 
                       d3 %>% mutate(agegroup="3mo"), 
                       d6 %>% mutate(agegroup="6mo"), 
                       d12 %>% mutate(agegroup="12mo"), 
                       d18 %>% mutate(agegroup="18mo"),
                       d24 %>% mutate(agegroup="24mo")) %>% 
  #drop rare age strata in specific studies (<25 measures at an age)
  #   -drops child at 6 (N=24) and 18 months (n=3) and vital (n=1)
  group_by(studyid, agegroup) %>% mutate(Nkids=n()) %>%
  filter(Nkids >= 25) %>% subset(., select=-c(Nkids, weight_centile, height_centile, birthweight_centile)) %>%
  ungroup() 

dim(d_age_sub)
dim(d_age_sub %>% distinct(studyid, subjido))

#save longform for plotting
table(d_age_sub$agegroup,  is.na(d_age_sub$growth_faltering))

saveRDS(d_age_sub, file=paste0(here::here(),"/results/imic_growth_outcomes_dataset_long.RDS"))

table(d_age_sub$agegroup,  is.na(d_age_sub$growth_faltering1_lbw))

#reshape to wide form
static_variables = c('studyid', 'subjid', 'subjido',"sex","gagebrth","birthwt",
                     # "growth_faltering1_lbw",  "growth_faltering2_1crossed",
                     # "growth_faltering3_2crossed",  "growth_faltering4_low_weight",
                     "height_centile_delta_3mo",  "height_centile_delta_6mo",  
                     "height_centile_delta_12mo", "height_centile_delta_18mo", "height_centile_delta_24mo", "weight_centile_delta_3mo",
                     "weight_centile_delta_6mo",  "weight_centile_delta_12mo", "weight_centile_delta_18mo", "weight_centile_delta_24mo")


colnames(d_age_sub)
d_wide<- d_age_sub %>%
  #pivot to wide 
  pivot_wider(id_cols = all_of(static_variables), 
              names_from = agegroup, 
              values_from=select(d_age_sub, -all_of(static_variables)) %>% names())
colnames(d_wide)

d_wide <- d_wide %>% subset(., select = -c(agedays_birth_6mo, 
                          height_centile_vel_birth_6mo, weight_centile_vel_birth_6mo,recovery_birth_6mo,
                          lencm_birth_6mo, wtkg_birth_6mo, hcircm_birth_6mo, muaccm_birth_6mo,
                          visitnum_birth_6mo, visit_birth_6mo))

#drop columns where all values are NA
d_wide <- d_wide[,colSums(is.na(d_wide))<nrow(d_wide)]


dim(d_wide)
dim(d_wide %>% distinct(studyid, subjido))
temp <- d_wide %>% select('studyid', 'subjid', 'subjido',starts_with("growth_faltering"))

#-------------------------------------------------------------------------------
# Calculate WHO Z-scores
#-------------------------------------------------------------------------------

d_wide = d_wide %>% 
  rowwise() %>% 
  mutate(
    len_vel_z_birth_3mo = as.numeric(calculate_WHO_velocity_zscore(sex=sex,age1=agedays_birth,age2=agedays_3mo,length1=lencm_birth,length2=lencm_3mo, who_data=who_data_length)),
    len_vel_z_birth_6mo = as.numeric(calculate_WHO_velocity_zscore(sex=sex,age1=agedays_birth,age2=agedays_6mo,length1=lencm_birth,length2=lencm_6mo, who_data=who_data_length)),
    len_vel_z_3mo_6mo = as.numeric(calculate_WHO_velocity_zscore(sex=sex,age1=agedays_3mo,age2=agedays_6mo,length1=lencm_3mo,length2=lencm_6mo, who_data=who_data_length)),
    len_vel_z_6mo_12mo = as.numeric(calculate_WHO_velocity_zscore(sex=sex,age1=agedays_6mo,age2=agedays_12mo,length1=lencm_6mo,length2=lencm_12mo, who_data=who_data_length)),
    
    wt_vel_z_birth_3mo = as.numeric(calculate_WHO_velocity_zscore(sex=sex,age1=agedays_birth,age2=agedays_3mo,length1=wtkg_birth*1000,length2=wtkg_3mo*1000, who_data=who_data_weight)),
    wt_vel_z_birth_6mo = as.numeric(calculate_WHO_velocity_zscore(sex=sex,age1=agedays_birth,age2=agedays_6mo,length1=wtkg_birth*1000,length2=wtkg_6mo*1000, who_data=who_data_weight)),
    wt_vel_z_3mo_6mo = as.numeric(calculate_WHO_velocity_zscore(sex=sex,age1=agedays_3mo,age2=agedays_6mo,length1=wtkg_3mo*1000,length2=wtkg_6mo*1000, who_data=who_data_weight)),
    wt_vel_z_6mo_12mo = as.numeric(calculate_WHO_velocity_zscore(sex=sex,age1=agedays_6mo,age2=agedays_12mo,length1=wtkg_6mo*1000,length2=wtkg_12mo*1000, who_data=who_data_weight))
    ) %>% ungroup() 

summary(d_wide$len_vel_z_birth_3mo)
summary(d_wide$len_vel_z_birth_6mo)
summary(d_wide$len_vel_z_3mo_6mo)
summary(d_wide$len_vel_z_6mo_12mo)

summary(d_wide$wt_vel_z_birth_3mo)
summary(d_wide$wt_vel_z_birth_6mo)
summary(d_wide$wt_vel_z_3mo_6mo)
summary(d_wide$wt_vel_z_6mo_12mo)

#-------------------------------------------------------------------------------
# calculate small-vulnerable-newborns from growth measures
# SVN is preterm, LBW, or SGA 
#-------------------------------------------------------------------------------

#preterm is <37 weeks GA
d_wide$preterm <- 1*(d_wide$gagebrth < 37*7)

#LBW is <2500 g
d_wide$lbw <- 1*(d_wide$birthwt < 2500)

#calc WAZ centile for GA
table(d_wide$sex)
table(is.na(d_wide$sex))
d_wide$waz_GA <- pnorm(igb_value2zscore(d_wide$gagebrth, d_wide$birthwt/1000, var = "wtkg", sex = d_wide$sex))*100
d_wide$sga <- 1*(d_wide$waz_GA < 10) #sga is <10th centile
d_wide$sga[is.na(d_wide$waz_GA)] <- NA

table(d_wide$sga)
table(d_wide$lbw)
table(d_wide$preterm)

d_wide$svn <- 1*(d_wide$sga==1 | d_wide$lbw==1 | d_wide$preterm==1)
table(d_wide$svn)
table(d_wide$studyid, d_wide$svn)
table(d_wide$studyid,is.na(d_wide$svn))
table(d_wide$studyid,is.na(d_wide$gagebrth))
table(d_wide$studyid,is.na(d_wide$birthwt))

#-------------------------------------------------------------------------------
# tabulate growth measures
#-------------------------------------------------------------------------------

#drop growth measures not calculated and order variables
colnames(d_wide)
d_wide <- d_wide %>% select(-matches("_delta_birth"),-matches("_delta_24mo"),
                            -matches('wast_birth'), -matches('stunt_birth'),  -matches('uwt_birth'), 
                            -matches('agegroup'), -starts_with('visit')) %>%
  select("studyid","subjid","subjido",starts_with('agedays_'),"svn", starts_with('haz_'), 
         starts_with('hcaz_'), starts_with('muaz_'), starts_with('baz_'), 
         matches('stunt'), matches('wast'), matches('uwt'), 
         matches('low_hc'), starts_with('growth_faltering'), starts_with('thriving'), starts_with('recovery'), 
         matches('centile'), -matches("waz_GA"),
         everything())


#tabulate variables
growth_tab = customSummary(df=d_wide, group_var='studyid')

#save dataset and table
saveRDS(growth_tab, file=paste0(here::here(),"/results/growth_summary_tab.RDS"))
write.csv(d_wide, file=paste0(here::here(),"/data/imic_full_growth_outcomes_dataset.csv"))


#subset to just outcome variables (and ages)
colnames(d_wide)

primary_outcomes <- c("haz_3mo", "haz_6mo", "whz_3mo",  "whz_6mo")
secondary_outcomes <- c("waz_3mo", "waz_6mo", "hcaz_3mo" ,"hcaz_6mo",
                        "len_vel_z_birth_6mo", "len_vel_z_3mo_6mo", 
                        "wt_vel_z_birth_6mo", "wt_vel_z_3mo_6mo",
                        "growth_faltering_birth_6mo")
growth_faltering_components <- c("growth_faltering1_lbw_birth_6mo",  "growth_faltering2_1crossed_birth_6mo",  
                                 "growth_faltering3_2crossed_birth_6mo",  "growth_faltering4_low_weight_birth_6mo")

analysis_df_vars <- c("studyid", "subjid", "subjido", "agedays_3mo", "agedays_6mo", "waz_birth", primary_outcomes, secondary_outcomes, growth_faltering_components, "svn")

#growth_faltering_components= ("growth_faltering1_lbw", "growth_faltering2_1crossed" , "growth_faltering3_2crossed", "growth_faltering4_low_weight")


dY <- d_wide %>% select(all_of(analysis_df_vars))

# double check for duplicates, should be 0
sum(duplicated(dY$subjido)) 

#save outcomes for Nolans analysis
write.csv(dY, file=paste0(here::here(),"/results/imic_growth_outcomes_dataset.csv"))


#now summarize with TEM exclusions

# saveRDS(impl_length_loss_trajectories, file=paste0(here::here(),"/results/growth_impl_length_loss_trajectories.RDS"))
# saveRDS(impl_headcir_loss_trajectories, file=paste0(here::here(),"/results/growth_impl_headcir_loss_trajectories.RDS"))

summary(impl_length_loss_trajectories$haz)
summary(impl_length_loss_trajectories$haz[is.na(impl_length_loss_trajectories$drop_flag)])


summary(impl_headcir_loss_trajectories$hcaz)
summary(impl_headcir_loss_trajectories$hcaz[is.na(impl_headcir_loss_trajectories$drop_flag)])





