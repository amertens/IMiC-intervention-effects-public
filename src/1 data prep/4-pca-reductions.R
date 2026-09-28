# =============================================================================
# src/1 data prep/4-pca-reductions.R
#
# Reduces each milk panel to one score per sample for Fig 1B: the first
# principal component of the panel's features (missing values median-imputed,
# near-zero-variance features dropped, features scaled), computed separately by
# study and visit and signed to correlate positively with the mean of the scaled
# features. Panels: macronutrients, micronutrients (ICP + fat-soluble vitamins),
# B-vitamins, HMOs, bioactive proteins, targeted metabolomics, microbiome (CLR
# ASVs) and untargeted metabolomics. Merges the scores with the covariates for
# src/2 analysis/3_adjusted_analysis_pca.R. Run it after
# 3_adjusted_analysis_microbiome.R, which writes the microbiome merge read here.
#
# Inputs:  data/clean milk data/{ELICIT,VITAL,MISAME}/*_{macronutrient,nutrient,ICP,HMO,MSDprotein,metabolite_b}.csv
#          data/milk/Allen_FSV_{ELICIT,VITAL,MISAME}.csv
#          data/clean milk data/merged_analysis_datasets_microbiome.RDS
#          data/{elicit,vital,misame}_untarget_metabolomics_pca.RDS (see the untargeted block below)
#          data/clean_baseline_covariates.RDS, data/clean_timevar_covariates.RDS
#          data/merged_analysis_datasets.RDS, metadata/milk_component.Rdata
# Outputs: data/pca_analysis_datasets.RDS
# [needs restricted data]
# =============================================================================

rm(list=ls())
source(paste0(here::here(),"/src/0-config.R"))
library(caret)

load(file=paste0(here::here(),"/metadata/milk_component.Rdata"))
d<-readRDS(paste0(here::here(),"/data/merged_analysis_datasets.RDS"))

load(file=paste0(here::here(),"/metadata/milk_component.Rdata"))
names(all_milk_components)
#"macro"        "micro"        "bvit"         "hmo"          "protein"      "metabolomics"

baseline <- readRDS(file=paste0(here::here(),"/data/clean_baseline_covariates.RDS"))
d_bmid <- readRDS(file=paste0(here::here(),"/data/clean_timevar_covariates.RDS")) #%>% distinct(bmid, visit)

#-------------------------------------------------------------------------------
# primary outcomes:  Macro- and micro-nutrients
#-------------------------------------------------------------------------------

#macronutrients
elicit_macro <- read.csv(paste0(here::here(),"/data/clean milk data/ELICIT/E_macronutrient.csv")) %>% rename(BMID=Unnamed..0, SUBJID=pid) %>% rename_with(tolower)
vital_macro <- read.csv(paste0(here::here(),"/data/clean milk data/VITAL/V_macronutrient.csv")) %>% rename(BMID=Unnamed..0, SUBJID=pid) %>% rename_with(tolower)
misame_macro <- read.csv(paste0(here::here(),"/data/clean milk data/MISAME/M_macronutrient.csv")) %>% rename(BMID=Unnamed..0, SUBJID=pid) %>% rename_with(tolower)

#micro-nutrients, b vitamins
elicit_bvit <- read.csv(paste0(here::here(),"/data/clean milk data/ELICIT/E_nutrient.csv")) %>% rename(BMID=Unnamed..0, SUBJID=pid) %>% rename_with(tolower)
vital_bvit <- read.csv(paste0(here::here(),"/data/clean milk data/VITAL/V_nutrient.csv")) %>% rename(BMID=Unnamed..0, SUBJID=pid) %>% rename_with(tolower)
misame_bvit <- read.csv(paste0(here::here(),"/data/clean milk data/MISAME/M_nutrient.csv")) %>% rename(BMID=Unnamed..0, SUBJID=pid) %>% rename_with(tolower)

elicit_micro <- read.csv(paste0(here::here(),"/data/clean milk data/ELICIT/E_ICP.csv")) %>% rename(BMID=Unnamed..0, SUBJID=pid) %>% rename_with(tolower)
vital_micro <- read.csv(paste0(here::here(),"/data/clean milk data/VITAL/V_ICP.csv")) %>% rename(BMID=Unnamed..0, SUBJID=pid) %>% rename_with(tolower)
misame_micro <- read.csv(paste0(here::here(),"/data/clean milk data/MISAME/M_ICP.csv")) %>% rename(BMID=Unnamed..0, SUBJID=pid) %>% rename_with(tolower)

# Fat-soluble vitamins (vitamin A, tocopherols), delivered as separate files;
# merged into the micronutrient panel below.
elicit_fsv <- read.csv(paste0(here::here(),"/data/milk/Allen_FSV_ELICIT.csv")) %>% rename(bmid=X) 
vital_fsv <- read.csv(paste0(here::here(),"/data/milk/Allen_FSV_VITAL.csv")) %>% rename(bmid=X) 
misame_fsv <- read.csv(paste0(here::here(),"/data/milk/Allen_FSV_MISAME.csv")) %>% rename(bmid=X) 

misame_fsv$bmid[!(misame_fsv$bmid %in% misame_micro$bmid)]
misame_micro$bmid[!(misame_micro$bmid %in% misame_fsv$bmid)]

dim(elicit_fsv)
dim(vital_fsv)
dim(misame_fsv)

elicit_fsv <- left_join(elicit_fsv, d_bmid %>% select(c("bmid", "visit","subjid")), by=c("bmid"))
vital_fsv <- left_join(vital_fsv, d_bmid %>% select(c("bmid", "visit","subjid")), by=c("bmid"))
misame_fsv <- left_join(misame_fsv, d_bmid %>% select(c("bmid", "visit","subjid")), by=c("bmid"))

# Drop two Misame FSV samples whose BMIDs match no sample in the other Misame
# milk panels (checked above; pn12-413 is a separate entry from pn12_413).
misame_fsv <- misame_fsv %>% filter(!(bmid %in% c("pn34_428", "pn12-413")))

elicit_micro <- full_join(elicit_micro, elicit_fsv, by=c("bmid", "visit","subjid")) 
vital_micro <- full_join(vital_micro, vital_fsv, by=c("bmid", "visit","subjid"))
misame_micro <- full_join(misame_micro, misame_fsv, by=c("bmid", "visit","subjid"))

# One Vital sample (A20010-LW) has no visit after the join; set it to visit 56.
vital_micro$visit[vital_micro$bmid=="A20010-LW"] <- 56

misame_micro[is.na(misame_micro$visit),]

misame_micro[misame_micro$bmid=="pn34_413",]

table(vital_micro$visit, useNA = "ifany")
table(misame_micro$visit, useNA = "ifany")
table(elicit_micro$visit, useNA = "ifany")

#-------------------------------------------------------------------------------
# secondary outcomes: HMO's, targeted proteins/bioactives
#-------------------------------------------------------------------------------

#HMO's
elicit_HMO <- read.csv(paste0(here::here(),"/data/clean milk data/ELICIT/E_HMO.csv")) %>% rename(BMID=Unnamed..0, SUBJID=pid) %>% rename_with(tolower)
vital_HMO <- read.csv(paste0(here::here(),"/data/clean milk data/VITAL/V_HMO.csv")) %>% rename(BMID=Unnamed..0, SUBJID=pid) %>% rename_with(tolower)
misame_HMO <- read.csv(paste0(here::here(),"/data/clean milk data/MISAME/M_HMO.csv")) %>% rename(BMID=Unnamed..0, SUBJID=pid) %>% rename_with(tolower)

#targeted proteins
elicit_protein <- read.csv(paste0(here::here(),"/data/clean milk data/ELICIT/E_MSDprotein.csv")) %>% rename(BMID=Unnamed..0, SUBJID=pid) %>% rename_with(tolower)
vital_protein <- read.csv(paste0(here::here(),"/data/clean milk data/VITAL/V_MSDprotein.csv")) %>% rename(BMID=Unnamed..0, SUBJID=pid) %>% rename_with(tolower)
misame_protein <- read.csv(paste0(here::here(),"/data/clean milk data/MISAME/M_MSDprotein.csv")) %>% rename(BMID=Unnamed..0, SUBJID=pid) %>% rename_with(tolower)

#-------------------------------------------------------------------------------
# tertiary outcomes: Targeted Metabolomics
#-------------------------------------------------------------------------------

#metabolomics - biocrate targeted
elicit_metabolomics <- read.csv(paste0(here::here(),"/data/clean milk data/ELICIT/E_metabolite_b.csv")) %>% rename(BMID=Unnamed..0, SUBJID=pid, ca_bio=CA, trp_bio=Trp) %>% rename_with(tolower)
vital_metabolomics <- read.csv(paste0(here::here(),"/data/clean milk data/VITAL/V_metabolite_b.csv")) %>% rename(BMID=Unnamed..0, SUBJID=pid, ca_bio=CA, trp_bio=Trp) %>% rename_with(tolower)
misame_metabolomics <- read.csv(paste0(here::here(),"/data/clean milk data/MISAME/M_metabolite_b.csv")) %>% rename(BMID=Unnamed..0, SUBJID=pid, ca_bio=CA, trp_bio=Trp) %>% rename_with(tolower)

#-------------------------------------------------------------------------------
# run PCAs
#-------------------------------------------------------------------------------

# First PC of one panel within one study/visit, signed so it correlates
# positively with the row mean of the scaled features.
run_pca <- function(df, modality, PC_num=1) {
  cat(df$visit[1],"\n")
  
  id <- df %>% select(bmid, subjid, arm, visit)
  df <- df %>% subset(., select=-c(bmid, subjid, arm, visit))
  df <- df %>% mutate(across(where(is.numeric), ~ if_else(is.na(.), median(., na.rm = TRUE), .)))
  #take average across columns
  mod_mean = as.numeric(rowMeans(scale(df), na.rm = T))
  if(length(nzv(df))>0){
    df <- df[,-nzv(df)]
  }
  
  df_pca <- prcomp(df, scale. = TRUE)$x[, PC_num]
  
  if(sign(cor(df_pca, mod_mean, use= "complete.obs")) < 0){
    df_pca = -1 * df_pca    
  }
  
  
  df_pca=data.frame(id, pca=df_pca)
  colnames(df_pca)[5] <- paste0(modality, "_pca")
  return(df_pca)
}

elicit_macro<-elicit_macro %>% group_by(visit) %>% do(run_pca(., "macro"))
vital_macro<-vital_macro %>% group_by(visit) %>% do(run_pca(., "macro"))
misame_macro<-misame_macro %>% group_by(visit) %>% do(run_pca(., "macro"))

elicit_bvit<-elicit_bvit %>% group_by(visit) %>% do(run_pca(., "bvit"))
vital_bvit<-vital_bvit %>% group_by(visit) %>% do(run_pca(., "bvit"))
misame_bvit<-misame_bvit %>% group_by(visit) %>% do(run_pca(., "bvit"))

elicit_micro<-elicit_micro %>% group_by(visit) %>% do(run_pca(., "micro"))
vital_micro<-vital_micro %>% group_by(visit) %>% do(run_pca(., "micro"))
misame_micro<-misame_micro %>% group_by(visit) %>% do(run_pca(., "micro"))

elicit_HMO<-elicit_HMO %>% group_by(visit) %>% do(run_pca(., "HMO"))
vital_HMO<-vital_HMO %>% group_by(visit) %>% do(run_pca(., "HMO"))
misame_HMO<-misame_HMO %>% group_by(visit) %>% do(run_pca(., "HMO"))

elicit_protein<-elicit_protein %>% group_by(visit) %>% do(run_pca(., "protein"))
vital_protein<-vital_protein %>% group_by(visit) %>% do(run_pca(., "protein"))
misame_protein<-misame_protein %>% group_by(visit) %>% do(run_pca(., "protein"))

elicit_metabolomics<-elicit_metabolomics %>% group_by(visit) %>% do(run_pca(., "metabolomics"))
vital_metabolomics<-vital_metabolomics %>% group_by(visit) %>% do(run_pca(., "metabolomics"))
misame_metabolomics<-misame_metabolomics %>% group_by(visit) %>% do(run_pca(., "metabolomics"))

#-------------------------------------------------------------------------------
# exploratory outcomes: microbiome and untargeted metabolomics
#-------------------------------------------------------------------------------

d_microbiome <- readRDS(paste0(here::here(),"/data/clean milk data/merged_analysis_datasets_microbiome.RDS"))
head(d_microbiome)

colnames(d_microbiome)

d_microbiome <- d_microbiome %>% subset(., select=-c(subjido, sex,                 mage,               
                                                     meducyrs,            mhtcm,               parity,              nperson,            
                                                     nrooms,              imp_water_src,       impfloor,            cookplac,           
                                                     dvseason,            dlvloc,              hhwealth,            hhfoodsecure,       
                                                     nperson_miss,        parity_miss,         meducyrs_miss, hhfoodsecure_miss,
                                                     Row.names, studyid, depth_below8000))
# elicit_
elicit_microbiome <- d_microbiome %>% filter(study=="Elicit") %>% subset(., select=-c(study))
misame_microbiome <- d_microbiome %>% filter(study=="Misame")  %>% subset(., select=-c(study))
vital_microbiome <- d_microbiome %>% filter(study=="Vital")  %>% subset(., select=-c(study))

elicit_microbiome<-elicit_microbiome %>% group_by(visit) %>% do(run_pca(., "microbiome"))
vital_microbiome<-vital_microbiome %>% group_by(visit) %>% do(run_pca(., "microbiome"))
misame_microbiome<-misame_microbiome %>% group_by(visit) %>% do(run_pca(., "microbiome"))

# Untargeted-metabolomics PCs. The commented block below produces
# data/{elicit,vital,misame}_untarget_metabolomics_pca.RDS (first PC of the
# untargeted features by study and visit), which are read right after it. It is
# kept commented out because it needs the full untargeted dataset loaded
# (d_metabolomics and Yvars from data/clean milk data/untargeted_metabolites.RData,
# written by src/2 analysis/3_adjusted_analysis_untargeted_metabolites.R) and is
# slow. To regenerate the three files, uncomment and run it once with that
# dataset loaded.
# load(here("data/clean milk data/untargeted_metabolites.RData"))
# dim(d_metabolomics)
# Yvars 
# d_metabolomics <- d_metabolomics %>% subset(., select=-c(subjido, sex,                 mage,               
#                                                          meducyrs,            mhtcm,               parity,              nperson,            
#                                                          nrooms,              imp_water_src,       impfloor,            cookplac,           
#                                                          dvseason,            dlvloc,              hhwealth,            hhfoodsecure,       
#                                                          nperson_miss,        parity_miss,         meducyrs_miss, hhfoodsecure_miss))
# 
# elicit_untarget_metabolomics <- d_metabolomics %>% filter(study=="Elicit") %>% subset(., select=-c(study))
# misame_untarget_metabolomics <- d_metabolomics %>% filter(study=="Misame")  %>% subset(., select=-c(study))
# vital_untarget_metabolomics <- d_metabolomics %>% filter(study=="Vital")  %>% subset(., select=-c(study))
# 
# elicit_untarget_metabolomics<-elicit_untarget_metabolomics %>% group_by(visit) %>% do(run_pca(., "untarget_metabolomics"))
# vital_untarget_metabolomics<-vital_untarget_metabolomics %>% group_by(visit) %>% do(run_pca(., "untarget_metabolomics"))
# misame_untarget_metabolomics<-misame_untarget_metabolomics %>% group_by(visit) %>% do(run_pca(., "untarget_metabolomics"))
# 
# colnames(elicit_untarget_metabolomics)[5] <- "untarget_metabolomics_pca"
# colnames(vital_untarget_metabolomics)[5] <- "untarget_metabolomics_pca"
# colnames(misame_untarget_metabolomics)[5] <- "untarget_metabolomics_pca"
# 
# 
# #save these datasets so later runs can skip this step
# saveRDS(elicit_untarget_metabolomics, file=paste0(here::here(),"/data/elicit_untarget_metabolomics_pca.RDS"))
# saveRDS(vital_untarget_metabolomics, file=paste0(here::here(),"/data/vital_untarget_metabolomics_pca.RDS"))
# saveRDS(misame_untarget_metabolomics, file=paste0(here::here(),"/data/misame_untarget_metabolomics_pca.RDS"))

elicit_untarget_metabolomics = readRDS(file=paste0(here::here(),"/data/elicit_untarget_metabolomics_pca.RDS")) %>% droplevels()
vital_untarget_metabolomics = readRDS(file=paste0(here::here(),"/data/vital_untarget_metabolomics_pca.RDS")) %>% droplevels()
misame_untarget_metabolomics = readRDS(file=paste0(here::here(),"/data/misame_untarget_metabolomics_pca.RDS")) %>% droplevels()

misame_untarget_metabolomics <- misame_untarget_metabolomics %>% mutate(subjid=as.numeric(subjid), 
                                                                        arm=factor(arm, levels = c( "Control","IFA/BEP", "BEP/IFA",  "BEP/BEP" )),
                                                                        arm=as.numeric(arm), visit=as.numeric(visit))

#-------------------------------------------------------------------------------
# Merge datasets
#-------------------------------------------------------------------------------

#merge misame
misame <- left_join(misame_macro, misame_micro, by=c("bmid","subjid","arm","visit"))
misame <- left_join(misame, misame_bvit, by=c("bmid","subjid","arm","visit"))
misame <- left_join(misame, misame_HMO, by=c("bmid","subjid","arm","visit"))
misame <- left_join(misame, misame_protein, by=c("bmid","subjid","arm","visit"))
misame <- left_join(misame, misame_metabolomics, by=c("bmid","subjid","arm","visit"))
misame <- left_join(misame, misame_microbiome %>% mutate(subjid=as.numeric(subjid), arm=as.numeric(arm), visit=as.numeric(visit)), by=c("bmid","subjid","arm","visit"))
misame <- left_join(misame, misame_untarget_metabolomics, by=c("bmid","subjid","arm","visit"))

unique(misame$subjid)
unique(misame_untarget_metabolomics$subjid)
unique(misame$arm)
unique(misame_untarget_metabolomics$arm)
class(misame$arm)
class(misame_untarget_metabolomics$arm)

head(misame)

#merge elicit
elicit <- left_join(elicit_macro, elicit_micro, by=c("bmid","subjid","arm","visit"))
elicit <- left_join(elicit, elicit_bvit, by=c("bmid","subjid","arm","visit"))
elicit <- left_join(elicit, elicit_HMO, by=c("bmid","subjid","arm","visit"))
elicit <- left_join(elicit, elicit_protein, by=c("bmid","subjid","arm","visit"))
elicit <- left_join(elicit, elicit_metabolomics, by=c("bmid","subjid","arm","visit"))
elicit <- left_join(elicit, elicit_microbiome %>% mutate(subjid=as.numeric(subjid), arm=as.numeric(arm), visit=as.numeric(visit)), by=c("bmid","subjid","arm","visit"))
elicit <- left_join(elicit, elicit_untarget_metabolomics %>% mutate(subjid=as.numeric(subjid), arm=as.numeric(arm), visit=as.numeric(visit)), by=c("bmid","subjid","arm","visit"))
head(elicit)

#merge vital
vital <- left_join(vital_macro, vital_micro, by=c("bmid","subjid","arm","visit"))
vital <- left_join(vital, vital_bvit, by=c("bmid","subjid","arm","visit"))
vital <- left_join(vital, vital_HMO, by=c("bmid","subjid","arm","visit"))
vital <- left_join(vital, vital_protein, by=c("bmid","subjid","arm","visit"))
vital <- left_join(vital, vital_metabolomics, by=c("bmid","subjid","arm","visit"))
vital <- left_join(vital, vital_microbiome %>% mutate(subjid=as.numeric(subjid), arm=as.numeric(arm), visit=as.numeric(visit)), by=c("bmid","subjid","arm","visit"))
vital <- left_join(vital, vital_untarget_metabolomics %>% mutate(subjid=as.numeric(subjid), arm=as.numeric(arm), visit=as.numeric(visit)), by=c("bmid","subjid","arm","visit"))

#-------------------------------------------------------------------------------
# Combine and save datasets and names
#-------------------------------------------------------------------------------

elicit <- elicit %>% mutate(study="Elicit")
misame <- misame %>% mutate(study="Misame")
vital <- vital %>% mutate(study="Vital")

bind_rows_auto_convert <- function(df1, df2) {
  # Identify columns that exist in both data frames but have different types
  common_cols <- intersect(names(df1), names(df2))
  mismatched_cols <- common_cols[sapply(common_cols, function(col) {
    !identical(typeof(df1[[col]]), typeof(df2[[col]]))
  })]
  
  # Convert mismatched columns to character
  for(col in mismatched_cols) {
    df1[[col]] <- as.character(df1[[col]])
    df2[[col]] <- as.character(df2[[col]])
  }
  
  # Bind rows
  res <- bind_rows(df1, df2)
}

d <- bind_rows_auto_convert(elicit, misame)
d <- bind_rows_auto_convert(d, vital)

d <- d %>% 
  mutate(milk_flag=1) %>%
  select(c("study","subjid", "arm", "visit",  "milk_flag", "bmid", ends_with("_pca"))) %>%
  rename(armcd=arm)

#merge baseline and time varying and merge 
d_bmid %>% distinct(studyid, subjid, subjido) %>% dim()
baseline %>% distinct(studyid, subjid, subjido) %>% dim()
dim(d_bmid)
dim(baseline)
d_bmid <- left_join(d_bmid, baseline, by = c("studyid","subjid","subjido"))

# Harmonize key types before joining with the milk scores
d_bmid <- d_bmid %>% mutate(subjid=as.numeric(subjid), armcd=as.numeric(armcd), visit=as.numeric(visit))

#create independent ID variable
d_bmid$clusterid <- d_bmid$siteid
d_bmid$clusterid[is.na(d_bmid$clusterid)] <- d_bmid$subjid[is.na(d_bmid$clusterid)]

#merge baseline and milk data
d_bmid$subjid <- as.numeric(d_bmid$subjid)
d_bmid$visit <- as.numeric(d_bmid$visit)
d <- left_join(d_bmid, d, by = c("subjid", "armcd", "visit", "bmid")) %>% filter(milk_flag==1) %>% subset(., select=-c(milk_flag))
head(d)

#clean Arm codes
d$arm<-as.character(d$arm)
d$arm[d$arm=="Placebo, Placebo"] <- "Control"
d$arm[d$arm=="IFA (preg) / IFA 6W postpartum"] <- "Control"
d$arm <- gsub("BEP\\+IFA \\(preg\\)","BEP",d$arm)
d$arm <- gsub("BEP\\+IFA 6W postpartum","BEP",d$arm)
d$arm <- gsub("Nutrient supplement","BEP",d$arm)
d$arm <- gsub("IFA 6W postpartum","IFA",d$arm)
d$arm <- gsub("IFA \\(preg\\)","IFA",d$arm)
d$arm <- gsub("Ex.BreastFeed","ExBf",d$arm)
d$arm <- gsub("Azithromycin","Az.",d$arm)
d$arm <- gsub("Nitazoxanide","",d$arm)
d$arm <- gsub("Nicotinamide","Nico",d$arm)
d$arm <- gsub(", Placebo","",d$arm)
d$arm <- gsub("Placebo, ","",d$arm)
d$arm <- gsub(" / ","/",d$arm)
d$arm <- gsub(" &",", ",d$arm)
d$arm <- gsub(", ","+",d$arm)
d$arm <- gsub("\\+$","",d$arm)
d$arm <- relevel(factor(d$arm), ref="Control")

table(d$arm)
table(d$armcd,d$arm)
colnames(d)

head(d)

saveRDS(d, file=paste0(here::here(),"/data/pca_analysis_datasets.RDS"))
