
rm(list=ls())
source(paste0(here::here(),"/src/0-config.R"))

#-------------------------------------------------------------------------------
# baseline covariates and ID variables
#-------------------------------------------------------------------------------

misame_id <- read.csv(paste0(here::here(),"/data/imic_harmonized/MISAME_3_IMiC_analysis.csv"))  %>% distinct( SUBJID, BMID, ARM, VISIT) %>% filter(BMID!="")  %>%
  mutate(BMID=gsub("_lama","",BMID))
vital_id <- read.csv(paste0(here::here(),"/data/imic_harmonized/VITAL_Lactation_IMiC_analysis.csv")) %>% distinct( SUBJID, BMID, ARM, VISIT) %>% filter(BMID!="")
elicit_id <- read.csv(paste0(here::here(),"/data/imic_harmonized/ELICIT_IMiC_analysis.csv")) %>% distinct( SUBJID, BMID, ARM, VISIT) %>% filter(BMID!="")

colnames(misame_id) <- tolower(colnames(misame_id))
colnames(vital_id) <- tolower(colnames(vital_id))
colnames(elicit_id) <- tolower(colnames(elicit_id))


baseline <- readRDS(file=paste0(here::here(),"/data/clean_baseline_covariates.RDS"))
d_bmid <- readRDS(file=paste0(here::here(),"/data/clean_timevar_covariates.RDS"))

#-------------------------------------------------------------------------------
# primary outcomes:  Macro- and micro-nutrients
#-------------------------------------------------------------------------------

#macronutrients
elicit_macro <- read.csv(paste0(here::here(),"/data/clean milk data/ELICIT/E_macronutrient.csv")) %>% rename(BMID=Unnamed..0, SUBJID=pid) %>% rename_with(tolower)
vital_macro <- read.csv(paste0(here::here(),"/data/clean milk data/VITAL/V_macronutrient.csv")) %>% rename(BMID=Unnamed..0, SUBJID=pid) %>% rename_with(tolower)
misame_macro <- read.csv(paste0(here::here(),"/data/clean milk data/MISAME/M_macronutrient.csv")) %>% rename(BMID=Unnamed..0, SUBJID=pid) %>% rename_with(tolower)
child_macro <- read.csv(paste0(here::here(),"/data/clean milk data/CHILD/C_macronutrient.csv")) %>% rename(BMID=Unnamed..0, SUBJID=pid) %>% rename_with(tolower)


median(elicit_macro$protein)
median(vital_macro$protein)
median(misame_macro$protein)
median(child_macro$protein)

elicit_macro_pca <- prcomp(elicit_macro[,-1] %>% mutate(across(where(is.numeric), ~ if_else(is.na(.), median(., na.rm = TRUE), .))), scale. = TRUE)$x[, 1]

elicit_macro_pca <- prcomp(elicit_macro[,-1] %>% mutate(across(where(is.numeric), ~ if_else(is.na(.), median(., na.rm = TRUE), .))), scale. = TRUE)$x[, 1]

elicit_macro_pca <- prcomp(elicit_macro[,-1] %>% mutate(across(where(is.numeric), ~ if_else(is.na(.), median(., na.rm = TRUE), .))), scale. = TRUE)$x[, 1]

#micro-nutrients, b vitamins
elicit_bvit <- read.csv(paste0(here::here(),"/data/clean milk data/ELICIT/E_nutrient.csv")) %>% rename(BMID=Unnamed..0, SUBJID=pid) %>% rename_with(tolower)
vital_bvit <- read.csv(paste0(here::here(),"/data/clean milk data/VITAL/V_nutrient.csv")) %>% rename(BMID=Unnamed..0, SUBJID=pid) %>% rename_with(tolower)
misame_bvit <- read.csv(paste0(here::here(),"/data/clean milk data/MISAME/M_nutrient.csv")) %>% rename(BMID=Unnamed..0, SUBJID=pid) %>% rename_with(tolower)

elicit_micro <- read.csv(paste0(here::here(),"/data/clean milk data/ELICIT/E_ICP.csv")) %>% rename(BMID=Unnamed..0, SUBJID=pid) %>% rename_with(tolower)
vital_micro <- read.csv(paste0(here::here(),"/data/clean milk data/VITAL/V_ICP.csv")) %>% rename(BMID=Unnamed..0, SUBJID=pid) %>% rename_with(tolower)
misame_micro <- read.csv(paste0(here::here(),"/data/clean milk data/MISAME/M_ICP.csv")) %>% rename(BMID=Unnamed..0, SUBJID=pid) %>% rename_with(tolower)

# #merge in Vit-A data added later
# d_visit <- d_bmid %>% distinct(bmid, visit)
# 
# d_visit[d_visit$bmid %in% c("pn34_428", "pn12-413"),]


elicit_fsv <- read.csv(paste0(here::here(),"/data/milk/Allen_FSV_ELICIT.csv")) %>% rename(bmid=X) 
vital_fsv <- read.csv(paste0(here::here(),"/data/milk/Allen_FSV_VITAL.csv")) %>% rename(bmid=X) 
misame_fsv <- read.csv(paste0(here::here(),"/data/milk/Allen_FSV_MISAME.csv")) %>% rename(bmid=X) 

#add in missing IDs
elicit_fsv <- left_join(elicit_fsv, d_bmid, by=c("bmid"))
vital_fsv <- left_join(vital_fsv, d_bmid, by=c("bmid"))
misame_fsv <- left_join(misame_fsv, d_bmid, by=c("bmid"))
                        
head(elicit_fsv)
head(vital_fsv)
head(misame_fsv)

#TEMPORARY
#Drop IDs missing from other datasets
misame_fsv <- misame_fsv %>% filter(!(bmid %in% c("pn34_428", "pn12-413")))




#-------------------------------------------------------------------------------
# secondary outcomes: HMO’s, targeted proteins/bioactives
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

#metabolomics - biocrates targeted
elicit_metabolomics <- read.csv(paste0(here::here(),"/data/clean milk data/ELICIT/E_metabolite_b.csv")) %>% rename(BMID=Unnamed..0, SUBJID=pid, ca_bio=CA, trp_bio=Trp) %>% rename_with(tolower)
vital_metabolomics <- read.csv(paste0(here::here(),"/data/clean milk data/VITAL/V_metabolite_b.csv")) %>% rename(BMID=Unnamed..0, SUBJID=pid, ca_bio=CA, trp_bio=Trp) %>% rename_with(tolower)
misame_metabolomics <- read.csv(paste0(here::here(),"/data/clean milk data/MISAME/M_metabolite_b.csv")) %>% rename(BMID=Unnamed..0, SUBJID=pid, ca_bio=CA, trp_bio=Trp) %>% rename_with(tolower)
#note ca and trp are renamed because they are variable names in other datasets too





#-------------------------------------------------------------------------------
# Save milk names by lab dataset
#-------------------------------------------------------------------------------

ls()

#get_milk_names <- function(x){colnames(x)[colnames(x)!="bmid"]}
get_milk_names <- function(x){colnames(x)[!(colnames(x) %in% c("bmid","subjid","arm","visit"))]}

all_milk_components <-elicit_milk_components <- misame_milk_components <- vital_milk_components <- list()

elicit_milk_components$macro <- get_milk_names(elicit_macro)
elicit_milk_components$micro <- c(get_milk_names(elicit_micro), "g.tocopherol", "a.tocopherol", "vitamin.a")
elicit_milk_components$bvit <- get_milk_names(elicit_bvit)
elicit_milk_components$hmo <- get_milk_names(elicit_HMO)
elicit_milk_components$protein <- get_milk_names(elicit_protein)
elicit_milk_components$metabolomics <- get_milk_names(elicit_metabolomics) 

vital_milk_components$macro <- get_milk_names(vital_macro)
vital_milk_components$micro <- c(get_milk_names(vital_micro), "g.tocopherol", "a.tocopherol", "vitamin.a")
vital_milk_components$bvit <- get_milk_names(vital_bvit)
vital_milk_components$hmo <- get_milk_names(vital_HMO)
vital_milk_components$protein <- get_milk_names(vital_protein)
vital_milk_components$metabolomics <- get_milk_names(vital_metabolomics) 

misame_milk_components$macro <- get_milk_names(misame_macro)
misame_milk_components$micro <- c(get_milk_names(misame_micro), "g.tocopherol", "a.tocopherol", "vitamin.a")
misame_milk_components$bvit <- get_milk_names(misame_bvit)
misame_milk_components$hmo <- get_milk_names(misame_HMO)
misame_milk_components$protein <- get_milk_names(misame_protein)
misame_milk_components$metabolomics <- get_milk_names(misame_metabolomics) 

all_milk_components$macro <- unique(c(elicit_milk_components$macro, vital_milk_components$macro, misame_milk_components$macro))
all_milk_components$micro <- unique(c(elicit_milk_components$micro, vital_milk_components$micro, misame_milk_components$micro))
all_milk_components$bvit <- unique(c(elicit_milk_components$bvit, vital_milk_components$bvit, misame_milk_components$bvit))
all_milk_components$hmo <- unique(c(elicit_milk_components$hmo, vital_milk_components$hmo, misame_milk_components$hmo))
all_milk_components$protein <- unique(c(elicit_milk_components$protein, vital_milk_components$protein, misame_milk_components$protein))
all_milk_components$metabolomics <- unique(c(elicit_milk_components$metabolomics, vital_milk_components$metabolomics, misame_milk_components$metabolomics))

#-------------------------------------------------------------------------------
# Merge datasets
#-------------------------------------------------------------------------------




#merge misame
misame <- full_join(misame_macro, misame_micro, by=c("bmid","subjid","arm","visit"))
misame <- full_join(misame, misame_bvit, by=c("bmid","subjid","arm","visit"))
misame <- full_join(misame, misame_HMO, by=c("bmid","subjid","arm","visit"))
misame <- full_join(misame, misame_protein, by=c("bmid","subjid","arm","visit"))
misame <- full_join(misame, misame_metabolomics, by=c("bmid","subjid","arm","visit"))
misame <- full_join(misame, misame_fsv, by=c("bmid","subjid","visit"))

#merge elicit
elicit <- full_join(elicit_macro, elicit_micro, by=c("bmid","subjid","arm","visit"))
elicit <- full_join(elicit, elicit_bvit, by=c("bmid","subjid","arm","visit"))
elicit <- full_join(elicit, elicit_HMO, by=c("bmid","subjid","arm","visit"))
elicit <- full_join(elicit, elicit_protein, by=c("bmid","subjid","arm","visit"))
elicit <- full_join(elicit, elicit_metabolomics, by=c("bmid","subjid","arm","visit"))
elicit <- full_join(elicit, elicit_fsv, by=c("bmid","subjid","visit"))

#merge vital
vital <- full_join(vital_macro, vital_micro, by=c("bmid","subjid","arm","visit"))
vital <- full_join(vital, vital_bvit, by=c("bmid","subjid","arm","visit"))
vital <- full_join(vital, vital_HMO, by=c("bmid","subjid","arm","visit"))
vital <- full_join(vital, vital_protein, by=c("bmid","subjid","arm","visit"))
vital <- full_join(vital, vital_metabolomics, by=c("bmid","subjid","arm","visit"))
vital <- full_join(vital, vital_fsv, by=c("bmid","subjid","visit"))


dim(misame)
dim(elicit)
dim(vital)

table(is.na(misame$visit))
table(is.na(elicit$visit))
table(is.na(vital$visit))

table(is.na(misame$arm))
table(is.na(elicit$arm))
table(is.na(vital$arm))

temp <- misame %>% filter(is.na(visit))

misame[is.na(misame$arm),]
misame$arm[misame$bmid=="pn34_107"] <- 1 #one bmid only

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
  mutate(milk_flag=1) %>% #use to mark presence of any data
  select(c("study","subjid", "arm", "visit", "bmid","milk_flag", everything())) %>%
    rename(armcd=arm)

#make sure all are numeric after merging
for(i in 6:ncol(d)){ d[,i] <- as.numeric(d[,i])}


#merge in baseline and time-varying covariates
dim(d)
dim(baseline)
dim(d_bmid)
head(d_bmid)
head(baseline)


#merge baseline and time varying and merge 
d_bmid %>% distinct(studyid, subjid, subjido) %>% dim()
baseline %>% distinct(studyid, subjid, subjido) %>% dim()
dim(d_bmid)
dim(baseline)
d_bmid <- left_join(d_bmid, baseline, by = c("studyid","subjid","subjido"))


#save bmids not in the milk data
d_bmid <- d_bmid %>% mutate(subjid=as.character(subjid), visit=as.numeric(visit), armcd=as.character(armcd))
d_unmerged_ids <- anti_join(d_bmid, d[,1:5], by = c("subjid", "armcd", "visit", "bmid"))
head(d_unmerged_ids)
dim(d_unmerged_ids)
table(d_unmerged_ids$studyid)
d_unmerged_ids$bmid
d_unmerged_ids$bmid %in% d$bmid

#create independent ID variable
d_bmid$clusterid <- d_bmid$siteid
d_bmid$clusterid[is.na(d_bmid$clusterid)] <- d_bmid$subjid[is.na(d_bmid$clusterid)]

#merge baseline and milk data
d <- d %>% subset(., select = -c(studyid,subjido))
d_bmid <- d_bmid %>% subset(., select = -c(agedays))
dim(d)
dim(d_bmid)
d <- left_join(d_bmid, d, by = c("subjid", "armcd", "visit", "bmid")) %>% filter(milk_flag==1) %>% subset(., select=-c(milk_flag))
dim(d)

table( d$arm, d$study)

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

table( d$arm, d$study)


table(d$arm)
table(d$armcd,d$arm)
colnames(d)

save(all_milk_components, elicit_milk_components, misame_milk_components, vital_milk_components, file=paste0(here::here(),"/metadata/milk_component.Rdata"))
# Normalize the fat-soluble vitamin A column to the lowercase name used by
# all_milk_components$micro. The Allen_FSV source files provide "vitamin.A"; the outcome
# list (and all downstream cleaning/labels) use "vitamin.a", so a case mismatch here breaks
# select(all_of(Yvars)) in the TMLE analyses. any_of() makes this a safe no-op if the column
# is already lowercase or absent.
d <- d %>% dplyr::rename_with(~ "vitamin.a", dplyr::any_of("vitamin.A"))

saveRDS(d, file=paste0(here::here(),"/data/merged_analysis_datasets.RDS"))

