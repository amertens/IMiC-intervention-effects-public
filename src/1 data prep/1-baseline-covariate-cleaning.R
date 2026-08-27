
rm(list=ls())
source(paste0(here::here(),"/src/0-config.R"))

#-------------------------------------------------------------------------------
#IMIC harmonized datasets
#-------------------------------------------------------------------------------

misame <- read.csv(paste0(here::here(),"/data/imic_harmonized/MISAME_3_IMiC_analysis.csv")) %>% mutate(SUBJIDO=as.character(SUBJIDO))
vital <- read.csv(paste0(here::here(),"/data/imic_harmonized/VITAL_Lactation_IMiC_analysis.csv")) %>% mutate(SUBJIDO=as.character(SUBJIDO))
elicit <- read.csv(paste0(here::here(),"/data/imic_harmonized/ELICIT_IMiC_analysis.csv")) %>% mutate(SUBJIDO=as.character(SUBJIDO))
child <- read.csv(paste0(here::here(),"/data/imic_harmonized/CHILD_IMiC_analysis.csv")) %>% mutate(SUBJIDO=as.character(SUBJIDO))
length(unique(child$SUBJIDO))

#lower case variable names
colnames(misame) <- tolower(colnames(misame))
colnames(vital) <- tolower(colnames(vital))
colnames(elicit) <- tolower(colnames(elicit))
colnames(child) <- tolower(colnames(child))


#-------------------------------------------------------------------------------
# Extract household wealth and food security from raw data
#-------------------------------------------------------------------------------

#Load data
misame_raw <- haven::read_sas(paste0(here::here(),"/data/raw study data/misame/misame3_wide.sas7bdat"))
elicit_raw <- haven::read_sas(paste0(here::here(),"/data/raw study data/elicit/elicit_metadata_imic_042021.sas7bdat"))
vital_raw <- haven::read_sas(paste0(here::here(),"/data/raw study data/vital/crf3b_baseline.sas7bdat"))

#misame
misame_codebook <- makeVlist(misame_raw)
misame_codebook[grepl("food",misame_codebook$label),]
write.csv(misame_codebook, file=paste0(here::here(),"/metadata/misame_raw_codebook.csv"))

misame_raw <- misame_raw %>% 
  rename(subjido=idnew, hhwealth=incl_hhasset1comp, hhfoodsecure=hh_hfias_score, gwtgain=gwg_gwg) %>% 
  select(subjido, hhwealth, hhfoodsecure, gwtgain, gwg_rate) %>%
  mutate(subjido=as.character(subjido),
         hhwealth=as.numeric(scale(hhwealth)))

#gestational weight gain and rate- only in misame
summary(misame_raw$gwtgain)
summary(misame_raw$gwg_rate)

#merge with main data
misame <- left_join(misame, misame_raw, by="subjido")


#elicit
elicit_codebook <- makeVlist(elicit_raw)
write.csv(elicit_codebook, file=paste0(here::here(),"/metadata/elicit_raw_codebook.csv"))


elicit_asset <- elicit_raw %>% 
  select(PID, CEF_HAVE_MATTRESS:CEF_AVG_MTHLY_INCOME) %>% 
  rename(subjido=PID, INCTOT=CEF_AVG_MTHLY_INCOME )
elicit_raw=assetPCA(elicit_asset)

#merge with main data
elicit <- left_join(elicit, elicit_raw, by="subjido")


#vital
vital_codebook <- makeVlist(vital_raw)

colnames(vital_raw)
vital_asset <- vital_raw %>% 
  rename(subjido=crf3b_q1, nhh=crf3b_q7, nrooms=crf3b_q8, cookroom=crf3b_q9, fueltype=crf3b_q10, floor=crf3b_q11, wall=crf3b_q12) %>%
  select(subjido, nhh, nrooms, cookroom, fueltype, floor, wall) %>%
  mutate(nhh=as.numeric(nhh),
         nrooms=as.numeric(nrooms),
         n_per_room=nhh/nrooms,
         cookroom=factor(cookroom), 
         fueltype=factor(fueltype), 
         floor=factor(floor), 
         wall=factor(wall))

table(vital_asset$nhh)
table(vital_asset$nrooms)
table(vital_asset$cookroom)

ret=vital_asset
reorder=T
vital_raw=assetPCA(vital_asset, reorder=T)

#merge with main data
vital <- left_join(vital, vital_raw, by="subjido")


#Z-score income as measure of wealth
child$hhwealth <- as.numeric(scale(child$inctot))
child$hhwealth_quart <- factor(ntile(child$hhwealth, 4))

summary(misame$hhwealth)
summary(elicit$hhwealth)
summary(child$hhwealth)
summary(vital$hhwealth)

#-------------------------------------------------------------------------------
# Check that pre-pregnancy BMI is being used
#-------------------------------------------------------------------------------

head(misame)
misame <- misame %>% arrange(subjido, agedays)
table(misame$visit, !is.na(misame$mbmi))

misame %>% select(subjido, agedays, visit, mbmi)
#Misame is Baseline / Enrolment

head(vital)
vital <- vital %>% arrange(subjido, agedays)
table(vital$visit, !is.na(vital$mbmi))
table(vital$visit, vital$ageimpfl)

vital %>% filter(ageimpfl==1,visit=="Baseline") %>% select(subjido, agedays, visit, mbmi)

head(elicit)
elicit <- elicit %>% arrange(subjido, agedays)
table(elicit$visit, !is.na(elicit$mbmi))
table(elicit$visit, elicit$ageimpfl)

elicit %>% filter(visit=="Enrolment Visit") %>% select(subjido, agedays, visit, mbmi)
#elicit appears to be a birth visit (post-pregancy BMI)


head(child)
child <- child %>% arrange(subjido, agedays)
table(child$visit, !is.na(child$mbmi))
table(child$visit, child$ageimpfl)

#NOTE! Dropping children here, make sure none are dropped, use later visits if necessary
child %>% filter(visit=="Prenatal R00") %>% select(subjido, agedays, visit, mbmi)

#try just picking first measure:
temp <- child %>% group_by(subjido) %>%
  arrange(agedays) %>% slice(1)
table(temp$visit)
table(temp$visit, !is.na(temp$mbmi))

#-------------------------------------------------------------------------------
# Make baseline dataset of adjustment covariates
#-------------------------------------------------------------------------------

misame %>% filter(visit=="Baseline / Enrolment") %>% select(subjido) %>% distinct() %>% dim()
misame %>% filter(visit=="Baseline / Enrolment") %>% distinct() %>% dim()

vital %>% filter(visit=="Baseline", ageimpfl==0) %>% select(subjido) %>% distinct() %>% dim()
vital %>% filter(visit=="Baseline", ageimpfl==0) %>% distinct() %>% dim()

elicit %>% filter(visit=="Enrolment Visit") %>% select(subjido) %>% distinct() %>% dim()
elicit %>% filter(visit=="Enrolment Visit") %>% distinct() %>% dim()

child %>% filter(visit=="Prenatal R00")  %>% select(subjido) %>% distinct() %>% dim()
child %>% filter(visit=="Prenatal R00")   %>% distinct() %>% dim()

d <- bind_rows(misame %>% filter(visit=="Baseline / Enrolment"),
               vital %>% filter(visit=="Baseline", ageimpfl==0),
               elicit %>% filter(visit=="Enrolment Visit"),
               child%>% group_by(subjido) %>%arrange(agedays) %>% slice(1))

d %>% group_by(studyid) %>% summarise(n=n_distinct(subjid), mean(agedays), min(agedays), max(agedays)) #check number of subjects per study

saveRDS(d, file=paste0(here::here(),"/data/raw_imic_risk_factors.RDS"))


#-------------------------------------------------------------------------------
# clean covariates
#-------------------------------------------------------------------------------

#Improved drinking-water sources are defined as those that are likely to be protected from outside contamination, 
#and from faecal matter in particular. 
#Improved water sources include household connections, public standpipes, 
#boreholes, protected dug wells, protected springs and rainwater collection. 
#Unimproved water sources include unprotected wells, unprotected springs, surface water (e.g. river, dam or lake), vendor-provided water, 
#bottled water (unless water for other uses is available from an improved source) and tanker truck-provided water. 
#https://www.who.int/data/nutrition/nlis/info/improved-sanitation-facilities-and-drinking-water-sources#:~:text=Improved%20drinking%2Dwater%20sources%20are,protected%20springs%20and%20rainwater%20collection.

imp_level = c("Well dug: Protected pit","Pump well or borehole", "Tube well or borehole",
              "Protected well", "Tap in yard", "Public tap/public drinking fountain", "Improved",
              "Piped", "Community tap water", "Boring", "Public tap/stand pipe",  "Piped to yard/plot" )

unimp_level = c("Surface water(river/dam/lake/pond/stream", "Unprotected well", "Bottle","Water tanker" ,
                "Well dug: Non Protected pit","")

clean_water_levels <- function(d, all_imp=F){
  if(all_imp){
    d$imp_water_src <- "Improved"
  }else{
    d <- d %>% mutate(
      imp_water_src = case_when(h2osrcp %in% imp_level ~ "Improved",
                                h2osrcp %in% unimp_level ~ "Unimproved",
                                h2osrcp == h2osrcp ~ "Missing"))
  }
  d$imp_water_src[d$studyid=="CHILD"] <- "Improved"
  d$imp_water_src <- factor(d$imp_water_src, levels=c("Improved","Unimproved","Missing"))
  d <- droplevels(d)
  return(d)
}

d <- clean_water_levels(d)
table(d$studyid, d$imp_water_src)


d <- d %>% mutate(impfloor=case_when(floor %in% c("Cement","Linoleum (plastic)","Tiles","Tiles, parquet","Wood") ~ 1,
                                     floor %in% c("Clay","Natural/mud") ~ 0))
d$impfloor[d$studyid=="CHILD"] <- 1

table(d$studyid, d$impfloor)

#------------------------------------------------------------------------------
#

#Child sex, month of data collection, exclusive breastfeeding, 
#maternal age, maternal education, maternal BMI, parity, number of household members, 
#income, household size, construction, and WASH conditions (roof and wall materials, 
#cooking place, water source

#trenton: Maternal age, body mass index, height, hemoglobin level, gestational age, 
#mid-upper arm circumference and parity at inclusion

# month of data collection, exclusive breastfeeding, 
#maternal age, maternal education, maternal BMI, parity, number of household members, 
#income, household size, construction, and WASH conditions (roof and wall materials, 
#cooking place, water source



colnames(d)

baseline_covars <- c("sex", "gagebrth", "gagecm", "brthyr","birthwt", "birthlen", "mage",  "parity","primparity" ,  
                     "nperson","nrooms" , "meducyrs", "imp_water_src","impfloor",  "mbmi", "mhtcm", "mmuaccm" ,
                     "mhgb", "delivery","dlvloc",  "dvseason", "inctot","cookplac","hhwealth",  "hhfoodsecure",
                     "mmuaccm","mhtcm", "gwtgain")

d <- d %>% select("studyid", "siteid", "subjid","subjido",  "arm", "armcd", !!(baseline_covars))

#study specific covariates:
#Misame: siteid (make factor)
#Vital lactation: gagecm (method of gagebirth), mmuaccm as alternative to mbmi



#-------------------------------------------------------------------------------
# Save baseline data for risk factor analysis:
#-------------------------------------------------------------------------------

#check missingness
table(is.na(d))

#save all for risk factor analysis:
saveRDS(d, file=paste0(here::here(),"/data/clean_imic_risk_factors.RDS"))



#-------------------------------------------------------------------------------
# Impute missingness for intervention effects analysis
# Drop CHILD and save dataset
#-------------------------------------------------------------------------------

#save just intervention studies for intervention effects analysis
d <- d %>% filter(studyid!="CHILD")

# keep just adjustment and subgroup variables
colnames(d)
d <- d %>% subset(., select = c("studyid", "subjid" ,"subjido" , "arm","armcd","sex","mage", "meducyrs",  
                                "parity", "nperson","nrooms", "imp_water_src", "impfloor", "cookplac",
                                "mbmi", "mhtcm", "mmuaccm" ,
                                "dvseason", "dlvloc","hhwealth", "hhfoodsecure"))

head(d)
for(i in 6:ncol(d)){
  cat(colnames(d)[i],"\n")
  #print(table(d$studyid, is.na(d[,i])))
  print(prop.table(table(d$studyid, !is.na(d[,i])),1)*100)
  print(class((d[,i])))
}

#add missing category to categorical variables
d$cookplac[is.na(d$cookplac)] <- "missing"
d$impfloor[is.na(d$impfloor)] <- "missing"

#add missingness indicators and median impute continious
d$nperson_miss <- ifelse(is.na(d$nperson),1,0)
d$parity_miss <- ifelse(is.na(d$parity),1,0)
d$meducyrs_miss <- ifelse(is.na(d$meducyrs),1,0)
d$hhfoodsecure_miss <- ifelse(is.na(d$hhfoodsecure),1,0)
table(d$studyid, d$hhfoodsecure_miss)

d <- d %>% group_by(studyid) %>%
  mutate(nperson=case_when(is.na(nperson)~median(nperson, na.rm=T),nperson==nperson~nperson),
         parity=case_when(is.na(parity)~median(parity, na.rm=T),parity==parity~parity),
         hhfoodsecure=case_when(is.na(hhfoodsecure)~median(hhfoodsecure, na.rm=T),hhfoodsecure==hhfoodsecure~hhfoodsecure),
         meducyrs=case_when(is.na(meducyrs)~median(meducyrs, na.rm=T),meducyrs==meducyrs~meducyrs)) %>%
  ungroup()

#impute  ELICIT and VITAL-Lactation hh food security with 0 (will get dropped in the analysis)
table(d$studyid[is.na(d$hhfoodsecure)])
d$hhfoodsecure[is.na(d$hhfoodsecure)] <- 0

#make sure each variable is the correct class
glimpse(d)
d <- d %>% mutate(
  nrooms=as.numeric(nrooms),
  cookplac=factor(cookplac),
  dvseason=factor(dvseason),
  dlvloc=factor(dlvloc),
  hhfoodsecure=as.numeric(hhfoodsecure),
  cookplac=factor(cookplac),
  cookplac=factor(cookplac),
  sex=factor(sex),
  impfloor=factor(impfloor)
)
head(d)
table(is.na(d))

dim(d)
d %>% distinct(studyid, subjid, subjido) %>% dim()
d %>% distinct() %>% dim()


saveRDS(d, file=paste0(here::here(),"/data/clean_baseline_covariates.RDS"))
write.csv(d, file=paste0(here::here(),"/data/clean_baseline_covariates.csv"))




#-------------------------------------------------------------------------------
# Make a merge key dataset to merge milk with baseline
#-------------------------------------------------------------------------------

#subset to breastmilk observations
misame <- misame %>% filter(bmid != "")  %>%
  mutate(bmid=gsub("_lama","",bmid))
elicit <- elicit %>% filter(bmid != "")
vital <- vital %>% filter(bmid != "")

#collapse late visits
vital$visit[vital$visit=="Follow up Day 42.02"] <- "Follow up Day 42"
vital$visit[vital$visit=="Follow up Day 56.02"] <- "Follow up Day 56"

table(misame$visit)
table(elicit$visit)
table(vital$visit)

#recode visit to match the numeric visit from Nolans lab datasets
#misame: 1, 2, 3
misame <- misame %>% mutate(visit = case_when( 
  visit == "BM Collection 14-21D" ~ 1,
  visit == "BM Collection 1M-2M" ~ 2,
  visit == "BM Collection 3M-4M" ~ 3))
#elicit 1, 5
elicit <- elicit %>% mutate(visit = case_when( 
  visit == "1 Month Visit" ~ 1,
  visit == "5 Months Visit" ~ 5))
#vital 40, 56
vital <- vital %>% mutate(visit = case_when( 
  visit == "Follow up Day 42" ~ 40,
  visit == "Follow up Day 56" ~ 56))

d_timevar <- bind_rows(misame, elicit, vital)

#code time-varying EBF
d_timevar$ebf <- 1*(d_timevar$dur_ebf >= d_timevar$agedays)

d_timevar <- d_timevar %>% select("studyid", "subjid", "subjido", "bmid", "visit", "agedays", "ebf") %>% distinct()
# need to add milk collection date


saveRDS(d_timevar, file=paste0(here::here(),"/data/clean_timevar_covariates.RDS"))




