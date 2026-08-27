# =============================================================================
# src/1 data prep/6-milq-adequecy.R
#
# Reads:  data/merged_analysis_datasets.RDS
#         data/MILQ standards/milq_age_specific_cutoffs_p05_p10_p90.csv
#         metadata/milk_component.Rdata
# Writes: data/milk_component_adequacy_MILQ.RDS
#         data/milq_age_specific_cutoffs_clean.RDS
#
# Paths above were recovered from this script's syntax tree and are
# repo-relative; they resolve from the repo root via here::here().
#
# Header generated from the code itself; it makes no claim about method.
# See README.md for run order and results/ARTIFACT_MANIFEST.csv for the
# exhibit each script feeds.
# =============================================================================


rm(list=ls())
source(paste0(here::here(),"/src/0-config.R"))

load(file=paste0(here::here(),"/metadata/milk_component.Rdata"))

milq <- read.csv(here("data/MILQ standards/milq_age_specific_cutoffs_p05_p10_p90.csv"))
unique(milq$nutrient)
milq[milq$nutrient=="Protein",]

knitr::kable(milq %>% filter(age_band=="3-4 m") %>% select(-age_band))

dfull <- readRDS(paste0(here::here(),"/data/merged_analysis_datasets.RDS"))
head(dfull)


unique(milq$nutrient)
colnames(dfull)
all_milk_components$macro
all_milk_components$micro
all_milk_components$bvit
#Find in IMiC data:
# "25(OH)D3"            "Vitamin D (ARA)"    

#rename to match imic component names
milq <- milq %>% mutate(
  nutrient_f=nutrient,
  nutrient=case_when(
    nutrient=="Protein"~"protein",
    nutrient=="Fat"~"fat",
    nutrient=="Carbohydrate"~"carbohydrate",
    nutrient=="Energy density"~"kcal.l",
    nutrient=="Sodium"~"na",
    nutrient=="Potassium"~"k",
    nutrient=="Magnesium"~"mg",
    nutrient=="Phosphorus"~"p",
    nutrient=="Calcium"~"ca",
    nutrient=="Copper"~"cu",
    nutrient=="Iron"~"fe",
    nutrient=="Zinc"~"zn",
    nutrient=="Selenium"~"se",
    nutrient=="Alpha-tocopherol"~"a.tocopherol",
    nutrient=="Gamma-tocopherol"~"g.tocopherol",
    nutrient=="Vitamin A (retinol)"~"vitamin.a",
    nutrient=="25(OH)D3"~"25(OH)D3",
    nutrient=="Vitamin B1"~"b1",
    nutrient=="Vitamin B2"~"b2",
    nutrient=="Vitamin B3"~"b3",
    nutrient=="Pantothenic acid"~"pa", 
    nutrient=="Vitamin B6"~"b6",
    nutrient=="Biotin"~"bio",
    nutrient=="Vitamin B12"~"b12",
    nutrient=="Choline"~"choline"
    ))

milq[milq$nutrient_f=="Protein",]

saveRDS(milq, file=paste0(here::here(),"/data/milq_age_specific_cutoffs_clean.RDS"))




imic_milq_nutrients <- c("na","k","mg","p","ca","cu","fe","zn","se","a.tocopherol","g.tocopherol","vitamin.a",
                         "b1","b2","b3","pa","b6","bio","b12","choline",
                         "protein","fat","carbohydrate","kcal.l")

colnames(dfull)



unique(milq$nutrient)
# Normalise vitamin A casing: the merged data column is `vitamin.A` (capital A) but
# imic_milq_nutrients + the MILQ rename (nutrient=="Vitamin A (retinol)" -> "vitamin.a")
# use lowercase. Without this, all_of("vitamin.a") errors and the pipeline can't run;
# and a naive fix that switched the select to `vitamin.A` would silently drop vitamin A
# at the `by="nutrient"` join (MILQ side is lowercase). Rename once so BOTH match.
dfull <- dfull %>% rename_with(~ "vitamin.a", any_of("vitamin.A"))
d <- dfull %>% select(study,studyid, arm, visit, sex, agedays, !!(Wvars), all_of(imic_milq_nutrients))
head(d)

milq[milq$nutrient=="fat",]
summary(d$fat)
milq[milq$nutrient=="protein",]
summary(d$protein)

#convert units of variables
d <- d %>%
  mutate(
    b12 = b12 * 1355 / 1e6,   # pmol/L  → µg/L
    b3  = b3 / 1000,  # µg/L → mg/L 
    cu  = cu / 1000,   # µg/L → mg/L 
    fe  = fe / 1000,   # µg/L → mg/L
    pa  = pa / 1000,   # µg/L → mg/L
    zn  = zn / 1000,   # µg/L → mg/L
    protein = protein * 10,
    fat = fat * 10,
    carbohydrate = carbohydrate * 10
  )

flag_out_of_range <- function(x, lo, hi) {
  ifelse(x < lo | x > hi, TRUE, FALSE)
}


d_milk <- d %>%
  mutate(
    # Micronutrients
    flag_b12_milq = flag_out_of_range(b12, 0.05, 1.5),
    flag_b3  = flag_out_of_range(b3, 0.3, 4.0),      # mg/L
    flag_pa  = flag_out_of_range(pa, 0.5, 5.0),      # mg/L
    flag_cu  = flag_out_of_range(cu, 0.05, 0.9),     # mg/L
    flag_fe  = flag_out_of_range(fe, 0.05, 1.0),     # mg/L
    flag_zn  = flag_out_of_range(zn, 0.2, 5.5),      # mg/L
    flag_se  = flag_out_of_range(se, 3, 25),         # µg/L
    
    # Macronutrients
    flag_protein = flag_out_of_range(protein, 4, 18),      # g/L
    flag_fat     = flag_out_of_range(fat, 10, 80),         # g/L
    flag_carb    = flag_out_of_range(carbohydrate, 50, 80) # g/L
  )

sapply(d_milk %>% select(starts_with("flag_")),
       function(x) mean(x, na.rm = TRUE)*100)
summary(d$b12)
summary(d$pa)
summary(d$b3)
table(d$b3[d$b3>10])

#check reasonable ranges after conversion to check correct conversion

#combine arms
table(d$arm)
d <- d %>% mutate(
  arm = case_when(
    arm=="Az." ~ "Control",
    arm=="BEP+ExBf+AZT" ~ "BEP",
    arm=="BEP+ExBf" ~ "BEP",
    arm=="Nico+Az." ~ "Nico",
    arm=="BEP/BEP" ~ "BEP",
    arm=="IFA/BEP" ~ "BEP",
    arm=="BEP/IFA" ~ "Control",
    arm==arm ~ arm
  )
)
d$arm <- factor(d$arm, levels=c("Control","BEP","Nico"))

d <- d %>% mutate(age_band=case_when(
  agedays < 4 ~ '0-3 d',
  agedays <= 17 ~ '4-17 d',
  agedays <= 31 ~ '18-31 d',
  agedays <= 60 ~ '1-2 m',
  agedays <= 90 ~ '2-3 m',
  agedays <= 120 ~ '3-4 m',
  agedays <= 150 ~ '4-5 m',
  agedays <= 180 ~ '5-6 m',
  TRUE ~ '5-6 m' # Beyond 6 months, RVs exist in tables but IMiC window is 1–6 mo; default to 5–6 m
))

head(d)


hm_component="fat"

milq_def_calc <- function(d, hm_component="na"){
  df <- d %>% select(study, studyid, arm, visit, sex, agedays, age_band,  !!(Wvars), !!(hm_component)) 
  colnames(df)[length(colnames(df))] <- "var" 
  df$nutrient <- hm_component
  
  df <- left_join(df, milq, by=c("age_band","nutrient"))
  head(df)
  
  df$def10 <- ifelse(df$var < df$P10, 1, 0)
  df$def5 <- ifelse(df$var < df$P05, 1, 0)
  return(df)
}



df = NULL
for(i in imic_milq_nutrients){
  temp <- milq_def_calc(d, hm_component=i)
  df = bind_rows(df, temp)
}
table(df$nutrient, df$def10)

#save milq deficiency data
head(df)
saveRDS(df, file=paste0(here::here(),"/data/milk_component_adequacy_MILQ.RDS"))
