# =============================================================================
# src/1 data prep/5-low-producers.R
#
# Reads:  data/clean milk data/CHILD/C_ICP.csv
#         data/merged_analysis_datasets.RDS
#         metadata/milk_component.Rdata
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
d<-readRDS(paste0(here::here(),"/data/merged_analysis_datasets.RDS"))

#add in CHILD
child_micro <- read.csv(paste0(here::here(),"/data/clean milk data/CHILD/C_ICP.csv")) %>% rename(BMID=Unnamed..0, SUBJID=pid) %>% rename_with(tolower) %>%
  mutate(studyid="CHILD", subjid=as.character(subjid)) %>% select(studyid, subjid, na)

d <- bind_rows(d, child_micro) %>% select(studyid, subjid, fat, kcal.l, na, mbmi, all_of(Wvars))

# based on the ICP data dictionary available on Synapse it shows Na is in mg/L; and that this is Na-23 which has a molar mass of 22.99 g/mol
# Unit conversion to mg/mmol = 22.99 mg/mmol
# mg/L divided by mg/mmol will get you the mmol/L units you want. 
#So just divide all values by 22.99 mg/mmol to get mmol/L for Na-23

summary(d$na/22.99 )

ggplot(d, aes(x=na, y=fat)) + geom_point() + scale_x_log10() + geom_smooth() + facet_wrap(~studyid) + theme_minimal()
ggplot(d, aes(x=na, y=kcal.l)) + geom_point() + scale_x_log10() + geom_smooth() + facet_wrap(~studyid) + theme_minimal()

ggplot(d %>% filter(!is.na(na)), aes(x=fat, group=(na/22.99 > 12), color=(na/22.99 > 12))) + geom_density() + facet_wrap(~studyid) + theme_minimal()
ggplot(d %>% filter(!is.na(na)), aes(x=fat, group=(na/22.99 > 16), color=(na/22.99 > 16))) + geom_density() + facet_wrap(~studyid) + theme_minimal()

ggplot(d %>% filter(!is.na(na)), aes(x=kcal.l, group=(na/22.99 > 12), color=(na/22.99 > 12))) + geom_density() + facet_wrap(~studyid) + theme_minimal()
ggplot(d %>% filter(!is.na(na)), aes(x=kcal.l, group=(na/22.99 > 16), color=(na/22.99 > 16))) + geom_density() + facet_wrap(~studyid) + theme_minimal()


d$na_mg.mmol<- d$na/22.99

d <- d %>% mutate(milk_production = case_when(
    d$na_mg.mmol < 12 ~ "Normal (<12 mmol/L Na)",
    d$na_mg.mmol >= 12 ~ "Low (>=12 mmol/L Na)",
    TRUE ~ NA
  ),
  milk_production=factor(milk_production, levels=c("Normal (<12 mmol/L Na)", "Low (>=12 mmol/L Na)")),
  milk_production_cat = case_when(
    d$na_mg.mmol < 12 ~ "Normal (<12 mmol/L Na)",
    d$na_mg.mmol > 16 ~ "Very low (>16 mmol/L Na)",
    d$na_mg.mmol >= 12 ~ "Low (12-16 mmol/L Na)",
    TRUE ~ NA
  ),
  milk_production_cat=factor(milk_production_cat, levels=c("Normal (<12 mmol/L Na)", "Low (12-16 mmol/L Na)", "Very low (>16 mmol/L Na)")),
  m_underweight=case_when(d$mbmi  < 18.5 ~ "Underweight", TRUE ~ "Normal")
  )
  
  head(d)
  table(d$milk_production)
  table(d$milk_production_cat)
  
  table(d$studyid, d$milk_production)
  table(d$studyid, d$milk_production_cat)

  prop.table(table(d$studyid, d$milk_production), margin=1)*100
  prop.table(table(d$studyid, d$milk_production_cat), margin=1)*100
  
  prop.table(table(d$m_underweight[d$studyid=="ELICIT"], d$milk_production[d$studyid=="ELICIT"]),  margin=1)*100
  prop.table(table(d$m_underweight[d$studyid=="MISAME-3"], d$milk_production[d$studyid=="MISAME-3"]),  margin=1)*100
  prop.table(table(d$m_underweight[d$studyid=="VITAL-Lactation"], d$milk_production[d$studyid=="VITAL-Lactation"]),  margin=1)*100

  #effect of intervention on milk production
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
  
  table(d$arm, d$milk_production, d$studyid)
  
  
  Wvars
  library(tmle)
  SL.library = c("SL.glm", "tmle.SL.dbarts2", "SL.glmnet")
  SL.library = c("SL.glm")
  
  #misame tmle
  d_misame <- d[d$studyid=="MISAME-3",] %>% select(milk_production, arm, all_of(Wvars)) %>% droplevels()
  d_misame<-d_misame[,-nzv(d_misame)]
  misame_res <- tmle(Y=ifelse(d_misame$milk_production=="Low (>=12 mmol/L Na)", 1, 0),  A=ifelse(d_misame$arm=="Control", 0, 1), W=d_misame[,-c(1:2)],  family="binomial", Q.SL.library=SL.library, g.SL.library=SL.library, g.Delta.SL.library=SL.library) 
  misame_res=summary(misame_res)
  misame_res= data.frame(studyid="MISAME-3",RR=misame_res$estimates$RR$psi, ci.lb=misame_res$estimates$RR$CI[1], ci.ub=misame_res$estimates$RR$CI[2], pval=misame_res$estimates$RR$pvalue)

  #elicit tmle
  d_elicit <- d[d$studyid=="ELICIT",] %>% select(milk_production, arm, all_of(Wvars)) %>% droplevels()
  d_elicit<-d_elicit[,-nzv(d_elicit)]
  elicit_res <- tmle(Y=ifelse(d_elicit$milk_production=="Low (>=12 mmol/L Na)", 1, 0),  A=ifelse(d_elicit$arm=="Control", 0, 1), W=d_elicit[,-c(1:2)],  family="binomial", Q.SL.library=SL.library, g.SL.library=SL.library, g.Delta.SL.library=SL.library)
  elicit_res=summary(elicit_res)
  elicit_res= data.frame(studyid="ELICIT", RR=elicit_res$estimates$RR$psi, ci.lb=elicit_res$estimates$RR$CI[1], ci.ub=elicit_res$estimates$RR$CI[2], pval=elicit_res$estimates$RR$pvalue)
  
  #vital tmle
  d_vital <- d[d$studyid=="VITAL-Lactation",] %>% select(milk_production, arm, all_of(Wvars)) %>% droplevels()
  d_vital<-d_vital[,-nzv(d_vital)]
  vital_res <- tmle(Y=ifelse(d_vital$milk_production=="Low (>=12 mmol/L Na)", 1, 0),  A=ifelse(d_vital$arm=="Control", 0, 1), W=d_vital[,-c(1:2)],  family="binomial", Q.SL.library=SL.library, g.SL.library=SL.library, g.Delta.SL.library=SL.library)
  vital_res=summary(vital_res)
  vital_res= data.frame(studyid="VITAL",RR=vital_res$estimates$RR$psi, ci.lb=vital_res$estimates$RR$CI[1], ci.ub=vital_res$estimates$RR$CI[2], pval=vital_res$estimates$RR$pvalue)
  
  tmle_res = bind_rows(misame_res, elicit_res, vital_res)
  
  ggplot(tmle_res, aes(x=studyid, y=RR)) + geom_point() + 
    geom_errorbar(aes(ymin=ci.lb, ymax=ci.ub)) + 
    geom_hline(yintercept=1, linetype="dashed") + 
    #log transformation
    scale_y_log10() +
    coord_flip() +
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(y="Risk Ratio", x="Study ID") + geom_text(aes(label=ifelse(pval<0.05, "*", "")), y=1.5)
  
  
  
#look at 
