# =============================================================================
# 9-subgroup-analysis.R
#
# Effect modification by maternal body mass index (BMI <= 18.5 vs > 18.5): pooled-arm
# TMLE effects on the primary, secondary and tertiary milk outcomes estimated within
# each BMI stratum by study and visit, plus a z-test for the difference between the
# two stratum-specific effects. Supports the Results statement that maternal BMI did
# not consistently modify macronutrient or micronutrient effects (Section 3 of the
# online resource).
#
# Inputs : data/merged_analysis_datasets.RDS, metadata/milk_component.Rdata
# Outputs: results/adjusted_subgroup_results.RDS
# [needs restricted data]
# =============================================================================

source(paste0(here::here(),"/src/0-config.R"))



# Function to calc interaction p-value and clean results
compare_treatment_effects <- function(subgroup_df){
  
  # Calculate the z-statistic
  subgroup_df <- subgroup_df %>%   mutate(se=(est - cil)/1.96)
  
  z_stat <- (subgroup_df$est[1] - subgroup_df$est[2]) / sqrt(subgroup_df$se[1]^2 + subgroup_df$se[2]^2)
  
  # Calculate the p-value (two-tailed test)
  p_value <- 2 * (1 - pnorm(abs(z_stat)))
  
  subgroup_df$interaction_pvalue <- p_value 
  
  # Print results
  return(subgroup_df)
  
}

calc_interaction_pvalue <- function(res_list){
  
  res <- Map(extract_res, res_list$res)
  res <- rbindlist(res, idcol = "study_visit_subgroup") 
  res$study <- str_split(res$study_visit_subgroup,"-", simplify=T)[,1]
  res$subgroup <- str_split(res$study_visit_subgroup,"_", simplify=T)[,2]
  res$visit <- str_split(str_split(res$study_visit_subgroup,"-", simplify=T)[,2],"_", simplify=T)[,1]
  
  results <- res %>% filter(measure=="ATE") %>% 
    group_by(study, visit, biomarker) %>%
    do(compare_treatment_effects(.))
  return(results)
  
   
 }


load(file=paste0(here::here(),"/metadata/milk_component.Rdata"))
d<-readRDS(paste0(here::here(),"/data/merged_analysis_datasets.RDS"))

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
table(d$arm)
d$arm <- factor(d$arm, levels=c("Control","BEP","Nico"))
levels(d$arm)
table(d$study, d$arm)
table(is.na(d$arm))

#Check for missingness in adjustment covariates.
missing_W <- d %>% select(all_of(Wvars)) %>% summarise_all(funs(sum(is.na(.))))
missing_W    


SL.lib  = c("SL.glm")


summary(d$mbmi)
table(d$studyid,d$mbmi <= 18.5)
d$m_underweight <- ifelse(d$mbmi <= 18.5,"underweight","normal or overweight")

res_primary <- d %>% group_by(study, visit, m_underweight) %>%
  do(res=run_bioTMLE(d=.,  Wvars = Wvars, bppar.debug=T,
                     Yvars=c(all_milk_components$macro, all_milk_components$micro, all_milk_components$bvit),
                     scale = T))
names(res_primary$res) <- paste0(res_primary$study, "-", res_primary$visit, "_", res_primary$m_underweight)
res_primary<-calc_interaction_pvalue(res_primary)


res_secondary <- d %>% group_by(study, visit, m_underweight) %>%
  do(res=run_bioTMLE(d=.,  Wvars = Wvars,
                     Yvars=c(all_milk_components$hmo, all_milk_components$protein),
                     scale = TRUE))
names(res_secondary$res) <- paste0(res_secondary$study, "-", res_secondary$visit, "_", res_secondary$m_underweight)

res_secondary<-calc_interaction_pvalue(res_secondary)

res_tertiary <- d %>% group_by(study, visit, m_underweight) %>%
  do(res=run_bioTMLE(d=.,  Wvars = Wvars,
                     Yvars=all_milk_components$metabolomics,
                     scale = TRUE))
names(res_tertiary$res) <- paste0(res_tertiary$study, "-", res_tertiary$visit, "_", res_tertiary$m_underweight)
res_tertiary<-calc_interaction_pvalue(res_tertiary)

saveRDS(list(res_primary=res_primary, 
             res_secondary=res_secondary,
             res_tertiary=res_tertiary),
        file=paste0(here::here(),"/results/adjusted_subgroup_results.RDS"))
