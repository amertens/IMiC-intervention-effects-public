# =============================================================================
# src/2 analysis/milq-deficiency-reduction-analysis.R
#
# Effect of the interventions on the prevalence of milk nutrient concentrations
# below the MILQ age-specific 10th percentile (def10). For each study, visit and
# nutrient with at least 5 below-P10 samples, fits an adjusted TMLE (tmle
# package, SuperLearner library) of each pooled intervention arm vs Control and
# reports the risk ratio with its 95% CI; p-values are BH-adjusted within study
# and visit. Feeds Fig S2.
#
# Inputs:  data/milk_component_adequacy_MILQ.RDS (from 1 data prep/6-milq-adequecy.R)
# Outputs: results/milq_deficiency_reduction_analysis_results.RDS
# [needs restricted data]
# =============================================================================

rm(list=ls())
source(paste0(here::here(),"/src/0-config.R"))
library(tmle)

df = readRDS(paste0(here::here(),"/data/milk_component_adequacy_MILQ.RDS"))
table(df$nutrient)
knitr::kable(table(df$nutrient, df$def10))
table(df$nutrient, df$def5)

table(df$nutrient, df$unit)

df %>% group_by(nutrient) %>% summarize(mn=mean(var, na.rm=T), med=median(var, na.rm=T), min=min(var, na.rm=T), P05=mean(P05), P10=mean(P10), P90=mean(P90)) %>% 
  mutate(flag=1*(P90 < med)) %>% arrange(-flag) %>% as.data.frame()

# Drop study/visit/nutrient cells with fewer than 5 samples below P10
all_nutrients <- unique(df$nutrient)
df <- df %>% group_by(studyid, visit, nutrient) %>% filter(sum(def10, na.rm=T) >= 5) %>% ungroup()
all_nutrients[!(all_nutrients %in% unique(df$nutrient))]

# TMLE risk ratio of the binary outcome (Yvars) for each arm vs the first
# (Control) level, adjusted for the near-zero-variance-screened Wvars
run_TMLE <- function(d, Wvars, g_lib = c("SL.glm"),Q_lib = c("SL.glm"),
                     Yvars=c(all_milk_components$macro, all_milk_components$micro, all_milk_components$bvit),
                     scale=FALSE, cv.folds=1, adjust_biomarker=FALSE){
  
  cat(d$study[1],"\n")
  cat(d$visit[1],"\n")
  cat(d$nutrient[1],"\n")
  
  d <- droplevels(d)
  confounderData <- d %>% select(all_of(Wvars))
  biomarkerData <- d %>% select(all_of(Yvars))
  
  #drop NZV columns
  if(colnames(confounderData)[2]!="dummy"){ #skip if unadjusted analysis
    
    if(length(nearZeroVar(confounderData))>0){
      confounderData<-confounderData[,-nearZeroVar(confounderData)]
    }
    ARM <- confounderData[,1]
    confounderData <- cbind(ARM, design_matrix(as.data.frame(confounderData[,-1])))
  }
  if(length(nearZeroVar(biomarkerData))>0){
    biomarkerData<-biomarkerData[,-nearZeroVar(biomarkerData)]
  }
  
  fullres=NULL
  if(ncol(biomarkerData)>0){
    
    
    if(scale){
      biomarkerData = as.matrix(as.data.frame(scale_continuous(biomarkerData)))
    }
    
    d <- bind_cols(confounderData, biomarkerData)
    Yvar = colnames(biomarkerData)[1]
    d <- d[complete.cases(d),]
    
    
    for(i in 2:length(levels(d$arm))){
      res=NULL
      df <- d %>% filter(arm %in% c(levels(d$arm)[1], levels(d$arm)[i]))
      if(nrow(df)>0){
        Wvars <- df %>% select(-arm, -all_of(Yvar)) %>% as.data.frame()
        Avar <- ifelse(df$arm==levels(df$arm)[i],1,0)
        Y <- df[[Yvar]]
        tmle_out <- tmle(Y=Y,A=Avar,W=Wvars, family="binomial",  Q.SL.library = Q_lib, g.SL.library = g_lib, g.Delta.SL.library = g_lib,  verbose=F)
        
        res = data.frame(N=nrow(df),
                         nY=sum(Y, na.rm=T),
                         Y=Yvar,
                         contrast=levels(d$arm)[i],
                         RR=tmle_out$estimates$RR$psi,
                         ci.lb=tmle_out$estimates$RR$CI[1],
                         ci.ub=tmle_out$estimates$RR$CI[2],
                         pvalue=tmle_out$estimates$RR$pvalue)      
      }
      
      fullres=bind_rows(fullres, res)   
      
    }
    
  }
  return(fullres)
}

SL.lib  = c("SL.mean","SL.glm","SL.glmnet","SL.xgboost")

res <- df %>% group_by(study, visit, nutrient) %>%
  do(res=run_TMLE(d=.,  Wvars = Wvars,  g_lib = SL.lib, Q_lib = SL.lib,
                  Yvars="def10",
                  scale = FALSE))
res
names(res$res) <- paste0(res$study, "-", res$visit, "_", res$nutrient)

res_df <- NULL
for(i in 1:length(res$res)){
  if(!is.null(res$res[[i]])){
    if(nrow(res$res[[i]])>0){
      res_df <- bind_rows(res_df, data.frame(studytime_hm=names(res$res)[i], res$res[[i]]))
    }
  }
}

res_df$studytime <- str_split_i(res_df$studytime_hm, "_", 1)
res_df$study <- str_split_i(res_df$studytime, "-", 1)
res_df$visit <- str_split_i(res_df$studytime, "-", 2)
res_df$biomarker <- str_split_i(res_df$studytime_hm, "_", 2)
head(res_df)

res_df = res_df %>% group_by(study, visit) %>% mutate(pval_adj=p.adjust(pvalue , method="BH"))

res_df$sig = 1*(res_df$pvalue < 0.05)
res_df$sigFDR = 1*(res_df$pval_adj < 0.05)

#save results
saveRDS(res_df, file=paste0(here::here(),"/results/milq_deficiency_reduction_analysis_results.RDS"))
