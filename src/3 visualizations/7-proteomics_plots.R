# =============================================================================
# src/3 visualizations/7-proteomics_plots.R
#
# Reads:  metadata/milk_component.Rdata
#         results/proteomics_intervention_effects_results.RDS
#         results/proteomics_intervention_effects_results_combined_arms.RDS
# Writes: figures/volcano_plot_proteomics.png
#         figures/volcano_plot_proteomics_arm_strat.png
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
d <- readRDS(paste0(here::here(),"/results/proteomics_intervention_effects_results_combined_arms.RDS"))

d$visit <- ""
names(d$res) <- paste0(d$study, "-", d$visit)

res_list=d
single_group=T
  
    res <- Map(extract_res, res_list$res)
    res <- rbindlist(res, idcol='studytime') %>% 
      as.data.frame() %>% filter(measure=="ATE")
    
    res$pval <- ci_to_pvalue(cil = res$cil, 
                             ciu = res$ciu)
    
    res = res %>% group_by(studytime, measure) %>% mutate(pval_adj=p.adjust(pval , method="BH") )
    #res$pval_adj <- p.adjust(res$pval , method="BH") 
    res$sig = 1*(res$pval < 0.05)
    res$sigFDR = 1*(res$pval_adj < 0.05)
    
  
  res <- res %>% mutate(
    biomarker= case_when(
      biomarker=="carbohydrate" ~ "Total Carbohydrate",
      biomarker==biomarker ~ biomarker))
  
  
  res <- clean_biomarker_labels(res)
  


#d$studytime <- 'Misame'
p<-plot_imic_volcano(res, labels=10)
p
ggsave(p, file=paste0(here::here(),"/figures/volcano_plot_proteomics.png"), width=18, height=10)


#-------------------------------------------------------------------------------
# Get arm-stratified results
#-------------------------------------------------------------------------------

d_strat <- readRDS(paste0(here::here(),"/results/proteomics_intervention_effects_results.RDS"))

d_strat$visit <- ""
names(d_strat$res) <- paste0(d_strat$study, "-", d_strat$visit)

res_list=d_strat
single_group=T

res <- Map(extract_res, res_list$res)
res <- rbindlist(res, idcol='studytime') %>% 
  as.data.frame()  %>% filter(measure=="ATE")

res$pval <- ci_to_pvalue(cil = res$cil, 
                         ciu = res$ciu)
#res$pval_adj <- p.adjust(res$pval , method="BH") 
res = res %>% group_by(studytime, measure) %>% mutate(pval_adj=p.adjust(pval , method="BH") )

res$sig = 1*(res$pval < 0.05)
res$sigFDR = 1*(res$pval_adj < 0.05)


res <- res %>% mutate(
  biomarker= case_when(
    biomarker=="carbohydrate" ~ "Total Carbohydrate",
    biomarker==biomarker ~ biomarker))


res <- clean_biomarker_labels(res)



#d$studytime <- 'Misame'
p_strat<-plot_imic_volcano(res, labels=10)
p_strat
ggsave(p, file=paste0(here::here(),"/figures/volcano_plot_proteomics_arm_strat.png"), width=18, height=10)


temp <- res %>% filter(studytime=="Vital-")
prop.table(table(temp$sigFDR)) * 100

# Re-extract from the non-combined-arms RDS (raw list form) and build a second
# volcano. `plot_imic_volcano` expects a tidy data.frame with pval_adj — extract
# from `d$res` first rather than passing the raw rowwise frame.
d_raw <- readRDS(paste0(here::here(),"/results/proteomics_intervention_effects_results.RDS"))
d_raw$visit <- ""
names(d_raw$res) <- paste0(d_raw$study, "-", d_raw$visit)
res2 <- Map(extract_res, d_raw$res) %>%
  rbindlist(idcol = "studytime") %>%
  as.data.frame() %>%
  filter(measure == "ATE") %>%
  mutate(pval = ci_to_pvalue(cil = cil, ciu = ciu)) %>%
  group_by(studytime, measure) %>%
  mutate(pval_adj = p.adjust(pval, method = "BH")) %>%
  ungroup() %>%
  mutate(sig = 1*(pval < 0.05),
         sigFDR = 1*(pval_adj < 0.05),
         outcome_group = "proteomics",
         label_f = biomarker)
p<-plot_imic_volcano(res2, labels=10)
ggsave(p, file=paste0(here::here(),"/figures/volcano_plot_proteomics.png"), width=18, height=10)




