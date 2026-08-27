# =============================================================================
# src/2 analysis/3_adjusted_analysis_pca_arm_strat.R
#
# Reads:  data/pca_analysis_datasets.RDS
#         metadata/milk_component.Rdata
#         results/pca_intervention_effects_results_arm_strat.RDS
# Writes: figures/pca_intervention_effects_arm_strat.png
#         results/pca_intervention_effects_results_arm_strat.RDS
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
d<-readRDS(paste0(here::here(),"/data/pca_analysis_datasets.RDS"))

head(d)

#Check for missingness in adjustment covariates.
missing_W <- d %>% select(all_of(Wvars)) %>% summarise_all(funs(sum(is.na(.))))
missing_W      

SL.lib  = c("SL.mean","SL.glm","SL.glmnet","SL.xgboost")
SL.lib  = c("SL.mean","SL.glm")


res <- d %>% group_by(study, visit) %>%
  do(res=run_bioTMLE(d=.,  Wvars = Wvars, bppar.debug=T, g_lib = SL.lib, Q_lib = SL.lib,
                     Yvars=c("macro_pca","micro_pca","bvit_pca","HMO_pca","protein_pca", "metabolomics_pca" ),
                     scale = TRUE,
                     bppar.type = BiocParallel::SnowParam()))
names(res$res) <- paste0(res$study, "-", res$visit)


res <- Map(extract_res, res$res) %>% rbindlist(., idcol='studytime') %>%
  mutate(label_f=biomarker ) %>%
  as.data.frame()
res

saveRDS(res, file=paste0(here::here(),"/results/pca_intervention_effects_results_arm_strat.RDS"))

res <- readRDS(file=paste0(here::here(),"/results/pca_intervention_effects_results_arm_strat.RDS"))


plotdf <- res 
plotdf <- plotdf %>% filter(measure=="ATE") %>%
  mutate(biomarker=str_to_title(gsub('_pca',"",biomarker)))

#split studytimeon "-"
plotdf <- plotdf %>% separate(studytime, into=c("study","time"), sep="-")
plotdf <- plotdf %>% mutate(study_intervention=paste0(study,"_", contrast))
head(plotdf)

#flip signs for consistent plotting
plotdf <- plotdf %>% mutate(cil=ifelse(est<0, -cil, cil),
                              ciu=ifelse(est<0, -ciu, ciu),
                              est=ifelse(est<0, -est, est))

p <-  ggplot(plotdf, aes(x=biomarker, y=est, color=label_f, group=time, shape=time)) + 
    geom_point(position = position_dodge(0.5)) +
    geom_linerange(aes(ymin=cil, ymax=ciu),position = position_dodge(0.5)) +
    geom_hline(yintercept = 0, linetype="dashed") +
    coord_flip() +
    facet_wrap(~study_intervention, scale="free") +
    ggtitle("Intervention Effects on the First Principal Component") +
    #scale_color_manual(values=tableau10) +
    theme(strip.background = element_blank(),
          axis.text = element_text(size = 6),
          strip.text = element_text(size = 8),
          legend.position.inside=c(0.1,0.8)) + 
    ylab("Difference") + xlab("Milk modality")
p

ggsave(p, file=paste0(here::here(),"/figures/pca_intervention_effects_arm_strat.png"), width=10, height=5)
