# =============================================================================
# src/3 visualizations/7-microbiome_plots_arm_strat.R
#
# Reads:  metadata/milk_component.Rdata
#         results/microbiome_diversity_intervention_effects_results_arm_strat.RDS
#         results/microbiome_intervention_effects_results_arm_strat.RDS
# Writes: figures/forest_plot_microbiome_diversity_arm_strat.png
#         figures/volcano_plot_microbiome_arm_strat.png
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
d <- readRDS(paste0(here::here(),"/results/microbiome_intervention_effects_results_arm_strat.RDS"))
d_diversity <- readRDS(paste0(here::here(),"/results/microbiome_diversity_intervention_effects_results_arm_strat.RDS"))
head(d)



res <- Map(extract_res, d$res) %>% rbindlist(., idcol='studytime') %>% 
  mutate(outcome_group='microbiome', label_f=biomarker ) %>% 
  as.data.frame()
res_diversity <- Map(extract_res, d_diversity$res) %>% rbindlist(., idcol='studytime') %>% 
  mutate(outcome_group='microbiome diversity') %>% 
  as.data.frame()

res$pval <- ci_to_pvalue(cil = res$cil, ciu = res$ciu)
res = res %>% group_by(outcome_group, measure) %>% mutate(pval_adj=p.adjust(pval , method="BH") )
res_diversity$pval <- ci_to_pvalue(cil = res_diversity$cil, ciu = res_diversity$ciu)
res_diversity = res_diversity %>% group_by(outcome_group, measure) %>% mutate(pval_adj=p.adjust(pval , method="BH") )


#res <- extract_bioTMLE_results(d)
res <- res %>% mutate(
  studytime=as.character(studytime),
  studytime = case_when(
    studytime=="Misame-1" ~ "Misame (14-21 days)",
    studytime=="Misame-2" ~ "Misame (1-2 mo.)",
    studytime=="Misame-3" ~ "Misame (3-4 mo.)",
    studytime=="Vital-40" ~ "Vital (1.5 mo.)",
    studytime=="Vital-56" ~ "Vital (2 mo.)",
    studytime=="Elicit-1" ~ "Elicit (1 mo.)",
    studytime=="Elicit-5" ~ "Elicit (5 mo.)"
  ),
  studytime=factor(studytime, levels=c("Misame (14-21 days)", "Misame (1-2 mo.)", "Misame (3-4 mo.)",
                                       "Vital (1.5 mo.)", "Vital (2 mo.)", "Elicit (1 mo.)", "Elicit (5 mo.)"))
)
res_diversity <- res_diversity %>% mutate(
  studytime=as.character(studytime),
  studytime = case_when(
    studytime=="Misame-1" ~ "Misame (14-21 days)",
    studytime=="Misame-2" ~ "Misame (1-2 mo.)",
    studytime=="Misame-3" ~ "Misame (3-4 mo.)",
    studytime=="Vital-40" ~ "Vital (1.5 mo.)",
    studytime=="Vital-56" ~ "Vital (2 mo.)",
    studytime=="Elicit-1" ~ "Elicit (1 mo.)",
    studytime=="Elicit-5" ~ "Elicit (5 mo.)"
  ),
  studytime=factor(studytime, levels=c("Misame (14-21 days)", "Misame (1-2 mo.)", "Misame (3-4 mo.)",
                                       "Vital (1.5 mo.)", "Vital (2 mo.)", "Elicit (1 mo.)", "Elicit (5 mo.)"))
)

res <- res %>% filter(measure=="ATE") %>%
  group_by(studytime, biomarker) %>% filter(abs(est)==max(abs(est))) %>% ungroup() %>% arrange(chi_pval_adj)
res_diversity <- res_diversity %>% filter(measure=="ATE") 



p_microbiome_plot1 <- plot_imic_volcano(res %>% filter(grepl("Elicit",studytime)))
p_microbiome_plot2 <- plot_imic_volcano(res %>% filter(grepl("Vital",studytime)))
p_microbiome_plot3 <- plot_imic_volcano(res %>% filter(grepl("Misame",studytime))) 
p_microbiome_plot <- cowplot::plot_grid(p_microbiome_plot3,p_microbiome_plot2,p_microbiome_plot1, ncol=1, rel_heights=c(1,1,1))
ggsave(p_microbiome_plot, file=paste0(here::here(),"/figures/volcano_plot_microbiome_arm_strat.png"), width=18, height=10)
#NOTE! Significant results exist if I keep stratified by treatment arm
#cross-check if the chi-pval matches the contrast-specific pval when there are only 2 arms
table(res$chi_pval < 0.05, 1*((res$cil>0 & res$ciu>0)|(res$cil<0 & res$ciu<0)))

res_diversity$biomarker <- gsub("_imp"," diversity",res_diversity$biomarker)
res_diversity$label_f <- res_diversity$biomarker
res_diversity$studytime <- factor(res_diversity$studytime, levels = rev(unique(res_diversity$studytime)))
p_diversity <-   ggplot(res_diversity, aes(x=studytime, y=est, color=contrast, shape=factor(sigFDR) )) + geom_point(position = position_dodge(0.5)) +
    geom_linerange(aes(ymin=cil, ymax=ciu), position = position_dodge(0.5)) +
    geom_hline(yintercept = 0, linetype="dashed") +
    coord_flip() +
    facet_wrap(~label_f, scale="free") +
    ggtitle("Intervention effect\n(BEP or Nico.)") +
    theme(strip.background = element_blank(),
          axis.text = element_text(size = 6),
          strip.text = element_text(size = 8),
          legend.position="bottom") + xlab("Study and timepoint") + ylab("Z-score difference")

p_diversity 




ggsave(p_diversity, file=paste0(here::here(),"/figures/forest_plot_microbiome_diversity_arm_strat.png"), width=9, height=6)
