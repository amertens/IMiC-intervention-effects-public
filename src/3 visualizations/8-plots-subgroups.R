# =============================================================================
# src/3 visualizations/8-plots-subgroups.R
#
# Reads:  results/adjusted_season_results.RDS
#         results/adjusted_subgroup_results.RDS
#         results/adjusted_subgroup_results_parity.RDS
# Writes: figure-data/subgroup_results.RDS
#         figures/season_plot.png
#         figures/subgroup_plot_primary.png
#         figures/subgroup_plot_tertiary.png
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

res <- readRDS(file=paste0(here::here(),"/results/adjusted_subgroup_results.RDS"))
res_season <- readRDS(file=paste0(here::here(),"/results/adjusted_season_results.RDS"))
res_parity <- readRDS(file=paste0(here::here(),"/results/adjusted_subgroup_results_parity.RDS"))

res <- bind_rows(
  res$res_primary %>% mutate(outcome_group="primary"),
  res$res_secondary %>% mutate(outcome_group="secondary"),
  res$res_tertiary %>% mutate(outcome_group="tertiary")
) %>% clean_biomarker_labels() %>% as.data.frame()

res_season <- bind_rows(
  res_season$res_primary %>% mutate(outcome_group="primary"),
  res_season$res_secondary %>% mutate(outcome_group="secondary"),
  res_season$res_tertiary %>% mutate(outcome_group="tertiary")
) %>% clean_biomarker_labels() %>% as.data.frame()

res_parity <- bind_rows(
  res_parity$res_primary %>% mutate(outcome_group="primary", visit=as.character(visit)),
  res_parity$res_secondary %>% mutate(outcome_group="secondary", visit=as.character(visit)),
  res_parity$res_tertiary %>% mutate(outcome_group="tertiary", visit=as.character(visit))
) %>% clean_biomarker_labels() %>% as.data.frame()

table(res$interaction_pvalue < 0.05)

#FDR correct the interaction Pvalues (remove duplicates)
head(res)
p_adj <- res %>% group_by(study, visit, biomarker, outcome_group) %>% 
  distinct(interaction_pvalue) %>% group_by(outcome_group) %>%
  mutate(interaction_pvalue_adj = p.adjust(interaction_pvalue, method="fdr")) %>% ungroup() %>%
  select(-interaction_pvalue) 
p_adj_season <- res_season %>% group_by(study, visit, biomarker, outcome_group) %>% 
  distinct(interaction_pvalue) %>% group_by(outcome_group) %>%
  mutate(interaction_pvalue_adj = p.adjust(interaction_pvalue, method="fdr")) %>% ungroup() %>%
  select(-interaction_pvalue)
p_adj_parity <- res_parity %>% group_by(study, visit, biomarker, outcome_group) %>% 
  distinct(interaction_pvalue) %>% group_by(outcome_group) %>%
  mutate(interaction_pvalue_adj = p.adjust(interaction_pvalue, method="fdr")) %>% ungroup() %>%
  select(-interaction_pvalue)

table(p_adj$interaction_pvalue_adj < 0.05)
table(p_adj$interaction_pvalue_adj < 0.2)
p_adj %>% filter(interaction_pvalue_adj < 0.05) 
p_adj %>% filter(interaction_pvalue_adj < 0.2) 

table(p_adj_season$interaction_pvalue_adj < 0.05)
table(p_adj_season$interaction_pvalue_adj < 0.2)
p_adj_season %>% filter(interaction_pvalue_adj < 0.05) 
p_adj_season %>% filter(interaction_pvalue_adj < 0.2) 


res <- left_join(res, p_adj, by=c("study","visit","biomarker","outcome_group"))
res_season <- left_join(res_season, p_adj_season, by=c("study","visit","biomarker","outcome_group"))
res_parity <- left_join(res_parity, p_adj_parity, by=c("study","visit","biomarker","outcome_group"))

plotdf <- res %>% filter(interaction_pvalue_adj < 0.2, measure=="ATE")


#save plot data
saveRDS(plotdf, file=paste0(here::here(),"/figure-data/subgroup_results.RDS"))


#maternal BMI plot
p <- ggplot(plotdf, aes(x=est, y=label_f , group = subgroup, color=subgroup)) + 
  geom_point(aes(shape=interaction_pvalue_adj < 0.05), position = position_dodge(width = 0.5)) + 
  geom_errorbarh(aes(xmin=cil, xmax=ciu), height=0.1, position = position_dodge(width = 0.5)) +
  geom_vline(xintercept=0.05, linetype="dashed") + 
  facet_grid(outcome_group~visit, scales="free_y") +
  theme_bw() + theme(legend.position="bottom") 


P_BMI_primary <- ggplot(plotdf %>% filter(outcome_group!="tertiary"), aes(x=est, y=label_f , group = subgroup, color=subgroup)) + 
  geom_point(aes(shape=interaction_pvalue_adj < 0.05), position = position_dodge(width = 0.5)) + 
  geom_errorbarh(aes(xmin=cil, xmax=ciu), height=0.1, position = position_dodge(width = 0.5)) +
  geom_vline(xintercept=0.05, linetype="dashed") + 
  facet_wrap(~paste0(study,"-",visit, " days"), scales="free_y") +
  xlab("ATE") + ylab("HM component") +
  guides(shape=F) +
  theme_bw() + theme(legend.position="bottom") 

plotdf_tertiary <- plotdf %>% filter(outcome_group=="tertiary", study=="Vital") %>%
                    arrange(visit, label_f) %>% 
                    mutate(label_f=factor(label_f, levels=unique(label_f)),
                           category=factor(category, levels=unique(category)))
P_BMI_tertiary <- ggplot(plotdf_tertiary, 
                         aes(x=est, y=label_f, group = subgroup, color=subgroup)) + 
  geom_point(aes(shape=interaction_pvalue_adj < 0.05), position = position_dodge(width = 0.5)) + 
  geom_errorbarh(aes(xmin=cil, xmax=ciu), height=0.1, position = position_dodge(width = 0.5)) +
  geom_vline(xintercept=0.05, linetype="dashed") + 
  #scale_y_discrete(labels=rev(plotdf_tertiary$category[plotdf$subgroup=="underweight"])) +
  facet_grid(~paste0(study,"-",visit, " days"), scales="free_y") +
  xlab("ATE") + ylab("HM component") +
  guides(shape=F) +
  theme_bw() + theme(legend.position="bottom", axis.text.y = element_text(size=8)) 



P_BMI_tertiary <- ggplot(plotdf_tertiary, 
                         aes(x=est, y=label_f, group = subgroup, color=subgroup)) + 
  geom_point(aes(shape=interaction_pvalue_adj < 0.05), position = position_dodge(width = 0.5)) + 
  geom_errorbarh(aes(xmin=cil, xmax=ciu), height=0.1, position = position_dodge(width = 0.5)) +
  geom_vline(xintercept=0.05, linetype="dashed") + 
  #scale_y_discrete(labels=rev(plotdf_tertiary$category[plotdf$subgroup=="underweight"])) +
  facet_wrap(~paste0(study,"-",visit, " days"), scales="free") +
  xlab("ATE") + ylab("HM component") +
  guides(shape=F) +
  theme_bw() + theme(legend.position="bottom", axis.text.y = element_text(size=8)) 


ggsave(P_BMI_primary, file=paste0(here::here(),"/figures/subgroup_plot_primary.png"), width=8, height=4)
ggsave(P_BMI_tertiary, file=paste0(here::here(),"/figures/subgroup_plot_tertiary.png"), width=7, height=9)

#misame plot
P_BMI_misame_primary <- ggplot(res %>% filter(measure=="ATE") %>% filter(study=="Misame", outcome_group=="primary"), aes(x=est, y=label_f , group = subgroup, color=subgroup)) + 
  geom_point(aes(shape=interaction_pvalue_adj < 0.05), position = position_dodge(width = 0.5)) + 
  geom_errorbarh(aes(xmin=cil, xmax=ciu), height=0.1, position = position_dodge(width = 0.5)) +
  geom_vline(xintercept=0.05, linetype="dashed") + 
  facet_wrap(~paste0(study,"-",visit, " days"), scales="free_y") +
  xlab("ATE") + ylab("HM component") +
  guides(shape=F) +
  theme_bw() + theme(legend.position="bottom") 

#season plot
season_plotdf <- res_season %>% filter(interaction_pvalue_adj < 0.2, measure=="ATE", outcome_group!="tertiary")
p_season <- ggplot(season_plotdf, 
                   aes(x=est, y=label_f , group = subgroup, color=subgroup)) + 
  geom_point(aes(shape=interaction_pvalue_adj < 0.05), position = position_dodge(width = 0.5)) + 
  geom_errorbarh(aes(xmin=cil, xmax=ciu), height=0.1, position = position_dodge(width = 0.5)) +
  geom_vline(xintercept=0.05, linetype="dashed") + 
  facet_wrap(~paste0(study,"-",visit, " month"), scales="free_y", ncol=2) +
  xlab("ATE") + ylab("HM component") +
  guides(shape=F) +
  theme_bw() + theme(legend.position="bottom") 

p_season

#save plot
ggsave(p_season, file=paste0(here::here(),"/figures/season_plot.png"), width=9, height=6)



table(p_adj_parity$interaction_pvalue_adj < 0.05)
table(p_adj_parity$interaction_pvalue_adj < 0.2)
p_adj_parity %>% filter(interaction_pvalue_adj < 0.05) 
p_adj_parity %>% filter(interaction_pvalue_adj < 0.2) 



parity_plotdf <- res_parity %>% 
  filter(interaction_pvalue_adj < 0.2, measure=="ATE", outcome_group!="tertiary",
         label_f!="Maternal Secretor Status")
p_parity <- ggplot(parity_plotdf, 
                   aes(x=est, y=label_f , group = subgroup, color=subgroup)) + 
  geom_point(aes(shape=interaction_pvalue_adj < 0.05), position = position_dodge(width = 0.5)) + 
  geom_errorbarh(aes(xmin=cil, xmax=ciu), height=0.1, position = position_dodge(width = 0.5)) +
  geom_vline(xintercept=0.05, linetype="dashed") + 
  facet_wrap(~paste0(study,"-",visit, " month"), scales="free_y", ncol=2) +
  xlab("ATE") + ylab("HM component") +
  guides(shape=F) +
  theme_bw() + theme(legend.position="bottom") 

p_parity
