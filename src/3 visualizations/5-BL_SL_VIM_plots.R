# =============================================================================
# src/3 visualizations/5-BL_SL_VIM_plots.R
#
# Reads:  results/SL_BL_individual_lab_vim_res.RDS
#         results/SL_BL_vim_res.RDS
# Writes: figures/SL_BL_vim_lab_plot.png
#         figures/SLvim_lab_plot_elicit.jpeg
#         figures/SLvim_lab_plot_misame.jpeg
#         figures/SLvim_lab_plot_neg_control.jpeg
#         figures/SLvim_lab_plot_vital.jpeg
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
library(cowplot)


#vim <- readRDS(paste0(here::here(),"/results/SL_BL_vim_res.RDS"))
vim_ind <- readRDS(paste0(here::here(),"/results/SL_BL_individual_lab_vim_res.RDS"))
vim_lab <- readRDS(paste0(here::here(),"/results/SL_BL_vim_res.RDS"))



# vim_lab <- vim_lab %>% mutate(group=str_to_title(group),
#                               group=case_when(group=="Hmo" ~ "HMOs",
#                                           group=="Macro" ~ "Macronutrients",
#                                           group=="Micro" ~ "Micronutrients",
#                                           group=="Protein" ~ "Targeted proteins",
#                                           group=="Metabolomics" ~ "Targeted metabolomics",
#                                           group=="Bvit" ~ "B-vitamins",
#                                           group==group ~ group),
#                               group=factor(group, levels=c("All","Macronutrients","Micronutrients","Targeted proteins","B-vitamins","HMOs","Targeted metabolomics")))
# unique(vim_lab$group)

vim_lab <- vim_lab %>% filter(!grepl("_miss",group))
vim_ind <- vim_ind %>% filter(!grepl("_miss",group))

SLvim_lab_plot <- ggplot(vim_lab, aes(x=group, y=cvMSE, group=studyid , color=studyid )) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin=ci.lb, ymax=ci.ub), position = position_dodge(width = 0.5)) +
  coord_flip() +
  facet_wrap(~studyid) +
  scale_color_manual(values = tableau10, drop = FALSE) +
  labs(color = "Collection time") +
  theme(strip.background = element_blank(),
        axis.text = element_text(size = 6),
        legend.position="bottom")  +
  xlab("Group of\npredictors used") + ylab("CV-MSE") 
SLvim_lab_plot


SLvim_ind_plot <- ggplot(vim_ind, aes(x=group, y=cvMSE, group=studyid , color=studyid )) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin=ci.lb, ymax=ci.ub), position = position_dodge(width = 0.5)) +
  coord_flip() +
  facet_wrap(~studyid) +
  scale_color_manual(values = tableau10, drop = FALSE) +
  labs(color = "Collection time") +
  theme(strip.background = element_blank(),
        axis.text = element_text(size = 6),
        legend.position="bottom")  +
  xlab("Group of\npredictors used") + ylab("CV-MSE") 
SLvim_ind_plot

#------------------------------------------------------------------------------
# By group of predictors
#------------------------------------------------------------------------------
table(vim_lab$visit_f)
vim_lab$visit_f <- factor(vim_lab$visit_f, levels=c("<1 month","1-2 months","2-5 months"))

SLvim_lab_plot1 <- plot_SLvim(vim_lab %>% filter(studyid=="ELICIT"))
SLvim_lab_plot2 <- plot_SLvim(vim_lab %>% filter(studyid=="VITAL-Lactation"))
SLvim_lab_plot3 <- plot_SLvim(vim_lab %>% filter(studyid=="MISAME-3"), legend_pos="bottom")
SLvim_lab_plot <- plot_grid(SLvim_lab_plot1, SLvim_lab_plot2, SLvim_lab_plot3, ncol=1, rel_heights=c(1,1,1.2))
SLvim_lab_plot

ggsave(SLvim_lab_plot, file=paste0(here::here(),"/figures/SL_BL_vim_lab_plot.png"), width=9, height=7)


#-------------------------------------------------------------------------------
# Presentation plots
#-------------------------------------------------------------------------------

p_elicit <- plot_SLvim(vim_lab %>% filter(studyid=="ELICIT", arm!="blank"), legend_pos="bottom")
p_vital <- plot_SLvim(vim_lab %>% filter(studyid=="VITAL-Lactation", arm!="blank"), legend_pos="bottom")
p_misame <- plot_SLvim(vim_lab %>% filter(studyid=="MISAME-3", arm!="blank"), legend_pos="bottom")



SLvim_lab_plot_neg_control <- ggplot(vim_lab_neg_control, aes(x=group, y=cvAUC, group=visit_f, color=visit_f)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin=ci.lb, ymax=ci.ub), position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0.5, linetype="dashed") +
  coord_flip() +
  facet_grid(studyid~arm_f) +
  scale_color_manual(values=tableau10[-1]) +
  #theme_ki() +
  labs(color = "Collection time") +
  theme(strip.background = element_blank(),
        axis.text = element_text(size = 6),
        legend.position="bottom")  +
  xlab("Group of predictors added") + ylab("CV-AUC") # + 
  #ggtitle("Predicting Azithromycin arm\nfrom breastmilk samples in Elicit")

ggsave(p_elicit, file=paste0(here::here(),"/figures/SLvim_lab_plot_elicit.jpeg"), width=6, height=4)
ggsave(p_vital, file=paste0(here::here(),"/figures/SLvim_lab_plot_vital.jpeg"), width=6, height=4)
ggsave(p_misame, file=paste0(here::here(),"/figures/SLvim_lab_plot_misame.jpeg"), width=6, height=4)
ggsave(SLvim_lab_plot_neg_control, file=paste0(here::here(),"/figures/SLvim_lab_plot_neg_control.jpeg"), width=5, height=3)

