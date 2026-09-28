# =============================================================================
# src/3 visualizations/5-SL_VIM_plots.R
#
# Formats the arm-classification AUC table from
# src/2 analysis/3b-SL_vim_individual_lab.R for plotting (visit and arm labels,
# predictor-group names, a blank Mumta-LW facet) and saves it as
# figure-data/SL_vim_plot_data.RDS, the input of
# figure-scripts/manuscript_figures/fig1-ml-vim-classifier.R (Fig 1A). It also
# draws working versions of the plot, per study and for the Elicit
# infant-azithromycin negative control; these are not printed exhibits.
#
# Inputs:  results/SL_individual_lab_vim_res.RDS
# Outputs: figure-data/SL_vim_plot_data.RDS
#          results/SLvim_lab_plots.RDS
#          figures/SLvim_lab_plot.png
#          figures/SLvim_lab_plot_{elicit,vital,misame,neg_control}.jpeg
# =============================================================================

rm(list=ls())
source(paste0(here::here(),"/src/0-config.R"))
library(cowplot)

vim_lab <- readRDS(paste0(here::here(),"/results/SL_individual_lab_vim_res.RDS"))

vim_lab <- vim_lab %>% mutate(visit=case_when(
  studyid=="ELICIT" & visit=="1"   ~ "1 month",
  studyid=="ELICIT" & visit=="5"   ~ "5 months",
  studyid=="VITAL-Lactation" & visit=="40"   ~ "1.5 months",
  studyid=="VITAL-Lactation" & visit=="56"   ~ "2 months",
  studyid=="MISAME-3" & visit=="1"   ~ "14-21 days",
  studyid=="MISAME-3" & visit=="2"   ~ "1-2 months",
  studyid=="MISAME-3" & visit=="3"   ~ "3-4 months"))
vim_lab <- vim_lab %>% mutate(visit_f=case_when(
  visit=="1 month" ~ "1-2 months",
  visit=="5 months"  ~ "2-5 months",
  visit=="1.5 months" ~ "1-2 months",
  visit== "2 months"  ~ "2-5 months",
  visit=="14-21 days" ~ "<1 month",
  visit== "1-2 months" ~ "1-2 months",
  visit=="3-4 months" ~ "2-5 months"),
  visit_f=factor(visit_f, levels=c("<1 month","1-2 months","2-5 months")))

unique(vim_lab$arm)
unique(vim_lab$visit)
unique(vim_lab$visit_f)
#"Az."          "Nico"         "Nico+Az."     "BEP/BEP"      "BEP/IFA"      "IFA/BEP"      "BEP+ExBf"     "BEP+ExBf+AZT" 

#add blank row to vital for blank facet
d <- data.frame(studyid="VITAL-Lactation",arm="blank",group="All",cvAUC=NA,ci.lb=0.5,ci.ub=0.5,visit="1.5 months",visit_f="1-2 months")
d2<-NULL
vim_lab <- bind_rows(vim_lab,d,d2)

vim_lab <- vim_lab %>% mutate(arm_f=case_when(
  arm=="blank" ~ "",
  arm=="Az." ~ "Infant Azithromycin",
  arm=="Nico" ~ "Nicotinamide",
  arm=="Nico+Az." ~ "Nicotinamide +\nInfant Azithromycin",
  arm=="BEP/BEP" ~ "Pre+Postnatal\nBEP",
  arm=="BEP/IFA" ~ "Prenatal BEP",
  arm=="IFA/BEP" ~ "Postnatal BEP",
  arm=="BEP+ExBf" ~ "Postnatal BEP",
  arm=="BEP+ExBf+AZT" ~ "Postnatal BEP +\nInfant Azithromycin"),
  arm_f=factor(arm_f, levels=c("","Infant Azithromycin","Prenatal BEP","Nicotinamide","Nicotinamide +\nInfant Azithromycin","Postnatal BEP","Pre+Postnatal\nBEP","Postnatal BEP +\nInfant Azithromycin")))
  
vim_lab$studyid[vim_lab$studyid=="VITAL-Lactation"] <- "Mumta-LW"
vim_lab$studyid[vim_lab$studyid=="MISAME-3"] <- "Misame-III"

table(vim_lab$studyid, vim_lab$arm_f)

vim_lab <- vim_lab %>% mutate(group=str_to_title(group),
                              group=case_when(group=="Hmo" ~ "HMOs",
                                          group=="Macro" ~ "Macronutrients",
                                          group=="Micro" ~ "Micronutrients",
                                          group=="Protein" ~ "Targeted proteins",
                                          group=="Metabolomics" ~ "Targeted metabolomics",
                                          group=="Bvit" ~ "B-vitamins",
                                          group==group ~ group),
                              group=factor(group, levels=rev(c("All","Macronutrients","Micronutrients","B-vitamins","HMOs","Targeted proteins","Targeted metabolomics"))))
unique(vim_lab$group)

#save plot data
saveRDS(vim_lab, file=paste0(here::here(),"/figure-data/SL_vim_plot_data.RDS"))
unique(vim_lab$group)

vim_lab_neg_control <- vim_lab %>% filter((studyid =="ELICIT" & arm == "Az."))
# The ELICIT infant-azithromycin arm stays in the plot data as a negative
# control: azithromycin was given to infants only after milk collection, so it
# cannot affect milk and the classifier should be at chance (cvAUC ~ 0.5).

#------------------------------------------------------------------------------
# By group of predictors
#------------------------------------------------------------------------------
table(vim_lab$visit_f)
vim_lab$visit_f <- factor(vim_lab$visit_f, levels=c("<1 month","1-2 months","2-5 months"))

SLvim_lab_plot1 <- plot_SLvim(vim_lab %>% filter(studyid=="ELICIT")) 
SLvim_lab_plot2 <- plot_SLvim(vim_lab %>% filter(studyid=="Mumta-LW"))
SLvim_lab_plot3 <- plot_SLvim(vim_lab %>% filter(studyid=="Misame-III"), legend_pos="bottom")
SLvim_lab_plot <- plot_grid(SLvim_lab_plot1, SLvim_lab_plot2, SLvim_lab_plot3, ncol=1, rel_heights=c(1,1,1.2))
SLvim_lab_plot

saveRDS(list(SLvim_lab_plot_elicit=SLvim_lab_plot1,
             SLvim_lab_plot_vital=SLvim_lab_plot2,
             SLvim_lab_plot_misame=SLvim_lab_plot3), 
        file=paste0(here::here(),"/results/SLvim_lab_plots.RDS"))
  
ggsave(SLvim_lab_plot, file=paste0(here::here(),"/figures/SLvim_lab_plot.png"), width=9, height=7)

#-------------------------------------------------------------------------------
# Per-study and negative-control plots (working figures, not printed exhibits)
#-------------------------------------------------------------------------------

p_elicit <- plot_SLvim(vim_lab %>% filter(studyid=="ELICIT", arm!="blank"), legend_pos="bottom")
p_vital <- plot_SLvim(vim_lab %>% filter(studyid=="Mumta-LW", arm!="blank"), legend_pos="bottom")
p_misame <- plot_SLvim(vim_lab %>% filter(studyid=="Misame-III", arm!="blank"), legend_pos="bottom")

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
