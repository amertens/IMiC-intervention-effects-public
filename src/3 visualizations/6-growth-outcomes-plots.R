

#To make for Nolan:

#-Descriptive figures (density plots) overlaying child growth measures by time point



rm(list=ls())
source(paste0(here::here(),"/src/0-config.R"))

growth_res <- readRDS(paste0(here::here(),"/results/growth_intervention_effects_results.RDS"))

growth_res$res[[1]]$res %>% filter(measure=="ATE", biomarker=="growth_faltering_birth_6mo")
growth_res$res[[2]]$res %>% filter(measure=="ATE", biomarker=="growth_faltering_birth_6mo")
growth_res$res[[3]]$res %>% filter(measure=="ATE", biomarker=="growth_faltering_birth_6mo")

primary_outcomes <- c("haz_birth","haz_3mo", "haz_6mo", "whz_birth","whz_3mo",  "whz_6mo")
secondary_outcomes <- c("waz_6mo", "hcaz_6mo",
                        "len_vel_z_birth_6mo", 
                        "wt_vel_z_birth_6mo", 
                        "growth_faltering_birth_6mo")

plotdf <- bind_rows(growth_res$res[[1]]$res %>% filter(measure=="ATE") %>% mutate(study="Elicit"),
                    growth_res$res[[2]]$res %>% filter(measure=="ATE") %>% mutate(study="Misame"),
                    growth_res$res[[3]]$res %>% filter(measure=="ATE") %>% mutate(study="Vital") )
plotdf_primary <- plotdf %>% filter(biomarker %in% primary_outcomes)
head(plotdf_primary)

plotdf_primary <- plotdf_primary %>%
  mutate(sig= factor(sig),
         Measure=str_split_i(biomarker, "_", 1),
         Measure=str_to_upper(Measure),
         age=str_split_i(biomarker, "_", 2),
         age=str_to_title(age),
         age=factor(age, levels=c("Birth", "3mo", "6mo")),
         contrast=gsub("AZT", "Az.", contrast))

plotdf_primary$pval <- ci_to_pvalue(cil = plotdf_primary$cil,  ciu = plotdf_primary$ciu)
plotdf_primary$sig <- factor(ifelse(plotdf_primary$pval<0.05, "Yes", "No"))

plotdf_primary$Measure <- as.character(plotdf_primary$Measure)
plotdf_primary$Measure[plotdf_primary$Measure=="HAZ"] <- "LAZ"
plotdf_primary$Measure[plotdf_primary$Measure=="WHZ"] <- "WLZ"
#plotdf_primary$contrast[plotdf_primary$contrast=="Nico+Az."] <- "Nico+Az.\\*" #debug
plotdf_primary$contrast=factor(plotdf_primary$contrast, levels=rev(c("BEP", "IFA/BEP","BEP/IFA", "BEP/BEP","Az.", "Nico", "Nico+Az.","BEP+ExBf", "BEP+ExBf+Az." )), labels=rev(c("Postnatal BEP", "Postnatal BEP","Prenatal BEP", "Pre+Postnatal BEP","Az.", "Nico", "Nico+Az","Postnatal BEP", "Postnatal BEP+Az." )))
             


primary_growth_plot <- ggplot(plotdf_primary, aes(x=contrast, y=est, group=Measure, color=Measure, shape=sig, alpha=sig)) + 
  geom_point(position = position_dodge(width = 0.5), size=2) +
  geom_linerange(aes(ymin=cil, ymax=ciu), position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, linetype="dashed") +
  coord_flip() +
  scale_color_manual(values=tableau10) +
  facet_grid(study~age, scale="free_y") +
  scale_alpha_manual(values = c(0.7,1)) + guides(alpha = "none") +
  scale_color_manual(values = c(tableau10[c(3,2,1,5)])) +
  scale_shape_manual(values = c(1,19), name = "Significant") + 
  theme_bw() +
  theme(strip.background = element_blank(),
        axis.text = element_text(size = 6),
        legend.position="bottom")  + 
  xlab("Intervention arm\n(compared to Control)") + ylab("Z-score difference") 
primary_growth_plot

saveRDS(plotdf_primary, file=paste0(here::here(),"/figure-data/primary_growth_plot.RDS"))
ggsave(primary_growth_plot, filename=paste0(here::here(),"/figures/primary_growth_plot.png"), width=8, height=6, dpi=300)



presentation_plot <- ggplot(plotdf_primary %>% filter(age!="Birth"), aes(x=contrast, y=est, group=Measure, color=Measure, shape=sig, alpha=sig)) + 
  geom_point(position = position_dodge(width = 0.5), size=2) +
  geom_linerange(aes(ymin=cil, ymax=ciu), position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, linetype="dashed") +
  coord_flip() +
  scale_color_manual(values=tableau10) +
  facet_grid(study~age, scale="free_y") +
  scale_alpha_manual(values = c(0.7,1)) + guides(alpha = "none") +
  scale_color_manual(values = c(tableau10[c(3,2,1,5)])) +
  scale_shape_manual(values = c(1,19), name = "Significant") + 
  theme_bw() +
  theme(strip.background = element_blank(),
        axis.text = element_text(size = 6),
        legend.position="bottom")  + 
  xlab("Intervention arm\n(compared to Control)") + ylab("Z-score difference") 
presentation_plot
  theme_imic() +
  theme(strip.background = element_blank(),
        axis.text = element_text(size = 6),
        legend.position="bottom")  + 
  xlab("Intervention arm (compared to Control)") + ylab("Z-score difference") 
primary_growth_plot


ggsave(primary_growth_plot, filename=paste0(here::here(),"/figures/primary_growth_plot.png"), width=8, height=6, dpi=300)

#------------------------------------------------------------------------------
# Old plots
#------------------------------------------------------------------------------

# plotdf <- plotdf %>% mutate(biomarker = factor(biomarker, levels=c("haz_birth", "haz_3mo", "haz_6mo", "waz_birth", "waz_3mo", "waz_6mo", "waz_delta_3mo", "waz_delta_6mo", "whz_birth", "whz_3mo", "whz_delta_3mo", "whz_delta_6mo", "whz_6mo")))
# 
# plotdf_haz <- plotdf %>% filter(grepl("haz", biomarker)) 
# combined_growth_plot <- ggplot(growth_res$fit$Elicit$res %>% filter(measure=="ATE"), 
#                                aes(x=contrast, y=est)) + geom_point() +
#   geom_linerange(aes(ymin=cil, ymax=ciu)) +
#   geom_hline(yintercept = 0, linetype="dashed") +
#   coord_flip() +
#   facet_wrap(~biomarker, scale="free_y", ncol=3) +
#   theme(strip.background = element_blank(),
#         axis.text = element_text(size = 6),
#         legend.position="none")  + xlab("Intervention arm") + ylab("Z-score difference") + theme_bw()
# 
# 
# 
# #TO DO: make a plot function that allows for resorting the treatment and faceting variables
# elicit_growth_plot <- ggplot(growth_res$fit$Elicit$res %>% filter(measure=="ATE"), 
#                              aes(x=contrast, y=est)) + geom_point() +
#   geom_linerange(aes(ymin=cil, ymax=ciu)) +
#   geom_hline(yintercept = 0, linetype="dashed") +
#   coord_flip() +
#   facet_wrap(~biomarker, scale="free_y", ncol=3) +
#   theme(strip.background = element_blank(),
#         axis.text = element_text(size = 6),
#         legend.position="none")  + xlab("Intervention arm") + ylab("Z-score difference") + theme_bw()
# 
# misame_growth_plot <- ggplot(growth_res$fit$Misame$res %>% filter(measure=="ATE"), 
#                              aes(x=contrast, y=est)) + geom_point() +
#   geom_linerange(aes(ymin=cil, ymax=ciu)) +
#   geom_hline(yintercept = 0, linetype="dashed") +
#   coord_flip() +
#   facet_wrap(~biomarker, scale="free_y", ncol=3) +
#   theme(strip.background = element_blank(),
#         axis.text = element_text(size = 6),
#         legend.position="none")  + xlab("Intervention arm") + ylab("Z-score difference")+ theme_bw()
# 
# vital_growth_plot <- ggplot(growth_res$fit$Vital$res %>% filter(measure=="ATE"), 
#                             aes(x=contrast, y=est)) + geom_point() +
#   geom_linerange(aes(ymin=cil, ymax=ciu)) +
#   geom_hline(yintercept = 0, linetype="dashed") +
#   coord_flip() +
#   facet_wrap(~biomarker, scale="free_y", ncol=3) +
#   theme(strip.background = element_blank(),
#         axis.text = element_text(size = 6),
#         legend.position="none")  + xlab("Intervention arm") + ylab("Z-score difference")+ theme_bw()
# 
# saveRDS(list(growth_res=growth_res, 
#              elicit_growth_plot=elicit_growth_plot,
#              vital_growth_plot=vital_growth_plot,
#              misame_growth_plot=misame_growth_plot),
#         file=paste0(here::here(),"/results/growth_intervention_effects_results.RDS"))
