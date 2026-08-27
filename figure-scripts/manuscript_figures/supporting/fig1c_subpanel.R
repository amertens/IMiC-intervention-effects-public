
rm(list=ls())
source(paste0(here::here(),"/src/0-config.R"))
# Source figure helpers (tableau10 palette, theme_imic(), save_figure_3way())
source(paste0(here::here(),"/figure-scripts/0_figure-functions.R"))
library(cowplot)

dfull <- readRDS(paste0(here::here(),"/data/merged_analysis_datasets.RDS"))
colnames(dfull)[1:30]
d <- dfull %>% select(studyid, study, subjid, visit, arm, agedays) %>%
  mutate(
    #standardize study names and visit times
    study=case_when(study=="Vital" ~ "Mumta",
                         study==study ~ study),
         
         
          studytime=paste0(study, "-", visit),
      studytime = case_when(
        studytime=="Misame-1" ~ "Misame (14-21 days)",
        studytime=="Misame-2" ~ "Misame (1-2 mo.)",
        studytime=="Misame-3" ~ "Misame (3-4 mo.)",
        studytime=="Mumta-40" ~ "Mumta (1.5 mo.)",
        studytime=="Mumta-56" ~ "Mumta (2 mo.)",
        studytime=="Elicit-1" ~ "Elicit (1 mo.)",
        studytime=="Elicit-5" ~ "Elicit (5 mo.)"
      ),
      studytime=factor(studytime, levels=rev(c("Elicit (1 mo.)", "Elicit (5 mo.)","Misame (14-21 days)", "Misame (1-2 mo.)", "Misame (3-4 mo.)",
                                           "Mumta (1.5 mo.)", "Mumta (2 mo.)")),
                       labels=rev(c("ELICIT (1 mo.)", "ELICIT (5 mo.)", "MISAME-III (14-21 days)", "MISAME-III (1-2 mo.)", "MISAME-III (3-4 mo.)",
                                      "Mumta-LW (1.5 mo.)", "Mumta-LW (2 mo.)"))),
         study=factor(study, levels=c( "Misame", "Mumta", "Elicit")))
head(d)

#make a ggplot violin plot with study_time on y-axis and agedays on x-axis, colored by study
p <- ggplot(d, aes(y=studytime, x=agedays/30.4167, fill=study, color=study)) +
  geom_jitter(width = 0, height = 0.15, alpha = 0.1) +
  geom_violin() +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  # scale_fill_manual(values=c( "#D55E00", "#CC79A7", "#F0E442", "#0072B2")) +
  # scale_color_manual(values=c( "#D55E00", "#CC79A7", "#F0E442", "#0072B2")) +
  scale_fill_manual(values=tableau10[c(4,2,1)]) +
  scale_color_manual(values=tableau10[c(4,2,1)]) +
  #geom_text() +
  labs(y="Study and Visit", x="Child Age in Days at Milk Sample Collection") +
  theme_bw() +
  theme(legend.position="none" ) 
p

ggsave(here::here("figures/fig1c.png"), p, width = 20, height = 10, units="cm", dpi = 600)




saveRDS(p, file= paste0(here::here(),"/figures/figure-data/figure1c_subplot.RDS"))




# ------------------------------------------------------------
# Recreate Fig 1C with formatting matched to simulated example
# ------------------------------------------------------------

p <- ggplot(d, aes(x = agedays/30.4167, y = studytime)) +
  geom_jitter(
    aes(color = study),
    height = 0.4,
    size   = 0.7,
    alpha  = 0.25) +
  geom_violin(
    fill      = NA,
    color     = "black",
    adjust    = 1,
    linewidth = 0.4) +
  geom_boxplot(
    fill      = NA,
    color     = "black",
    width     = 0.2,
    linewidth = 0.4,
    outlier.shape = NA) +
  geom_text(
    data = d %>%
      group_by(studytime) %>%
      summarize(
        n = n()),
    aes(
      x = max(d$agedays/30.4167, na.rm=T) + 1,
      y = studytime,
      label = paste0("n=", n)),
    size = 2.5,
    hjust = 0) +
  scale_color_manual(
    values = c(
      "Misame" = tableau10[4],
      "Mumta"  = tableau10[2],
      "Elicit" = tableau10[1])) +
  labs(
    x = "Child Age (Months) at Milk Sample Collection",
    y = "Study and Visit") +
  theme_bw(base_size = 7) +
  theme(
    axis.text        = element_text(size = 7),
    axis.title       = element_text(size = 7),
    plot.margin = margin(t = 10, r = 40, b = 10, l = 10),
    legend.position  = "none")

p

ggsave(here::here("figures/fig1c.png"), p, width = 20, height = 10, units="cm", dpi = 600)
# (Removed Google Drive export — file is committed to repo at figures/fig1c.png)


