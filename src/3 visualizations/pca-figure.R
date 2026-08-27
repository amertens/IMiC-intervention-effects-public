# =============================================================================
# src/3 visualizations/pca-figure.R
#
# Reads:  results/pca_intervention_effects_results.RDS
# Writes: figure-data/pca_intervention_effects_results.RDS
#         figures/graphical_abstract_pca.png
#         figures/pca_intervention_effects.png
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

df <- readRDS(file=paste0(here::here(),"/results/pca_intervention_effects_results.RDS"))
saveRDS(df, file=paste0(here::here(),"/figure-data/pca_intervention_effects_results.RDS"))

p <-  ggplot(df, aes(x=label_f, y=est, color=label_f)) + geom_point() +
  geom_linerange(aes(ymin=cil, ymax=ciu)) +
  geom_hline(yintercept = 0, linetype="dashed") +
  coord_flip() +
  facet_wrap(studytime~contrast, scale="free") +
  ggtitle("Intervention Effects on the First Principal Component") +
  #scale_color_manual(values=tableau10) +
  theme(strip.background = element_blank(),
        axis.text = element_text(size = 6),
        strip.text = element_text(size = 8),
        legend.position="none") + 
  ylab("Average treatment effect") + xlab("Milk modality")
p
ggsave(p, file=paste0(here::here(),"/figures/pca_intervention_effects.png"), width=10, height=5)


#graphical abstract versions
head(df)
table(df$studytime)
plotdf <- df %>% filter(studytime %in% c("Elicit-1","Misame-1","Vital-40")) %>%
  mutate(studytime = recode(studytime,
                            "Elicit-1"="Elicit",
                            "Misame-1"="Misame",
                            "Vital-40"="Mumta-LW"))
p_ga <-  ggplot(plotdf, aes(x=label_f, y=est, color=label_f)) + 
  geom_point(size=3) +
  geom_linerange(aes(ymin=cil, ymax=ciu), size=1) +
  geom_hline(yintercept = 0, linetype="dashed") +
  coord_flip() +
  facet_grid(~studytime) +
  #ggtitle("Intervention Effects on\nthe First Principal Component") +
  #scale_color_manual(values=tableau10) +
  theme(strip.background = element_blank(),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 12),
        legend.position="none") + 
  ylab("Average treatment effect") + xlab("Milk modality")
p_ga
ggsave(p_ga, file=here("figures/graphical_abstract_pca.png"), width=4.5, height=2.1)
