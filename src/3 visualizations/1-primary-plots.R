# =============================================================================
# src/3 visualizations/1-primary-plots.R
#
# Reads:  metadata/milk_component.Rdata
#         results/adjusted_intervention_effects_results_clean.RDS
# Writes: results/plotdf_misame_macro.csv
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

res <- readRDS(file=paste0(here::here(),"/results/adjusted_intervention_effects_results_clean.RDS"))
head(res)

plotdf <- res %>% filter(study=="Misame", str_to_lower(biomarker) %in% misame_milk_components$macro, measure=="ATE") %>%
  mutate(visit = factor(visit, levels=c("14-21 days", "1-2 mo.","3-4 mo.")),
         contrast = factor(contrast, levels=c("IFA/BEP", "BEP/IFA", "BEP/BEP")))
write.csv(plotdf, file=paste0(here::here(),"/results/plotdf_misame_macro.csv"), row.names=F)


plotdf2 <- res %>% filter(study=="Misame", str_to_lower(biomarker) %in% misame_milk_components$macro, measure=="MN") %>%
  mutate(visit = factor(visit, levels=c("14-21 days", "1-2 mo.","3-4 mo.")),
         contrast = factor(contrast, levels=c("Control","IFA/BEP", "BEP/IFA", "BEP/BEP")))
ggplot(plotdf2, aes(x=contrast, y=est)) + geom_point() +
  geom_linerange(aes(ymin=cil, ymax=ciu)) +
  geom_hline(yintercept = 0, linetype="dashed") +
  facet_wrap(visit~biomarker, scale="fixed") +
  ggtitle("Scaled means by intervention arm") +
  theme(strip.background = element_blank(),
        axis.text = element_text(size = 6),
        strip.text = element_text(size = 8),
        legend.position="none")  + ylab("Intervention arm") + 
  xlab("Z-score difference")
