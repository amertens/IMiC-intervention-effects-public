# figS6-microbiome-diversity.R
# =============================================================================
# Figure S6: intervention effects on human milk microbiome alpha diversity
# (observed richness and Shannon diversity), by study and visit.
#
# IN-REPO RECONSTRUCTION of a panel that the supplement previously embedded as an
# extracted image (Manuscript/qmd/extracted/media_supplement/media/image6.png),
# which had no generator under version control.
#
# Ported from src/3 visualizations/7-microbiome_plots.R.
#
# NOTE ON THE COLLISION: that script builds the diversity forest TWICE -- first
# the published variant (coloured by study, point shape by FDR significance), then
# a plain monochrome variant -- and BOTH write figures/forest_plot_microbiome_
# diversity.png. The monochrome one runs second, so the checked-in file is the
# variant that does NOT match the published figure. This script reproduces the
# published (coloured, shape-coded) variant under its own filename.
#
# Input:  results/microbiome_diversity_intervention_effects_results.RDS
# Output: figures/figureS6_microbiome_diversity.{pdf,eps,png}
#
# Run from repo root: Rscript figure-scripts/manuscript_figures/figS6-microbiome-diversity.R
# =============================================================================
source(paste0(here::here(), "/src/0-config.R"))

d_diversity <- readRDS(paste0(here::here(),
  "/results/microbiome_diversity_intervention_effects_results.RDS"))

res_diversity <- Map(extract_res, d_diversity$res) %>%
  rbindlist(idcol = "studytime") %>%
  mutate(outcome_group = "microbiome diversity") %>%
  as.data.frame()

# these results carry no p-value; derive it from the CI, then FDR-correct
res_diversity$pval <- ci_to_pvalue(cil = res_diversity$cil, ciu = res_diversity$ciu)
res_diversity <- res_diversity %>%
  group_by(outcome_group, measure) %>%
  mutate(pval_adj = p.adjust(pval, method = "BH")) %>%
  ungroup()

res_diversity <- res_diversity %>% mutate(
  studytime = as.character(studytime),
  studytime = case_when(
    studytime == "Misame-1" ~ "MISAME-III (14-21 days)",
    studytime == "Misame-2" ~ "MISAME-III (1-2 mo.)",
    studytime == "Misame-3" ~ "MISAME-III (3-4 mo.)",
    studytime == "Vital-40" ~ "Mumta-LW (1.5 mo.)",
    studytime == "Vital-56" ~ "Mumta-LW (2 mo.)",
    studytime == "Elicit-1" ~ "ELICIT (1 mo.)",
    studytime == "Elicit-5" ~ "ELICIT (5 mo.)"),
  study = case_when(
    grepl("MISAME", studytime) ~ "MISAME-III",
    grepl("Mumta",  studytime) ~ "Mumta-LW",
    grepl("ELICIT", studytime) ~ "ELICIT"),
  studytime = factor(studytime,
    levels = c("MISAME-III (14-21 days)", "MISAME-III (1-2 mo.)", "MISAME-III (3-4 mo.)",
               "Mumta-LW (1.5 mo.)", "Mumta-LW (2 mo.)", "ELICIT (1 mo.)", "ELICIT (5 mo.)")))

res_diversity <- res_diversity %>% filter(measure == "ATE")
res_diversity$biomarker <- gsub("_imp", " diversity", res_diversity$biomarker)
res_diversity$label_f   <- res_diversity$biomarker
res_diversity$studytime <- factor(res_diversity$studytime,
                                  levels = rev(unique(res_diversity$studytime)))

p <- ggplot(res_diversity, aes(x = studytime, y = est,
                               color = study, shape = factor(sigFDR))) +
  geom_point() +
  geom_linerange(aes(ymin = cil, ymax = ciu)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = tableau10[c(1, 3, 2)]) +
  coord_flip() +
  facet_wrap(~ label_f, scale = "free") +
  ggtitle("Intervention effect\n(BEP or Nicotinamide)") +
  # theme_imic (Helvetica, Reviewer-2 font floors) is already this session's global
  # default (theme_set() in 0_figure-functions.R) -- the explicit overrides removed
  # here had been defeating it (axis.text=6 below the 7pt floor; strip.background
  # blanked + strip.text=8, below the 10pt/bold/grey90 spec). panel.border re-added
  # to match the other harmonized figures (harmonized 2026-08-26).
  theme(panel.border    = element_rect(colour = "black", fill = NA, linewidth = 0.3),
        legend.position = "none") +
  xlab("Study and timepoint") + ylab("Z-score difference")

save_figure_3way(p, "figureS6_microbiome_diversity",
                 width = 9, height = 4,
                 dir   = paste0(here::here(), "/figures"))

message("Wrote figures/figureS6_microbiome_diversity.{pdf,eps,png}")
