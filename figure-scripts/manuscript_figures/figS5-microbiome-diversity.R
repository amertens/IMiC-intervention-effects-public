# =============================================================================
# figS5-microbiome-diversity.R
#
# Builds Fig S5: intervention effects (Z-score difference vs control, 95% CI) on
# human milk microbiome alpha diversity (observed richness and Shannon diversity), by
# study and visit, coloured by study and shaped by FDR significance. The results
# carry no p-value, so it is derived from the CI and BH-corrected here.
#
# Inputs:  results/microbiome_diversity_intervention_effects_results.RDS
#            (src/2 analysis/3_adjusted_analysis_microbiome.R)
# Outputs: figures/figureS5_microbiome_diversity.{pdf,eps,png}
# Run from repo root: Rscript figure-scripts/manuscript_figures/figS5-microbiome-diversity.R
# =============================================================================
source(paste0(here::here(), "/src/0-config.R"))
source(paste0(here::here(), "/figure-scripts/manuscript_figures/study_colors.R"))  # shared study colours

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
  # shared study colours, indexed by name. Study is also named on the axis.
  scale_color_manual(values = unname(imic_study_cols[c("ELICIT", "MISAME-III", "Mumta-LW")]),
                     guide = "none") +
  # key worded as in Fig 2; drop = FALSE keeps the FDR-significant (triangle) key
  scale_shape_manual(values = c("0" = 16, "1" = 17), labels = c("Not significant", "FDR-significant"),
                     name = "Statistical Significance", drop = FALSE) +
  coord_flip() +
  facet_wrap(~ label_f, scale = "free") +
  ggtitle("Intervention effect\n(BEP or Nicotinamide)") +
  # theme_imic() is the session default (theme_set() in 0_figure-functions.R); add a
  # panel border and plain (unbolded) titles, as in the main figures
  theme(panel.border    = element_rect(colour = "black", fill = NA, linewidth = 0.3),
        plot.title      = element_text(face = "plain", hjust = 0.5),
        strip.text      = element_text(face = "plain"),
        legend.position = "bottom") +
  xlab("Study and timepoint") + ylab("Z-score difference")

save_figure_3way(p, "figureS5_microbiome_diversity",
                 width = 9, height = 4,
                 dir   = paste0(here::here(), "/figures"))

message("Wrote figures/figureS5_microbiome_diversity.{pdf,eps,png}")
