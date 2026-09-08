# figS4-triglyceride-means.R
# =============================================================================
# Figure S5: treatment-specific means (95% CI) of MISAME-III milk triglycerides,
# unadjusted and adjusted for total fat content.
#
# IN-REPO RECONSTRUCTION of a panel that the supplement previously embedded as an
# extracted image (Manuscript/qmd/extracted/media_supplement/media/image5.png),
# which had no generator under version control.
#
# Ported verbatim from figure-scripts/archive/fig_s4_triglyc_forest-plot_combined.R
# (debug/dput lines dropped, export switched to save_figure_3way). Layout is a
# 4 x 6 grid: intervention arm (rows) x visit-and-adjustment (columns), panels a-x.
#
# Inputs:  results/adjusted_intervention_effects_results_clean.RDS       (unadjusted)
#          results/fat_adjusted_metabolomics_intervention_effects_results.RDS
#          data/merged_analysis_datasets.RDS  (to drop sparse binary-recoded TGs)
# Output:  figures/figureS4_triglyceride_means.{pdf,eps,png}
#
# Run from repo root: Rscript figure-scripts/manuscript_figures/figS4-triglyceride-means.R
# =============================================================================
source(paste0(here::here(), "/src/0-config.R"))
load(file = paste0(here::here(), "/metadata/milk_component.Rdata"))

dfull <- readRDS(paste0(here::here(), "/data/merged_analysis_datasets.RDS"))

# triglycerides recoded to binary under sparse detection are excluded throughout
tgs <- dfull %>% filter(arm == "Control", study == "Misame", visit == "1") %>%
  select(all_milk_components$metabolomics) %>% select(starts_with("tg."))
cols_with_2_values <- names(tgs)[sapply(tgs, function(x) length(unique(x)) %in% c(1, 2))]

sigcat_of <- function(ate) {
  ate %>% mutate(
    sigcat = case_when(sigFDR == 1                ~ "Significant after FDR",
                       sigFDR == 0 & sig == 1     ~ "Significant before FDR",
                       sigFDR == 0 & sig == 0     ~ "Not Significant"),
    sigcat = factor(sigcat, levels = rev(c("Significant after FDR",
                                           "Significant before FDR",
                                           "Not Significant"))))
}

ARM_LEVELS <- c("Control", "BEP/IFA", "IFA/BEP", "BEP/BEP")
ARM_LABELS <- c("Control", "Prenatal BEP", "Postnatal BEP", "Pre+Postnatal BEP")

# ---- unadjusted --------------------------------------------------------------
res_full <- readRDS(paste0(here::here(), "/results/adjusted_intervention_effects_results_clean.RDS")) %>%
  filter(outcome_group == "tertiary", category == "Triglycerides")
res <- res_full %>% filter(measure == "MN")  %>% select(study, studytime, contrast, est, cil, ciu, biomarker)
ate <- res_full %>% filter(measure == "ATE") %>% sigcat_of() %>%
  select(study, studytime, contrast, biomarker, pval, pval_adj, sig, sigFDR, sigcat)

res <- left_join(res, ate, by = c("study", "studytime", "contrast", "biomarker")) %>%
  filter(!is.na(sigcat) | contrast == "Control") %>%
  filter(!(str_to_lower(biomarker) %in% cols_with_2_values))
res$sigcat[is.na(res$sigcat)] <- "Not Significant"   # control arm carries no ATE

res$contrast <- factor(res$contrast, levels = ARM_LEVELS)
res <- res %>% arrange(studytime, contrast, est)
res$ordered_biomarker <- factor(str_to_lower(res$biomarker),
                                levels = unique(str_to_lower(res$biomarker)))
res_misame <- res %>% filter(study == "Misame")
res_misame$contrast <- factor(res_misame$contrast, levels = ARM_LEVELS, labels = ARM_LABELS)

# ---- fat-adjusted ------------------------------------------------------------
res_full <- readRDS(paste0(here::here(), "/results/fat_adjusted_metabolomics_intervention_effects_results.RDS")) %>%
  filter(study == "Misame", category == "Triglycerides")
res <- res_full %>% filter(measure == "MN")  %>% select(studytime, contrast, est, cil, ciu, biomarker)
ate <- res_full %>% filter(measure == "ATE") %>% sigcat_of() %>%
  select(studytime, contrast, biomarker, pval, pval_adj, sig, sigFDR, sigcat)

res <- left_join(res, ate, by = c("studytime", "contrast", "biomarker")) %>%
  filter(!is.na(sigcat) | contrast == "Control") %>%
  filter(!(str_to_lower(biomarker) %in% cols_with_2_values))
res$sigcat[is.na(res$sigcat)] <- "Not Significant"

res$contrast <- factor(res$contrast, levels = ARM_LEVELS)
res      <- res      %>% arrange(studytime, contrast,  est)
res_full <- res_full %>% arrange(studytime, contrast, -est)
# order the x-axis on the fat-adjusted panels by the same biomarker sequence
res$ordered_biomarker <- factor(str_to_lower(res$biomarker),
                                levels = unique(str_to_lower(res_full$biomarker)))
res$contrast <- factor(res$contrast, levels = ARM_LEVELS, labels = ARM_LABELS)
res_fat <- res

# ---- combine and plot --------------------------------------------------------
res <- bind_rows(res_misame %>% mutate(adjustment = "Unadjusted"),
                 res_fat    %>% mutate(adjustment = "Fat-adjusted"))

# "MISAME-III" dropped from the strip text (every facet is MISAME-III, per the
# caption) and "Unadjusted"/"Fat-adjusted" shortened, so the now-larger/bold
# theme_imic() strip text (10pt floor, up from the previous unfloored 8pt) still
# fits this dense 6-column grid without clipping (harmonized 2026-08-26).
res$facet <- gsub("Unadjusted", "Unadj.",
                  gsub("Fat-adjusted", "Fat-adj.",
                       gsub("Misame\\s*", "", paste0(res$studytime, "\n", res$adjustment))))
res$facet <- factor(res$facet,
                    levels = c("(14-21 days)\nUnadj.", "(14-21 days)\nFat-adj.",
                               "(1-2 mo.)\nUnadj.",    "(1-2 mo.)\nFat-adj.",
                               "(3-4 mo.)\nUnadj.",    "(3-4 mo.)\nFat-adj."))

panel_labels <- expand.grid(facet    = sort(unique(res$facet)),
                            contrast = sort(unique(res$contrast)))
panel_labels$panel_num <- letters[1:nrow(panel_labels)]

p <- ggplot(res, aes(x = ordered_biomarker, y = est,
                     color = factor(sigcat), alpha = factor(sigcat))) +
  geom_linerange(aes(ymin = cil, ymax = ciu)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(contrast ~ facet) +
  geom_text(data = panel_labels, aes(label = panel_num),
            color = "black", alpha = 1,
            x = res$ordered_biomarker[1], y = 0.9,
            hjust = -1, size = 4, fontface = "bold") +
  scale_alpha_manual(values = c(0.3, 0.5, 0.8)) +
  scale_color_manual(values = c("grey50", tableau10)) +
  # theme_imic (Helvetica, Reviewer-2 font floors) harmonized 2026-08-26 (was bare
  # theme_bw() with axis.text=6, BELOW the 7pt floor, and a blanked/undersized strip).
  # panel.border/grid overrides kept: this is a dense 4x6 grid with an unlabelled
  # x-axis (too many biomarkers to show text), so the default theme_imic() grid would
  # just add clutter.
  theme_imic(base_size = 8) +
  theme(axis.text.x      = element_blank(),
        legend.position  = "bottom",
        legend.title     = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.3)) +
  ylab("Z-scored mean") + xlab("Milk component")

save_figure_3way(p, "figureS4_triglyceride_means",
                 width = 9, height = 7,
                 dir   = paste0(here::here(), "/figures"))

message("Wrote figures/figureS4_triglyceride_means.{pdf,eps,png}")
