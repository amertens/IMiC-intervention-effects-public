# figS2-growth-outcomes.R
# =============================================================================
# Figure S2: intervention effects on child growth (LAZ, WLZ) within the IMiC
# substudies of the three trials.
#
# IN-REPO RECONSTRUCTION of a panel that the supplement previously embedded as an
# extracted image (Manuscript/qmd/extracted/media_supplement/media/image1.png),
# which had no generator under version control.
#
# Ported verbatim from figure-scripts/archive/figure-s1-intervention-effects.R,
# which is the variant that matches the published figure: full arm labels and the
# Birth rows dropped for ELICIT and Mumta-LW (all their arms start after birth).
#
# NOTE ON THE ARCHIVE COLLISION: that script and src/3 visualizations/6-growth-
# outcomes-plots.R BOTH wrote figures/primary_growth_plot.png, and 6-growth-
# outcomes-plots.R does NOT apply the Birth filter or the full arm labels. Last
# writer won, so the checked-in primary_growth_plot.png is the wrong variant.
# This script writes its own distinct filename so neither can clobber it.
#
# Input:  results/growth_intervention_effects_results.RDS
# Output: figures/figureS2_growth_outcomes.{pdf,eps,png}
#
# Run from repo root: Rscript figure-scripts/manuscript_figures/figS2-growth-outcomes.R
# =============================================================================
source(paste0(here::here(), "/src/0-config.R"))   # brings in ci_to_pvalue(), theme/palette helpers

growth_res <- readRDS(paste0(here::here(), "/results/growth_intervention_effects_results.RDS"))

primary_outcomes <- c("haz_birth", "haz_3mo", "haz_6mo",
                      "whz_birth", "whz_3mo", "whz_6mo")

plotdf <- bind_rows(
  growth_res$res[[1]]$res %>% filter(measure == "ATE") %>% mutate(study = "Elicit"),
  growth_res$res[[2]]$res %>% filter(measure == "ATE") %>% mutate(study = "Misame"),
  growth_res$res[[3]]$res %>% filter(measure == "ATE") %>% mutate(study = "Vital"))

plotdf_primary <- plotdf %>%
  filter(biomarker %in% primary_outcomes) %>%
  mutate(sig      = factor(sig),
         Measure  = str_to_upper(str_split_i(biomarker, "_", 1)),
         age      = str_to_title(str_split_i(biomarker, "_", 2)),
         age      = factor(age, levels = c("Birth", "3mo", "6mo")),
         contrast = gsub("AZT", "Az.", contrast))

# significance from the CI (the results object carries no p-value for these)
plotdf_primary$pval <- ci_to_pvalue(cil = plotdf_primary$cil, ciu = plotdf_primary$ciu)
plotdf_primary$sig  <- factor(ifelse(plotdf_primary$pval < 0.05, "Yes", "No"))

plotdf_primary$study[plotdf_primary$study == "Vital"] <- "Mumta-LW"
plotdf_primary$study <- factor(plotdf_primary$study,
                               levels = c("Elicit", "Misame", "Mumta-LW"),
                               labels = c("ELICIT", "MISAME-III", "Mumta-LW"))

plotdf_primary$Measure <- as.character(plotdf_primary$Measure)
plotdf_primary$Measure[plotdf_primary$Measure == "HAZ"] <- "LAZ"
plotdf_primary$Measure[plotdf_primary$Measure == "WHZ"] <- "WLZ"

plotdf_primary$contrast <- gsub("Nico", "Nicotinamide",  plotdf_primary$contrast)
plotdf_primary$contrast <- gsub("Az.",  "Azithromycin",  plotdf_primary$contrast)
plotdf_primary$contrast <- gsub("\\+",  "\\+\n",         plotdf_primary$contrast)
plotdf_primary$contrast <- factor(
  plotdf_primary$contrast,
  levels = rev(c("IFA/BEP", "BEP/IFA", "BEP/BEP", "Azithromycin", "Nicotinamide",
                 "Nicotinamide+\nAzithromycin", "BEP+\nExBf", "BEP+\nExBf+\nAzithromycin")),
  labels = rev(c("Postnatal BEP", "Prenatal BEP", "Pre+Postnatal BEP", "Azithromycin",
                 "Nicotinamide", "Nicotinamide+\nAzithromycin", "Postnatal BEP",
                 "Postnatal BEP\nAzithromycin")))

# ELICIT and Mumta-LW randomized after birth, so their Birth rows are not estimable
plotdf_primary <- plotdf_primary %>%
  filter(!(study %in% c("ELICIT", "Mumta-LW") & age == "Birth"))

p <- ggplot(plotdf_primary,
            aes(x = contrast, y = est, group = Measure,
                color = Measure, shape = sig, alpha = sig)) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_linerange(aes(ymin = cil, ymax = ciu), position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  facet_grid(study ~ age, scale = "free_y") +
  scale_alpha_manual(values = c(0.7, 1)) + guides(alpha = "none") +
  scale_color_manual(values = tableau10[c(3, 2, 1, 5)]) +
  scale_shape_manual(values = c(1, 19), name = "Significant") +
  # theme_imic (Helvetica, Reviewer-2 font floors) harmonized 2026-08-26 (was bare
  # theme_bw() with an explicit axis.text=6, BELOW the 7pt floor). theme_imic's own
  # bold grey90 strip style (dormant everywhere else in this repo, since no other
  # theme_imic figure uses a real facet_grid with strips) now applies here for the
  # first time -- matches the "panel/strip >= 10pt, bold, gray90" spec verbatim.
  theme_imic(base_size = 8) +
  theme(panel.border    = element_rect(colour = "black", fill = NA, linewidth = 0.3),
        legend.position = "bottom") +
  xlab("Intervention arm\n(compared to Control)") + ylab("Z-score difference")

saveRDS(plotdf_primary, file = paste0(here::here(), "/figure-data/figureS2_growth_plot_data.RDS"))
save_figure_3way(p, "figureS2_growth_outcomes",
                 width  = science_dims$full_page$width,   # 7.25 in
                 height = 4.2,
                 dir    = paste0(here::here(), "/figures"))

message("Wrote figures/figureS2_growth_outcomes.{pdf,eps,png}")
