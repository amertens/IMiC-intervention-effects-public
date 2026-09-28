# =============================================================================
# figS1-growth-outcomes.R
#
# Builds Fig S1: intervention effects on child growth (LAZ, WLZ; Z-score difference
# vs control, 95% CI) at birth, 3 and 6 months within the IMiC substudies of the
# three trials, one row per study and arm. Birth rows are dropped for ELICIT and
# Mumta-LW, whose arms all start after birth. Significance is derived from the CI.
#
# Inputs:  results/growth_intervention_effects_results.RDS (src/2 analysis/5-growth-outcomes.R)
# Outputs: figures/figureS1_growth_outcomes.{pdf,eps,png}
# Run from repo root: Rscript figure-scripts/manuscript_figures/figS1-growth-outcomes.R
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
         age      = factor(age, levels = c("Birth", "3mo", "6mo"),
                           labels = c("Birth", "3 mo.", "6 mo.")),   # "mo." as in the main figures
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
  # spaced "+" as in Fig 1 (and the "+" the Mumta-LW BEP+azithromycin label was missing)
  labels = rev(c("Postnatal BEP", "Prenatal BEP", "Pre + Postnatal BEP", "Azithromycin",
                 "Nicotinamide", "Nicotinamide +\nAzithromycin", "Postnatal BEP",
                 "Postnatal BEP +\nAzithromycin")))

# ELICIT and Mumta-LW randomized after birth, so their Birth rows are not estimable
plotdf_primary <- plotdf_primary %>%
  filter(!(study %in% c("ELICIT", "Mumta-LW") & age == "Birth"))

p <- ggplot(plotdf_primary,
            aes(x = contrast, y = est, group = Measure,
                color = Measure, shape = sig, alpha = sig)) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_linerange(aes(ymin = cil, ymax = ciu, linetype = Measure), position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  facet_grid(study ~ age, scale = "free_y") +
  scale_alpha_manual(values = c(0.7, 1)) + guides(alpha = "none") +
  # colourblind-safe: Okabe-Ito blue/orange plus a CI linetype per measure
  scale_color_manual(values = c("#0072B2", "#E69F00", "#009E73", "#CC79A7")) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash")) +
  scale_shape_manual(values = c(1, 19), name = "Significant") +
  # theme_imic() (Helvetica, font-size floors), with a panel border and plain
  # (unbolded) strip titles as in the main figures
  theme_imic(base_size = 8) +
  theme(panel.border    = element_rect(colour = "black", fill = NA, linewidth = 0.3),
        strip.text      = element_text(face = "plain"),
        legend.position = "bottom") +
  xlab("Intervention arm\n(compared to Control)") + ylab("Z-score difference")

save_figure_3way(p, "figureS1_growth_outcomes",
                 width  = science_dims$full_page$width,   # 7.25 in
                 height = 4.2,
                 dir    = paste0(here::here(), "/figures"))

message("Wrote figures/figureS1_growth_outcomes.{pdf,eps,png}")
