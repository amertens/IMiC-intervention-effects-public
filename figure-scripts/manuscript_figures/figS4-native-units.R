# figS4-native-units.R
# =============================================================================
# Native-unit companion to Fig. S4 (MISAME-III milk triglyceride means by arm, without
# and with adjustment for total fat). Fig. S4 plots Z-scored arm means; this plots them
# in concentration units, on a log scale because species span ~0.02 to ~2,000 uM.
# Requested 2026-09-25 (Reviewer 1, comment 4e).
#
#  * Unadjusted (left block): per-arm model-based means in uM straight from the unscaled
#    per-arm TMLE run (results/adjusted_intervention_effects_unscaled_results_clean.RDS).
#  * Per unit fat (right block): the "fat-adjusted" analysis
#    (src/2 analysis/10_adjusted_analysis_metabolomics_fat_scaled.R) divides every
#    component by total fat and was only run Z-scored. run_bioTMLE(scale = TRUE) Z-scores
#    each outcome within study x visit with base scale(), so the arm means are put back
#    in native units as  mean_native = Z * SD + mean, using the SD and mean of TG / fat
#    in the same study x visit data. That back-transform is VALIDATED below on the
#    unadjusted analysis, where both the Z-scored and the native-unit runs exist.
#    Units: uM per g/100 mL of milk fat.
# Significance colours (vs control) as in Fig. S4.
#
# Output: figures/figureS4_native_units.{pdf,eps,png}
# Run from repo root: Rscript figure-scripts/manuscript_figures/figS4-native-units.R
# =============================================================================
suppressMessages({library(dplyr); library(ggplot2); library(cowplot)})
root <- here::here()
source(file.path(root, "figure-scripts/0_figure-functions.R"))
load(file.path(root, "metadata/milk_component.Rdata"))

dfull <- readRDS(file.path(root, "data/merged_analysis_datasets.RDS"))
tg_cols <- grep("^tg\\.", names(dfull), value = TRUE)
# as in Fig. S4: drop triglycerides recoded to binary under sparse detection
ctl1 <- dfull %>% filter(arm == "Control", study == "Misame", visit == "1")
binary_tg <- tg_cols[sapply(ctl1[tg_cols], function(x) length(unique(x)) %in% c(1, 2))]
keep_tg <- setdiff(tg_cols, binary_tg)

ARM_LEVELS <- c("Control", "BEP/IFA", "IFA/BEP", "BEP/BEP")
ARM_LABELS <- c("Control", "Prenatal BEP", "Postnatal BEP", "Pre+Postnatal BEP")
VISITS <- c("1" = "14-21 days", "2" = "1-2 mo.", "3" = "3-4 mo.")

# per study x visit mean and SD of each outcome, exactly as scale() computed them
scale_params <- function(df) {
  df %>% filter(study == "Misame") %>% group_by(visit) %>%
    summarise(across(all_of(keep_tg), list(m = ~ mean(.x, na.rm = TRUE), s = ~ sd(.x, na.rm = TRUE))),
              .groups = "drop") %>%
    tidyr::pivot_longer(-visit, names_to = c("biomarker", ".value"), names_pattern = "(.*)_(m|s)$") %>%
    mutate(visit = unname(VISITS[as.character(visit)]))   # data codes visits 1/2/3; results use labels
}
arm_means <- function(res) {
  res <- res %>% filter(study == "Misame", tolower(biomarker) %in% keep_tg) %>%
    mutate(biomarker = tolower(biomarker))
  mn  <- res %>% filter(measure == "MN") %>% select(visit, contrast, biomarker, est, cil, ciu)
  ate <- res %>% filter(measure == "ATE") %>% select(visit, contrast, biomarker, pval, pval_adj)
  left_join(mn, ate, by = c("visit", "contrast", "biomarker"))
}
tier <- function(pval, pval_adj) factor(case_when(
  !is.na(pval_adj) & pval_adj < 0.05 ~ "Significant after FDR",
  !is.na(pval) & pval < 0.05          ~ "Significant before FDR",
  TRUE                                ~ "Not Significant"),
  levels = c("Not Significant", "Significant before FDR", "Significant after FDR"))

# ---- unadjusted: native-unit run ----------------------------------------------
unadj <- arm_means(readRDS(file.path(root, "results/adjusted_intervention_effects_unscaled_results_clean.RDS")))

# ---- validate the back-transform on the unadjusted analysis -------------------
z_un <- arm_means(readRDS(file.path(root, "results/adjusted_intervention_effects_results_clean.RDS")))
p_un <- scale_params(dfull)
val <- z_un %>% inner_join(p_un, by = c("visit", "biomarker")) %>%
  mutate(back = est * s + m) %>%
  inner_join(unadj %>% select(visit, contrast, biomarker, native = est),
             by = c("visit", "contrast", "biomarker")) %>%
  mutate(rel = abs(back - native) / abs(native))
cat(sprintf("back-transform check (unadjusted, %d arm means): median rel. diff %.2g, 95th pct %.2g, max %.2g\n",
            nrow(val), median(val$rel), quantile(val$rel, 0.95), max(val$rel)))
stopifnot(median(val$rel) < 0.01)

# ---- per unit fat: back-transformed Z-scored run --------------------------------
dfat <- dfull
dfat[keep_tg] <- dfat[keep_tg] / dfat$fat                   # same division as script 10
p_fat <- scale_params(dfat)
fat <- arm_means(readRDS(file.path(root, "results/fat_adjusted_metabolomics_intervention_effects_results.RDS"))) %>%
  inner_join(p_fat, by = c("visit", "biomarker")) %>%
  mutate(est = est * s + m, cil = cil * s + m, ciu = ciu * s + m) %>%
  select(-m, -s)

# ---- plot ------------------------------------------------------------------------
# one TG order for every panel: ascending control-arm concentration at 14-21 days
ord <- unadj %>% filter(contrast == "Control", visit == "14-21 days") %>% arrange(est) %>% pull(biomarker)
prep <- function(x) x %>% mutate(
  tier = tier(pval, pval_adj),
  arm  = factor(contrast, levels = ARM_LEVELS, labels = ARM_LABELS),
  vis  = factor(visit, levels = VISITS),
  tg   = factor(biomarker, levels = ord))
FLOOR <- 0.01
panel <- function(x, ylab) {
  ggplot(prep(x) %>% filter(!is.na(arm), !is.na(tg)),
         aes(x = tg, y = est, colour = tier, alpha = tier)) +
    geom_linerange(aes(ymin = pmax(cil, FLOOR), ymax = ciu), linewidth = 0.25) +
    geom_point(size = 0.35) +
    facet_grid(arm ~ vis) +
    scale_y_log10(labels = scales::label_number(drop0trailing = TRUE, big.mark = ","),
                  oob = scales::squish) +
    scale_colour_manual(values = c("Not Significant" = "grey50", "Significant before FDR" = tableau10[1],
                                   "Significant after FDR" = tableau10[2]), drop = FALSE, name = NULL) +
    scale_alpha_manual(values = c(0.3, 0.5, 0.85), drop = FALSE, name = NULL) +
    labs(x = "Triglyceride species (ordered by control-arm concentration at 14-21 days)", y = ylab) +
    theme_bw(base_size = 8) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          strip.background = element_blank(), strip.text = element_text(face = "bold", size = 7.5),
          legend.position = "none",
          plot.margin = margin(t = 16, r = 4, b = 2, l = 2))   # room for the A/B block labels
}
pa <- panel(unadj, "Mean concentration, µM (log scale)")
pb <- panel(fat,   "Mean per unit fat, µM per g/100 mL fat (log scale)")
g   <- ggplotGrob(pa + theme(legend.position = "bottom") +
                    guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))))
leg <- g$grobs[[which(g$layout$name %in% c("guide-box-bottom", "guide-box"))[1]]]
fig <- plot_grid(plot_grid(pa, pb, nrow = 1, labels = c("A  Unadjusted", "B  Per unit milk fat"),
                           label_size = 9, hjust = 0, label_x = 0.01),
                 leg, ncol = 1, rel_heights = c(1, 0.06))
save_figure_3way(fig, "figureS4_native_units", width = 10, height = 7,
                 dir = file.path(root, "figures"))
cat("wrote figures/figureS4_native_units.{pdf,eps,png} |", length(ord), "triglycerides\n")
