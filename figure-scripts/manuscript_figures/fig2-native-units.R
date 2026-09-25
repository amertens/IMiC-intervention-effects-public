# fig2-native-units.R
# =============================================================================
# Native-unit companion to Fig. 2 (primary-outcome forest). Fig. 2 plots Z-scored
# average treatment effects so every nutrient shares one axis; here each nutrient gets
# its own small panel on its own native-unit axis (units differ, so they cannot share
# one), with the seven trial-visits as rows. Requested 2026-09-25 (Reviewer 1, comment 4e
# asked where native-unit estimates are shown).
#
# Estimates = the native-unit ATEs in Table S1's machine-readable file
# (results/tables/table_s1_primary_secondary_native_units_wide.csv: ATE, 95% CI, P, q),
# i.e. the same numbers the Results text quotes. Significance tiers as in Fig. 2:
# grey = not significant, blue open = P < 0.05 only, orange filled = q < 0.05.
#
# Units come from metadata/Milk_Component_Spec_IMiC_*.csv (UNIT column) and the unit
# conversions in src/1 data prep/6-milq-adequecy.R; Cr, Mo, Mn and As have no recorded
# unit and are labelled ug/L from their magnitudes (flagged in UNIT_INFERRED).
#
# Output: figures/figure2_native_units.{pdf,eps,png}
# Run from repo root: Rscript figure-scripts/manuscript_figures/fig2-native-units.R
# =============================================================================
suppressMessages({library(dplyr); library(ggplot2); library(cowplot)})
root <- here::here()
source(file.path(root, "figure-scripts/0_figure-functions.R"))

UNIT <- c(
  kcal.l = "kcal/L", fat = "g/100 mL", protein = "g/100 mL", cho = "g/100 mL",
  vitamin.a = "mg/L", a.tocopherol = "mg/L", g.tocopherol = "mg/L",
  ca = "mg/L", mg = "mg/L", k = "mg/L", na = "mg/L", p = "mg/L",
  cu = "µg/L", fe = "µg/L", zn = "µg/L", se = "µg/L",
  mn = "µg/L", mo = "µg/L", cr = "µg/L", as = "µg/L",
  b12 = "pmol/L")
UNIT_INFERRED <- c("mn", "mo", "cr", "as")
SHORT <- c(kcal.l = "Energy", fat = "Fat", protein = "Protein", cho = "Carbohydrate")          # no unit in the spec files
unit_of <- function(code) {
  u <- unname(UNIT[tolower(code)])
  ifelse(is.na(u), "µg/L", u)                  # every remaining primary outcome is a B-vitamin (ug/L)
}

# MISAME-III codes carbohydrate as "total carbohydrate", the other trials as "cho"; Fig. 2
# shows them as one row, and both are g/100 mL (control means 7.1 and 6.8).
norm_code <- function(x) { x <- tolower(x); ifelse(x %in% c("total carbohydrate", "carbohydrate"), "cho", x) }

# Panel membership and order: exactly Fig. 2's three panels, top-to-bottom order of
# Fig. 2 (largest mean Z-scored effect first).
panel_of <- function(f, panel) {
  x <- read.csv(file.path(root, "results/subsetted results", f))
  x %>% mutate(biomarker = norm_code(biomarker)) %>%
    group_by(biomarker) %>% summarise(z = mean(est), .groups = "drop") %>%
    arrange(desc(z)) %>% mutate(panel = panel, ord = row_number())
}
panels <- bind_rows(panel_of("primary_macro.csv", "A"),
                    panel_of("primary_micro.csv", "B"),
                    panel_of("primary_bvit.csv",  "C"))

nat <- read.csv(file.path(root, "results/tables/table_s1_primary_secondary_native_units_wide.csv")) %>%
  filter(outcome_class == "primary") %>%
  mutate(code = norm_code(biomarker))

# consistency check: these ATEs are the est_unscaled values behind Fig. 2
chk <- bind_rows(lapply(c("primary_macro.csv", "primary_micro.csv", "primary_bvit.csv"), function(f)
  read.csv(file.path(root, "results/subsetted results", f)))) %>%
  transmute(code = norm_code(biomarker), study, visit, est_unscaled) %>%
  inner_join(nat %>% select(code, study, visit, ATE), by = c("code", "study", "visit"))
stopifnot(nrow(chk) > 0, max(abs(chk$ATE - chk$est_unscaled) / pmax(1e-9, abs(chk$ATE))) < 1e-6)

STUDY <- c(Vital = "Mumta-LW", Elicit = "ELICIT", Misame = "MISAME-III")
ROWS  <- c("Mumta-LW 1.5 mo.", "Mumta-LW 2 mo.", "ELICIT 1 mo.", "ELICIT 5 mo.",
           "MISAME-III 14-21 days", "MISAME-III 1-2 mo.", "MISAME-III 3-4 mo.")

d <- nat %>%
  inner_join(panels %>% mutate(code = norm_code(biomarker)) %>% select(code, panel, ord), by = "code") %>%
  mutate(row  = factor(paste(STUDY[study], visit), levels = rev(ROWS)),
         tier = factor(case_when(pval_adj < 0.05 ~ "Sig",
                                 pval < 0.05     ~ "Sig before FDR",
                                 TRUE            ~ "Not Significant"),
                       levels = c("Not Significant", "Sig before FDR", "Sig")),
         name = paste0(ifelse(code %in% names(SHORT), SHORT[code], abbr_label(code)),
                       " (", unit_of(code), ")",
                       ifelse(code %in% UNIT_INFERRED, "*", "")))

forest <- function(dd) {
  lv <- dd %>% distinct(ord, name) %>% arrange(ord) %>% pull(name)
  bold <- function(x) paste0("bold(", plotmath_label(x), ")")   # plotmath strips ignore strip face
  dd$facet <- factor(bold(dd$name), levels = bold(lv))
  ggplot(dd, aes(x = ATE, y = row, colour = tier, shape = tier)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
    geom_errorbarh(aes(xmin = ATE_cil, xmax = ATE_ciu), height = 0, linewidth = 0.4) +
    geom_point(size = 1.4) +
    facet_wrap(~ facet, ncol = 5, scales = "free_x", labeller = label_parsed) +
    scale_colour_manual(values = c("Not Significant" = "grey60", "Sig before FDR" = tableau10[1],
                                   "Sig" = tableau10[2]), drop = FALSE, name = NULL) +
    scale_shape_manual(values = c("Not Significant" = 1, "Sig before FDR" = 1, "Sig" = 19),
                       drop = FALSE, name = NULL) +
    scale_x_continuous(n.breaks = 3, labels = scales::label_number(big.mark = ",")) +
    expand_limits(x = 0) +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 7) +
    theme(strip.background = element_blank(), strip.text = element_text(size = 6.5, face = "bold"),
          panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
          axis.text = element_text(size = 6), legend.position = "none",
          panel.spacing.x = unit(3, "mm"), plot.margin = margin(2, 6, 2, 2))
}

pA <- forest(filter(d, panel == "A")); pB <- forest(filter(d, panel == "B")); pC <- forest(filter(d, panel == "C"))
nrow_of <- function(p) ceiling(length(unique(p$data$facet)) / 5)
g   <- ggplotGrob(pA + theme(legend.position = "bottom", legend.text = element_text(size = 7)) +
                    guides(colour = guide_legend(override.aes = list(size = 2))))
leg <- g$grobs[[which(g$layout$name %in% c("guide-box-bottom", "guide-box"))[1]]]
note <- ggdraw() + draw_label(paste0("Average treatment effect in native units (95% CI). ",
  "* unit inferred from measured values (no unit recorded for Cr, Mo, Mn, As)."),
  size = 6.5, x = 0.01, hjust = 0)
fig <- plot_grid(pA, pB, pC, leg, note, ncol = 1, labels = c("A", "B", "C", "", ""),
                 rel_heights = c(nrow_of(pA) + 0.25, nrow_of(pB), nrow_of(pC), 0.35, 0.25) )

save_figure_3way(fig, "figure2_native_units", width = 7.25,
                 height = 0.95 * (nrow_of(pA) + nrow_of(pB) + nrow_of(pC)) + 0.9,
                 dir = file.path(root, "figures"))
cat("wrote figures/figure2_native_units.{pdf,eps,png} |", n_distinct(d$code), "outcomes,",
    nrow(d), "estimates; ATE check vs Fig 2 est_unscaled: max rel diff",
    signif(max(abs(chk$ATE - chk$est_unscaled) / pmax(1e-9, abs(chk$ATE))), 2), "\n")
