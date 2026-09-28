# =============================================================================
# fig6D-crosscompartment.R
#
# Builds Fig 6D: the BEP-associated up-regulated metabolite features in MISAME-III
# that are linked across at least two compartments (maternal plasma, maternal VAMS,
# milk, infant VAMS), with their KEGG pathways. Each row is one metabolite; features
# are grouped by their own putative name (unannotated features by "m/z <mass>"), the
# name-based rule of ct_name_track() in src/metaboanalyst/R/compartment-tracking.R,
# the same rule Table S11 uses. Points are the average treatment effect in each
# compartment with horizontal 95% CIs, coloured by compartment and shaped by
# BEP-supplement detection. The y-axis label carries the metabolite name and its most
# enriched KEGG pathway with that pathway's P-value, from a local MetaboAnalystR
# pathway run, so pathway membership is decided by KEGG name matching.
#
# Inputs:  results/compartment_tracking/linked_crosscompartment.csv  (plot rows)
#          results/compartment_tracking/metabolite_pathways.csv      (compound -> pathway hits)
#            (both written by src/metaboanalyst/run-compartment-pathway.R, which runs after
#            src/metaboanalyst/build-compartment-query-list.R)
#          results/adjusted_combined_arms_intervention_effects_untargeted_results_clean_ATE.RDS
#            (milk 95% CIs; src/2 analysis/clean_results.R)
# Outputs: figures/figure6_panelD_crosscompartment.png
# [needs an on-request file] the milk untargeted ATE RDS is a feature-level result
# file that is not shipped (size).
# Run from repo root: Rscript figure-scripts/manuscript_figures/fig6D-crosscompartment.R
# =============================================================================
suppressMessages({ library(dplyr); library(stringr); library(ggplot2); library(forcats); library(readr) })
source("figure-scripts/0_figure-functions.R")   # theme_imic() (Helvetica, font-size floors)
source("figure-scripts/manuscript_figures/fig6_layout.R")   # FIG6_* geometry + print-scale text sizes

# Y-axis labels: metabolite name in title case, black and bold; its pathway line in italics; both right-aligned against the axis. Drawn as
# ggtext markdown/HTML. Title case capitalises each all-lowercase word, including after a
# hyphen ("Sedoheptulose 1-Phosphate", "4-Hydroxy-2-Oxoglutarate"), and leaves mixed-case
# tokens alone ("dCMP", "(R)-", "m/z").
title_case_name <- function(x) {
  vapply(strsplit(x, "(?<=[ -])", perl = TRUE), function(parts) {
    paste0(vapply(parts, function(w)
      if (grepl("^[a-z]{2,}", w) && !grepl("[A-Z/]", w))
        paste0(toupper(substr(w, 1, 1)), substring(w, 2)) else w, ""), collapse = "")
  }, "")
}
md_name <- function(x) paste0("<b>", sub("^m/z ", "<i>m/z</i> ", title_case_name(x)), "</b>")

LINKED <- "results/compartment_tracking/linked_crosscompartment.csv"
PATHS  <- "results/compartment_tracking/metabolite_pathways.csv"
OUT    <- "figures/figure6_panelD_crosscompartment.png"
# Milk untargeted ATE results -- the only compartment whose raw_p is absent from the
# linked table; used to recover milk 95% CIs (carries est/cil/ciu per feature x visit).
MILK_ATE <- "results/adjusted_combined_arms_intervention_effects_untargeted_results_clean_ATE.RDS"
# linked `timepoint` codes -> milk-ATE `visit` labels (Misame milk only).
TP2VISIT <- c("1421d" = "14-21 days", "pn12" = "1-2 mo.", "pn34" = "3-4 mo.",
              "pn56" = "5-6 mo.", "1mo" = "1 mo.", "5mo" = "5 mo.")

# Compartment ordering as in the original R Markdown analysis. The colours stay
# >= 27.9 CIEDE2000 apart under simulated protan/deutan/tritan vision, and each
# compartment keeps a fixed dodge order within every row as a colour-independent cue;
# shape is reserved for supplement detection.
COMP_LEVELS <- c("Maternal plasma", "Maternal VAMS", "Milk", "Infant VAMS")
COMP_COLS   <- c("Maternal plasma" = "#004488", "Maternal VAMS" = "#56B4E9",
                 "Milk" = "#E69F00", "Infant VAMS" = "#000000")
# Supplement-detection tiers (`supplement_label`); fixed order -> fixed shapes.
SUPP_LEVELS <- c("Supplement-abundant", "Detected in supplement",
                 "Not reliably detected in supplement")
SUPP_SHAPES <- c("Supplement-abundant" = 17, "Detected in supplement" = 15,
                 "Not reliably detected in supplement" = 16)

build_panelD <- function(linked_path = LINKED, paths_path = PATHS, out_png = OUT) {
  linked <- read_csv(linked_path, show_col_types = FALSE) %>%
    filter(direction == "Upregulated", !is.na(effect_size))

  # --- pathway labels, from the local pathway-run hits ---
  # one "Pathway (P = ...)" line per metabolite. A metabolite can map to several KEGG
  # pathways; each is annotated with only its most enriched (lowest-P) pathway
  # (MAX_PATHS = 1) so the wrapped, multi-line y-axis labels do not overrun neighbouring
  # rows. (Citrate, (R)-Glycerate etc. map to 3+; the strongest is shown.)
  MAX_PATHS <- 1L
  paths <- read_csv(paths_path, show_col_types = FALSE) %>%
    group_by(tracking_name) %>% slice_min(raw_p, n = MAX_PATHS, with_ties = FALSE) %>%
    ungroup() %>%
    # wrap long pathway names onto multiple rows so the y-axis labels stay narrow
    # non-breaking spaces keep "P = 0.0631" together when the line wraps
    mutate(pathway_label = str_wrap(paste0(pathway, " (P\u00a0=\u00a0",
                                  format.pval(raw_p, digits = 2, eps = 0.001), ")"), width = 32))
  axis_labels <- paths %>% arrange(raw_p) %>% group_by(tracking_name) %>%
    summarise(pathway_text = paste(unique(pathway_label), collapse = "\n"), .groups = "drop") %>%
    mutate(feature_pathway_label = paste0(md_name(tracking_name), "<br><i>",
                                          gsub("\n", "<br>", pathway_text), "</i>"))

  # attach labels; metabolites with no mapped pathway keep just their name.
  d <- linked %>%
    left_join(axis_labels, by = "tracking_name") %>%
    mutate(feature_pathway_label = ifelse(is.na(feature_pathway_label),
                                          md_name(as.character(tracking_name)), feature_pathway_label),
           # order rows by each metabolite's strongest (max) upregulated effect
           feature_pathway_label = fct_reorder(feature_pathway_label, effect_size, .fun = max),
           compartment      = factor(compartment, levels = COMP_LEVELS),
           supplement_label = factor(supplement_label, levels = SUPP_LEVELS))

  # One dot per compartment: a compound significant in the same compartment at several
  # timepoints (e.g. N-Acetylserotonin in Infant VAMS at pn12/pn34/pn56) would otherwise
  # draw several same-colour dots. Keep its strongest (max effect) appearance per
  # compartment -- this panel is about cross-compartment transfer, not the within-
  # compartment timepoint trajectory.
  d <- d %>% group_by(feature_pathway_label, tracking_name, compartment) %>%
    slice_max(effect_size, n = 1, with_ties = FALSE) %>% ungroup()

  # --- 95% CIs -------------------------------------------------------------------
  # Non-milk compartments carry raw_p, so derive the SE from the z-statistic
  # (se = |effect|/z; CI = effect +/- 1.96 se). Milk has no raw_p in this table, so
  # recover its cil/ciu from the milk untargeted ATE results, matched by feature id
  # (case-insensitive) and visit.
  #
  # Numerical note: do not write this as qnorm(1 - raw_p/2). This panel's strongest
  # features have raw_p down to 3.6e-45; below ~2.2e-16 (.Machine$double.eps),
  # `1 - raw_p/2` rounds to exactly 1, qnorm(1) is Inf, se becomes NA and
  # geom_errorbar silently drops the interval. Taking the upper tail directly avoids
  # that cancellation and is exact for arbitrarily small p.
  # raw_p is clamped away from 0 so an exactly-zero p yields a finite (very large) z
  # rather than a zero-width interval.
  z  <- stats::qnorm(pmax(d$raw_p, .Machine$double.xmin) / 2, lower.tail = FALSE)
  se <- ifelse(is.finite(z) & z > 0, abs(d$effect_size) / z, NA_real_)
  d$cil <- d$effect_size - 1.96 * se
  d$ciu <- d$effect_size + 1.96 * se
  milk_ci <- readRDS(MILK_ATE)
  milk_ci <- milk_ci[milk_ci$contrast == "BEP", c("biomarker", "visit", "cil", "ciu")]
  milk_ci <- milk_ci %>% transmute(mkey = tolower(biomarker), visit = visit,
                                   m_cil = cil, m_ciu = ciu) %>%
    distinct(mkey, visit, .keep_all = TRUE)
  d <- d %>%
    mutate(mkey = tolower(feature), visit = unname(TP2VISIT[timepoint])) %>%
    left_join(milk_ci, by = c("mkey", "visit")) %>%
    mutate(cil = ifelse(compartment == "Milk" & !is.na(m_cil), m_cil, cil),
           ciu = ifelse(compartment == "Milk" & !is.na(m_ciu), m_ciu, ciu)) %>%
    select(-mkey, -visit, -m_cil, -m_ciu)

  # Points dodged per compartment with horizontal 95% CIs and no connecting line;
  # pathway-annotated y-axis, compartment colours, supplement shapes.
  pd <- position_dodge(width = 0.6)
  p <- ggplot(d, aes(x = effect_size, y = feature_pathway_label,
                     color = compartment, group = compartment)) +
    geom_vline(xintercept = 0, color = "grey70", linewidth = 0.4) +
    geom_errorbarh(aes(xmin = cil, xmax = ciu), position = pd, height = 0,
                   linewidth = 0.4, alpha = 0.9) +
    geom_point(aes(shape = supplement_label), size = 2.4, alpha = 0.95, position = pd) +
    scale_color_manual(values = COMP_COLS, drop = FALSE) +
    # shorten the long supplement labels so the legend fits on a single line
    # drop = TRUE: only show shape-legend keys actually present in the plotted
    # data (this panel's current data never uses "Supplement-abundant", so it
    # would otherwise show an empty/unused "Abundant" key)
    scale_shape_manual(values = SUPP_SHAPES, drop = TRUE,
                       labels = c("Supplement-abundant" = "Abundant",
                                  "Detected in supplement" = "Detected",
                                  "Not reliably detected in supplement" = "Not detected")) +
    labs(x = "Average Treatment Effect (95% CI)", y = NULL,
         color = "Compartment", shape = "BEP supplement") +
    # The 3-key shape legend fits on one row; the 4-key compartment legend does not
    # at this panel width (its last key is clipped by the right edge), so it wraps to
    # two rows (2x2).
    guides(shape = guide_legend(order = 1, nrow = 1, override.aes = list(size = 2)),
           color = guide_legend(order = 2, nrow = 2, override.aes = list(size = 2, linewidth = 0.6))) +
    theme_imic(base_size = 9) +   # shared theme (Helvetica, font-size floors)
    theme(# name (bold) + pathway (italic) labels, right-aligned against the axis; sizes
          # are print-scale corrected (fig6_layout.R) so D prints at the A-C/Figs 3, 5 sizes
          axis.text.y = ggtext::element_markdown(size = fig6_pt(6.5, FIG6_D_SCALE), colour = "black",
                                                 lineheight = 1.0, hjust = 1, halign = 1,
                                                 margin = margin(r = 3)),   # gap from the border
          axis.text.x = element_text(size = fig6_pt(FIG6_TICK_PT, FIG6_D_SCALE)),
          axis.title  = element_text(size = fig6_pt(FIG6_TITLE_PT, FIG6_D_SCALE)),
          # no horizontal row lines (theme_imic()); full box border (theme_imic() drops
          # it by default), matching the other Fig 6 panels
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3),
          legend.position = "bottom", legend.box = "vertical",
          legend.box.just = "left", legend.direction = "horizontal",
          # much smaller legend than the theme_imic floors
          legend.title = element_text(size = fig6_pt(FIG6_LEGEND_PT, FIG6_D_SCALE)),
          legend.text  = element_text(size = fig6_pt(FIG6_LEGEND_PT, FIG6_D_SCALE)),
          legend.key.size = unit(3, "mm"), legend.spacing.y = unit(0, "mm"),
          legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(2, "mm"),
          plot.margin = margin(6, 8, 6, 6))

  # Right column of the full-page (7.25 in) Fig 6, tall enough for the 20 multi-line
  # pathway-annotated rows (fig6-composite.R scales it to the column height).
  ggsave(out_png, p, width = FIG6_D_W_IN, height = FIG6_D_H_IN, units = "in", dpi = 300,
         bg = "white", device = ragg::agg_png)
  # Guard against the silent-drop failure above ever returning: a missing interval
  # is invisible in the rendered panel, so it has to be reported at build time.
  n_missing_ci <- sum(is.na(d$cil) | is.na(d$ciu))
  if (n_missing_ci > 0) {
    warning(n_missing_ci, " of ", nrow(d), " points have no 95% CI and will draw ",
            "as bare points:\n",
            paste0("  - ", d$tracking_name[is.na(d$cil) | is.na(d$ciu)], " (",
                   d$compartment[is.na(d$cil) | is.na(d$ciu)], ")", collapse = "\n"),
            call. = FALSE)
  }
  cat("  points with a 95% CI:", nrow(d) - n_missing_ci, "/", nrow(d), "\n")
  cat("  x range plotted (incl. CIs): ",
      sprintf("%.3f .. %.3f\n", min(d$cil, d$effect_size, na.rm = TRUE),
              max(d$ciu, d$effect_size, na.rm = TRUE)))
  cat("wrote", out_png, "| compounds:", nlevels(d$feature_pathway_label),
      "| pathway-annotated:", sum(!is.na(match(levels(d$feature_pathway_label),
        axis_labels$feature_pathway_label))), "\n")
  invisible(p)
}

if (sys.nframe() == 0) invisible(build_panelD())
