# fig6D-crosscompartment.R
# =============================================================================
# Manuscript Figure 6, Panel D: cross-compartment tracking of BEP-associated
# UPREGULATED metabolite features AND their KEGG pathways.
#
# Web-independent port of Trenton's "Compartment Tracking.Rmd" pathway-annotated
# figure ("Figure X. Cross-Compartment Tracking of BEP-Associated Metabolite
# Features and Pathways"). Each row is one upregulated metabolite that transfers
# across >= 2 distinct COMPARTMENTS. Features are grouped into one metabolite by
# Trenton's NAME-BASED rule (ct_trenton_track() in
# src/trenton-ports/compartment-tracking.R: tracking_name = the feature's own
# putative name, or "m/z <round4>"); dots are its average treatment effect in
# each compartment (coloured by compartment, shaped by BEP-supplement detection).
# The y-axis label carries the metabolite name plus the KEGG pathway(s) it
# drives, with each pathway's enrichment P-value.
#
# METHOD FIDELITY -- this panel previously used a union-find connected-component
# grouping (src/metaboanalyst/run-compartment-pathway.R), which is a
# GENERALISATION of Trenton's method, not his method -- it merges isobars he
# kept separate and put this panel on a DIFFERENT grouping rule than Table S11
# and Fig S10 (Manuscript/PROVENANCE_AUDIT.md §6.4). This panel now reads the
# same name-based rule, built by src/metaboanalyst/run-compartment-pathway-trenton.R.
#
# WEB INDEPENDENCE (the point of the original rebuild, unchanged here): Trenton's
# Rmd read a DOWNLOADED MetaboAnalyst "pathway_results.csv" and a HAND-CODED
# metabolite->pathway tribble. Both are replaced by a LOCAL MetaboAnalystR KEGG
# pathway analysis -- src/metaboanalyst/run-compartment-pathway-trenton.R --
# whose outputs we consume:
#   results/compartment_tracking/trenton_linked_crosscompartment.csv (linked plot data, >=2 compartments)
#   results/compartment_tracking/metabolite_pathways_trenton.csv     (ORA hit membership = the tribble)
# Pathway names/P-values match the web tool's KEGG library (metpa); the pathway SET
# differs slightly from his hand-curated tribble because the local KEGG name-matching,
# not manual selection, decides membership (documented in the figure comparison doc).
#
# The earlier BOTH-DIRECTION concordance dot-plot (95% CIs, ordered by #compartments)
# is preserved at fig6D-concordance-alt.R.
#
# Inputs (in-repo, produced by src/metaboanalyst/run-compartment-pathway-trenton.R):
#   results/compartment_tracking/trenton_linked_crosscompartment.csv
#   results/compartment_tracking/metabolite_pathways_trenton.csv
# Output: figures/figure6_panelD_crosscompartment.png
#
# Run from repo root: Rscript figure-scripts/manuscript_figures/fig6D-crosscompartment.R
# (build order: appendix-compartment-tracking.R -> run-compartment-pathway-trenton.R
#  -> appendix-compartment-tracking.R again -> this script; see the trenton script's header)
# =============================================================================
suppressMessages({ library(dplyr); library(stringr); library(ggplot2); library(forcats); library(readr) })
source("figure-scripts/0_figure-functions.R")   # theme_imic() = Science-submission theme (Helvetica, font floors)

LINKED <- "results/compartment_tracking/trenton_linked_crosscompartment.csv"
PATHS  <- "results/compartment_tracking/metabolite_pathways_trenton.csv"
OUT    <- "figures/figure6_panelD_crosscompartment.png"
# Milk untargeted ATE results -- the ONLY compartment whose raw_p is absent from the
# linked table; used to recover milk 95% CIs (carries est/cil/ciu per feature x visit).
MILK_ATE <- "results/adjusted_combined_arms_intervention_effects_untargeted_results_clean_ATE.RDS"
# linked `timepoint` codes -> milk-ATE `visit` labels (Misame milk only).
TP2VISIT <- c("1421d" = "14-21 days", "pn12" = "1-2 mo.", "pn34" = "3-4 mo.",
              "pn56" = "5-6 mo.", "1mo" = "1 mo.", "5mo" = "5 mo.")

# Compartment ordering from Trenton's Compartment Tracking.Rmd. Palette is
# colourblind-safe since 2026-09-23: his tableau map put Maternal VAMS (green) and
# Infant VAMS (red) at CIEDE2000 dE 0.7 under simulated deuteranopia. The new colours
# stay >= 27.9 apart under protan/deutan/tritan simulation, and each compartment also
# gets its own CI LINETYPE (plus a fixed dodge order within every row) as a
# colour-independent cue; shape stays reserved for supplement detection.
COMP_LEVELS <- c("Maternal plasma", "Maternal VAMS", "Milk", "Infant VAMS")
COMP_COLS   <- c("Maternal plasma" = "#004488", "Maternal VAMS" = "#56B4E9",
                 "Milk" = "#E69F00", "Infant VAMS" = "#000000")
COMP_LTY    <- c("Maternal plasma" = "solid", "Maternal VAMS" = "dashed",
                 "Milk" = "dotted", "Infant VAMS" = "twodash")
# Supplement-detection tiers (his `supplement_label`); fixed order -> fixed shapes.
SUPP_LEVELS <- c("Supplement-abundant", "Detected in supplement",
                 "Not reliably detected in supplement")
SUPP_SHAPES <- c("Supplement-abundant" = 17, "Detected in supplement" = 15,
                 "Not reliably detected in supplement" = 16)

build_panelD <- function(linked_path = LINKED, paths_path = PATHS, out_png = OUT) {
  linked <- read_csv(linked_path, show_col_types = FALSE) %>%
    filter(direction == "Upregulated", !is.na(effect_size))

  # --- pathway labels (his metaboanalyst_axis_labels), from the LOCAL ORA hits ---
  # one "Pathway (P = ...)" line per metabolite. A metabolite can map to several KEGG
  # pathways; we annotate each with only its MOST enriched (lowest-P) pathway (MAX_PATHS=1)
  # so the wrapped, multi-line y-axis labels do not overrun neighbouring rows in the tall
  # 29-row panel. (Citrate/(R)-Glycerate etc. map to 3+; the strongest is shown.)
  MAX_PATHS <- 1L
  paths <- read_csv(paths_path, show_col_types = FALSE) %>%
    group_by(tracking_name) %>% slice_min(raw_p, n = MAX_PATHS, with_ties = FALSE) %>%
    ungroup() %>%
    # wrap long pathway names onto multiple rows so the y-axis labels stay narrow
    mutate(pathway_label = str_wrap(paste0(pathway, " (P = ",
                                  format.pval(raw_p, digits = 2, eps = 0.001), ")"), width = 32))
  axis_labels <- paths %>% arrange(raw_p) %>% group_by(tracking_name) %>%
    summarise(pathway_text = paste(unique(pathway_label), collapse = "\n"), .groups = "drop") %>%
    mutate(feature_pathway_label = paste0(tracking_name, "\n", pathway_text))

  # attach labels; metabolites with no mapped pathway keep just their name.
  d <- linked %>%
    left_join(axis_labels, by = "tracking_name") %>%
    mutate(feature_pathway_label = ifelse(is.na(feature_pathway_label),
                                          as.character(tracking_name), feature_pathway_label),
           # order rows by each metabolite's strongest (max) upregulated effect
           feature_pathway_label = fct_reorder(feature_pathway_label, effect_size, .fun = max),
           compartment      = factor(compartment, levels = COMP_LEVELS),
           supplement_label = factor(supplement_label, levels = SUPP_LEVELS))

  # ONE dot per compartment: a compound significant in the same compartment at several
  # timepoints (e.g. N-Acetylserotonin in Infant VAMS at pn12/pn34/pn56) would otherwise
  # draw several same-colour dots. Keep its strongest (max effect) appearance per
  # compartment -- this panel is about cross-COMPARTMENT transfer, not the within-
  # compartment timepoint trajectory.
  d <- d %>% group_by(feature_pathway_label, tracking_name, compartment) %>%
    slice_max(effect_size, n = 1, with_ties = FALSE) %>% ungroup()

  # --- 95% CIs (restore Andrew's original design) -----------------------------
  # Non-milk compartments carry raw_p, so derive the SE from the z-statistic
  # (se = |effect|/z; CI = effect +/- 1.96 se). Milk has no raw_p in this table, so
  # recover its cil/ciu from the milk untargeted ATE results, matched by feature id
  # (case-insensitive) and visit.
  #
  # NUMERICAL NOTE -- do not "simplify" this back to qnorm(1 - raw_p/2).
  # This panel's strongest features have raw_p down to 3.6e-45. Once raw_p falls
  # below ~2.2e-16 (.Machine$double.eps), `1 - raw_p/2` rounds to exactly 1 and
  # qnorm(1) is Inf, so se became NA and geom_errorbar SILENTLY DROPPED the whole
  # interval -- octenoylcarnitine, citrate, 4-hydroxy-2-oxoglutarate and perillic
  # acid all lost visible bars 0.14-0.23 ATE units wide. Taking the upper tail
  # directly avoids that cancellation and is exact for arbitrarily small p.
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

  # Andrew's original geometry: points DODGED per compartment with horizontal 95% CIs
  # and NO connecting line (Trenton's slope line removed), keeping his pathway-annotated
  # y-axis + compartment colours + supplement shapes + the same axis/theme formatting.
  pd <- position_dodge(width = 0.6)
  p <- ggplot(d, aes(x = effect_size, y = feature_pathway_label,
                     color = compartment, group = compartment)) +
    geom_vline(xintercept = 0, color = "grey70", linewidth = 0.4) +
    geom_errorbarh(aes(xmin = cil, xmax = ciu, linetype = compartment), position = pd, height = 0,
                   linewidth = 0.4, alpha = 0.9) +
    geom_point(aes(shape = supplement_label), size = 2.4, alpha = 0.95, position = pd) +
    scale_color_manual(values = COMP_COLS, drop = FALSE) +
    scale_linetype_manual(values = COMP_LTY, drop = FALSE) +   # merges with the colour legend
    # shorten the long supplement labels so the legend fits on a single line
    # drop = TRUE: only show shape-legend keys actually present in the plotted
    # data (this panel's current data never uses "Supplement-abundant", so it
    # would otherwise show an empty/unused "Abundant" key)
    scale_shape_manual(values = SUPP_SHAPES, drop = TRUE,
                       labels = c("Supplement-abundant" = "Abundant",
                                  "Detected in supplement" = "Detected",
                                  "Not reliably detected in supplement" = "Not detected")) +
    labs(x = "Average Treatment Effect (95% CI)", y = NULL,
         color = "Compartment", linetype = "Compartment", shape = "BEP supplement") +
    # The 3-key shape legend fits on one row; the 4-key compartment legend does NOT
    # at this panel width -- on one row its last key ("Infant VAMS") is clipped by
    # the right edge. Wrap it to two rows (2x2) so every key stays inside the canvas.
    guides(shape = guide_legend(order = 1, nrow = 1, override.aes = list(size = 2)),
           color = guide_legend(order = 2, nrow = 2, override.aes = list(size = 2, linewidth = 0.6)),
           linetype = guide_legend(order = 2, nrow = 2)) +
    theme_imic(base_size = 9) +   # Science-submission theme (Helvetica, font floors)
    theme(# 29 wrapped, multi-line rows: 6 pt y label keeps the tall panel legible
          axis.text.y = element_text(size = 6, lineheight = 0.8, hjust = 0),
          panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
          # full box border (theme_imic() drops it by default), matching the other
          # Fig 6 panels (A/B/C) rather than the bottom-only axis.line.x this had before.
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3),
          legend.position = "bottom", legend.box = "vertical",
          legend.box.just = "left", legend.direction = "horizontal",
          # much smaller legend than the theme_imic floors
          legend.title = element_text(size = 6.5), legend.text = element_text(size = 6),
          legend.key.size = unit(3, "mm"), legend.spacing.y = unit(0, "mm"),
          legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(2, "mm"),
          plot.margin = margin(6, 8, 6, 6))

  # Right column of the full-page (7.25 in) Fig 6. Full page-height (9.5 in max per
  # 0_figure-functions.R) for the 29 multi-line pathway rows.
  ggsave(out_png, p, width = 5.0, height = 9.9, units = "in", dpi = 300,
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
