# =============================================================================
# figS11-crosscompartment-volcanoes.R
#
# Supplementary Figure S11: volcano plots of BEP intervention effects on the
# metabolome across compartments and visits (MISAME-III).
#
# WHY THIS SCRIPT EXISTS
#   The shipped Fig. S11 (figures/cross_compartment/trenton_xcompartment_volcanoes.png)
#   is, like Fig. S10, an ORPHAN BINARY from Trenton's analysis with no generator in
#   this repo -- nothing could re-derive it after the 2026-07-13 blood TMLE re-run.
#   This script rebuilds it from the in-repo bioTMLE results.
#
# FIDELITY TO THE SHIPPED FIGURE
#   Reproduced: the 12-panel compartment x visit grid and its layout (maternal
#   plasma 3, maternal VAMS 2, human milk 3, infant VAMS 4), free per-panel scales,
#   -log10(raw P) vs ATE, the three-tier colour encoding (not significant / nominal
#   P < 0.05 / BH FDR < 0.05), the zero reference line, and repelled labels on
#   putatively annotated UP-regulated features.
#
#   Fidelity limit, stated rather than hidden: the shipped figure labels some
#   features this repo cannot name. Its label source was Trenton's own mummichog
#   annotation pass; the in-repo equivalent is
#   results/supplement_status_fdr_features.csv (script 54), which carries putative
#   names only for FDR-significant features. Panels whose FDR-significant features
#   are unnamed there will show fewer labels than the shipped PNG. No label is
#   invented to close that gap.
#
# Inputs:
#   results/blood_compartment_adjusted_combined_arms_intervention_effects_results_clean.RDS
#   results/adjusted_combined_arms_intervention_effects_untargeted_results_clean_ATE.RDS
#   results/supplement_status_fdr_features.csv        (putative names)
# Output:
#   figures/cross_compartment/figureS11_crosscompartment_volcanoes.png
#
# Run from repo root:
#   Rscript figure-scripts/manuscript_figures/figS11-crosscompartment-volcanoes.R
# =============================================================================
suppressMessages({ library(data.table); library(ggplot2); library(ggrepel)
                   library(patchwork) })
source("figure-scripts/0_figure-functions.R")   # theme_imic()

BLOOD <- "results/blood_compartment_adjusted_combined_arms_intervention_effects_results_clean.RDS"
MILK  <- "results/adjusted_combined_arms_intervention_effects_untargeted_results_clean_ATE.RDS"
NAMES <- "results/supplement_status_fdr_features.csv"
OUT   <- "figures/cross_compartment/figureS11_crosscompartment_volcanoes.png"

# Labels are allocated PER DIRECTION, not per panel. Labelling only up-regulated
# features (as the original figure did) misrepresents several panels: the
# FDR-significant split is 334 down vs 87 up at human milk 1-2 mo, 327 down vs 266 up
# at maternal plasma 1-2 mo, and 220 down vs 151 up at human milk 3-4 mo. A reader
# scanning only up-regulated labels would infer a one-sided intervention effect.
# Both directions get the same allowance and identical label styling, so neither is
# visually privileged.
MAX_LABELS_PER_DIR <- 4L
LABEL_CHARS <- 30L    # putative names run to 40+ chars in milk; truncate so labels fit
TIER_COLS  <- c("Not significant" = "grey25",
                "P < 0.05"        = "#4C9BD4",
                "FDR < 0.05"      = "#F5A623")

# panel key -> (display title, source compartment, source visit code)
# NB: the column is `panel_key`, not `key` -- data.table() treats a `key` argument
# as the sort key and would error on these values.
# Order is the shipped figure's reading order; the layout below reproduces its
# ragged grid (3 / 2 / 3 / 4) rather than packing the panels into full rows.
PANELS <- data.table(
  panel_key = c("MP_incl","MP_tri3","MP_pn12","MV_tri3","MV_pn56",
            "MK_1421","MK_pn12","MK_pn34",
            "IV_acco","IV_pn12","IV_pn34","IV_pn56"),
  # Two-line panel titles: compartment on the first row, visit on the second. On one
  # row the wider titles ("Maternal Plasma · Third Trimester" next to "Maternal
  # Plasma · 1–2 Months") run into each other at this panel width.
  title = c("Maternal Plasma\nInclusion","Maternal Plasma\nThird Trimester",
            "Maternal Plasma\n1–2 Months","Maternal VAMS\nThird Trimester",
            "Maternal VAMS\n5–6 Months","Human Milk\n14–21 Days",
            "Human Milk\n1–2 Months","Human Milk\n3–4 Months",
            "Infant VAMS\nBirth","Infant VAMS\n1–2 Months",
            "Infant VAMS\n3–4 Months","Infant VAMS\n5–6 Months"),
  src   = c(rep("blood", 5), rep("milk", 3), rep("blood", 4)),
  dataset = c(rep("MaternalPlasma", 3), rep("VamsPostnatalMaternal", 2),
              rep("Misame", 3), rep("VamsPostnatalInfant", 4)),
  visit = c("incl","tri3","pn12","tri3","pn56",
            "14-21 days","1-2 mo.","3-4 mo.",
            "acco","pn12","pn34","pn56"),
  # the compartment / timepoint spelling used by supplement_status_fdr_features.csv
  nm_compartment = c(rep("Maternal plasma", 3), rep("Maternal VAMS", 2),
                     rep("Milk", 3), rep("Infant VAMS", 4)),
  nm_timepoint   = c("incl","tri3","pn12","tri3","pn56",
                     "1421d","pn12","pn34","acco","pn12","pn34","pn56"))

load_inputs <- function() {
  b <- as.data.table(readRDS(BLOOD))[measure == "ATE" & contrast == "BEP",
        .(src = "blood", dataset, visit, feature = biomarker, est, pval, pval_adj)]
  m <- as.data.table(readRDS(MILK))[contrast == "BEP" & study == "Misame" &
        grepl("^rlc_", tolower(biomarker)),
        .(src = "milk", dataset = "Misame", visit, feature = biomarker, est, pval, pval_adj)]
  d <- rbind(b, m)
  d <- d[is.finite(est) & is.finite(pval)]
  d[, tier := fifelse(!is.na(pval_adj) & pval_adj < 0.05, "FDR < 0.05",
              fifelse(pval < 0.05, "P < 0.05", "Not significant"))]
  d[, tier := factor(tier, levels = names(TIER_COLS))]
  d[, fkey := tolower(feature)]

  nm <- fread(NAMES)[!is.na(putative_name) & trimws(putative_name) != "" &
                       putative_name != "Acid",
                     .(nm_compartment = compartment, nm_timepoint = timepoint,
                       fkey = tolower(feature), putative_name)]
  nm <- unique(nm, by = c("nm_compartment", "nm_timepoint", "fkey"))
  list(d = d, nm = nm)
}

one_panel <- function(i, d, nm) {
  p <- PANELS[i]
  dd <- d[src == p$src & dataset == p$dataset & visit == p$visit]
  if (!nrow(dd)) return(NULL)
  dd[, y := -log10(pmax(pval, .Machine$double.xmin))]

  lab <- merge(dd[tier == "FDR < 0.05"],
               nm[nm_compartment == p$nm_compartment & nm_timepoint == p$nm_timepoint,
                  .(fkey, putative_name)], by = "fkey")
  if (nrow(lab)) {
    lab[, lab_dir := ifelse(est > 0, "up", "down")]
    lab <- lab[order(pval)][, head(.SD, MAX_LABELS_PER_DIR), by = lab_dir]
    lab[, putative_name := ifelse(nchar(putative_name) > LABEL_CHARS,
                                  paste0(substr(putative_name, 1, LABEL_CHARS - 1), "…"),
                                  putative_name)]
  }

  g <- ggplot(dd, aes(est, y)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, colour = "grey20") +
    geom_point(aes(colour = tier), size = 0.55, alpha = 0.75) +
    scale_colour_manual(values = TIER_COLS, drop = FALSE, name = NULL) +
    labs(title = p$title, x = NULL, y = NULL) +
    theme_imic(base_size = 8) +
    theme(# bold, left-justified title, matching the other volcano-grid figures (3A/5A)
          plot.title = element_text(size = 8, face = "bold", hjust = 0,
                                    lineheight = 1.05, margin = margin(b = 2)),
          panel.grid.major.y = element_line(colour = "grey93", linewidth = 0.25),
          axis.text = element_text(size = 6.5), legend.position = "none",
          # black panel border, matching the other volcano-grid figures (3A/5A/6A)
          # and the submitted S9 -- theme_imic() drops the border by default.
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3),
          plot.margin = margin(3, 4, 3, 3))
  if (nrow(lab))
    g <- g + geom_text_repel(data = lab, aes(label = putative_name), size = 1.9,
                             colour = "grey15", segment.colour = "grey55",
                             segment.size = 0.25, min.segment.length = 0,
                             max.overlaps = Inf, box.padding = 0.35,
                             point.padding = 0.1, seed = 1L)
  # carried so the build can report the up/down label split and the asymmetry
  # this fixed cannot silently return
  attr(g, "imic_label_split") <- c(up = sum(lab$est > 0), down = sum(lab$est < 0))
  g
}

build_figS11 <- function(out_png = OUT) {
  inp <- load_inputs()
  panels <- lapply(seq_len(nrow(PANELS)), one_panel, d = inp$d, nm = inp$nm)
  missing <- PANELS$panel_key[vapply(panels, is.null, logical(1))]
  if (length(missing))
    warning("no data for panel(s): ", paste(missing, collapse = ", "), call. = FALSE)
  names(panels) <- PANELS$panel_key
  panels <- Filter(Negate(is.null), panels)

  # the shipped figure's ragged grid: 3 / 2 / 3 / 4 panels per row
  design <- "ABC#\nDE##\nFGH#\nIJKL"
  # a shared legend, taken from a throwaway plot that carries all three tiers
  legend_src <- ggplot(data.table(x = 1:3, y = 1:3, tier = factor(names(TIER_COLS),
                                  levels = names(TIER_COLS))),
                       aes(x, y, colour = tier)) +
    geom_point(size = 2) + scale_colour_manual(values = TIER_COLS, name = NULL) +
    theme_imic(base_size = 8) + theme(legend.position = "bottom")

  # Compose as [panel grid] over [shared legend strip]. Collecting guides across a
  # ragged patchwork design is version-fragile, so the legend is extracted once from
  # `legend_src` and laid out as its own row.
  grid <- wrap_plots(panels, design = design) & theme(legend.position = "none")
  # Shared axis titles: per-panel titles would repeat 12 times, so they are drawn once
  # as text grobs framing the grid (the shipped figure does the same).
  ylab <- wrap_elements(grid::textGrob(expression(-log[10]("Raw P-value")), rot = 90,
                                       gp = grid::gpar(fontsize = 9, fontfamily = "Helvetica")))
  xlab <- wrap_elements(grid::textGrob("Average Treatment Effect",
                                       gp = grid::gpar(fontsize = 9, fontfamily = "Helvetica")))
  body <- (ylab | grid) + plot_layout(widths = c(0.028, 1))
  out  <- body / xlab / extract_legend(legend_src) +
    plot_layout(heights = c(1, 0.022, 0.030)) +
    plot_annotation(theme = theme(plot.margin = margin(4, 6, 2, 4)))

  dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)
  ggsave(out_png, out, width = 7.25, height = 9.5, units = "in", dpi = 300,
         bg = "white", device = ragg::agg_png)
  cat("wrote", out_png, "| panels:", length(panels), "\n")
  cat("  labels drawn per panel (up / down):\n")
  for (k in names(panels)) {
    L <- attr(panels[[k]], "imic_label_split")
    if (!is.null(L))
      cat(sprintf("    %-34s %d / %d\n",
                  gsub("\n", " · ", PANELS$title[match(k, PANELS$panel_key)]), L[1], L[2]))
  }
  invisible(out)
}

# Extract a ggplot's legend as a standalone plot object (no cowplot dependency).
# ggplot2 3.5 renames the guide grob ("guide-box-bottom"), so match on the prefix.
extract_legend <- function(g) {
  gt  <- ggplotGrob(g)
  idx <- which(startsWith(vapply(gt$grobs, function(x) x$name %||% "", character(1)),
                          "guide-box"))
  if (!length(idx)) return(patchwork::plot_spacer())
  keep <- idx[vapply(idx, function(i) !inherits(gt$grobs[[i]], "zeroGrob"), logical(1))]
  if (!length(keep)) return(patchwork::plot_spacer())
  patchwork::wrap_elements(full = gt$grobs[[keep[1]]])
}
`%||%` <- function(a, b) if (is.null(a)) b else a

if (sys.nframe() == 0) invisible(build_figS11())
