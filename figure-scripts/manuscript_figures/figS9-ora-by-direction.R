# figS9-ora-by-direction.R
# =============================================================================
# SUPPLEMENT : direction-split over-representation (ORA) enrichment volcanoes.
#
# WHAT THIS IS
#   The main-text enrichment panels (see `figure4-panelB-msea.R`) collapse
#   up- and down-regulated features onto ONE signed enrichment-ratio axis.
#   For the supplement, Andrew x Trenton want the ORA of UP-regulated and
#   DOWN-regulated features shown SEPARATELY: instead of signing one axis by
#   direction, we facet by direction. Each facet is a study-coloured volcano:
#     x = enrichment ratio (hits / expected)  -- UNsigned (direction is facet)
#     y = -log10(raw p)
#   with the nominal p = 0.05 line (red dashed), the enrichment-ratio = 2
#   reference (blue dashed), and each significant pathway labelled.
#
# LABELLING RULE (matches the existing panels)
#   Per direction: if any pathway is FDR-significant (fdr_native < 0.05) in that
#   direction, label those; otherwise fall back to the nominally significant
#   ones (raw_p < 0.05). Where a direction has many nominal hits but some FDR
#   hits (e.g. tertiary down-regulated), only the FDR-significant pathways are
#   labelled to keep the panel legible.
#
# INPUT (already computed; NOT recomputed here)
#   results/metaboanalyst/primary_combined/primary_combined_supplementary_table.csv
#   results/metaboanalyst/tertiary_combined/tertiary_combined_supplementary_table.csv
#   (one row per cell x pathway; carries a `direction` (up/down) column.)
#
# OUTPUT
#   figures/figureS_primary_ora_by_direction.png
#   figures/figureS_tertiary_ora_by_direction.png
#   results/metaboanalyst/primary_combined/primary_ora_upregulated.csv
#   results/metaboanalyst/primary_combined/primary_ora_downregulated.csv
#   results/metaboanalyst/tertiary_combined/tertiary_ora_upregulated.csv
#   results/metaboanalyst/tertiary_combined/tertiary_ora_downregulated.csv
#   (each split CSV = the supplementary table filtered to that direction,
#    columns unchanged, written atomically.)
#
# Run from the repo root:
#   Rscript "figure-scripts/manuscript_figures/figS9-ora-by-direction.R"
# =============================================================================

suppressMessages({
  library(dplyr); library(stringr); library(ggplot2); library(ggrepel); library(ggforce)
})
source("figure-scripts/manuscript_figures/study_colors.R")   # canonical study colours (match graphical abstract)

# Tableau-20 palette (Trenton's), one colour per study.
tableau20 <- c(
  "#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
  "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC",
  "#86BCD6", "#FFBE7D", "#FF5850", "#A0CBE8", "#8CD17D",
  "#B6992D", "#499894", "#FABFD2", "#D37295", "#B7B7B7")

ALPHA  <- 0.05   # nominal significance threshold
ER_REF <- 2      # enrichment-ratio reference (Trenton's blue dashed line)

# ---------------------------------------------------------------------------
# Atomic CSV write: write to a sibling tempfile, then rename over the target.
# ---------------------------------------------------------------------------
write_csv_atomic <- function(df, path) {
  tmp <- tempfile(tmpdir = dirname(path), fileext = ".csv.tmp")
  write.csv(df, tmp, row.names = FALSE, na = "")
  if (file.exists(path)) file.remove(path)
  file.rename(tmp, path)
  invisible(path)
}

# ---------------------------------------------------------------------------
# Split the supplementary table by direction and write the four CSVs.
# Columns are left UNCHANGED (the online supplement DT-renders them as-is).
# ---------------------------------------------------------------------------
split_direction_csvs <- function(supp, dir_out, stem) {
  paths <- c()
  for (d in c("up", "down")) {
    sub  <- supp[tolower(supp$direction) == d, , drop = FALSE]
    tag  <- if (d == "up") "upregulated" else "downregulated"
    path <- file.path(dir_out, sprintf("%s_ora_%s.csv", stem, tag))
    write_csv_atomic(sub, path)
    cat(sprintf("  wrote %s  (%d rows)\n", path, nrow(sub)))
    paths <- c(paths, path)
  }
  paths
}

# ---------------------------------------------------------------------------
# Build the two-facet (up | down) direction-split volcano for one outcome set.
# ---------------------------------------------------------------------------
build_direction_figure <- function(supp, out_png, title) {

  set.seed(123)  # reproducible jitter for overlapping non-significant points

  plot_df <- supp %>%
    mutate(
      study_label = recode(study,
                           Elicit = "ELICIT",
                           Misame = "MISAME-III",
                           Vital  = "Mumta-LW",
                           .default = study),
      direction_f = factor(tolower(direction), levels = c("up", "down"),
                           labels = c("Upregulated (ORA)",
                                      "Downregulated (ORA)")),
      enrichment_ratio = hits / expected,
      logP        = -log10(raw_p),
      is_fdr_sig  = !is.na(fdr_native) & fdr_native < ALPHA,
      is_nom_sig  = !is.na(raw_p) & raw_p < ALPHA,
      point_color = if_else(is_nom_sig, study_label, "Not Significant"),
      jitter_x = if_else(is_nom_sig, enrichment_ratio,
                         enrichment_ratio + runif(n(), -0.4, 0.4)),
      jitter_y = if_else(is_nom_sig, logP,
                         pmax(0, logP + runif(n(), -0.12, 0.12)))
    )

  # Per-direction label set: prefer FDR-significant; else nominal. Collapse to
  # one label per pathway x study (drop repeats across visits).
  label_df <- plot_df %>%
    group_by(direction_f) %>%
    group_modify(function(g, ...) {
      base <- if (any(g$is_fdr_sig)) g[g$is_fdr_sig, ] else g[g$is_nom_sig, ]
      base
    }) %>%
    ungroup() %>%
    group_by(direction_f, study_label, pathway) %>%
    slice_max(logP, n = 1, with_ties = FALSE) %>%
    ungroup()

  # Ellipses around "repeated" pathways (Trenton's MSEA/pathway figure feature):
  # within each direction facet, a pathway nominally significant (raw p < 0.05)
  # in >= 2 analysis cells gets a black outline around its significant points.
  repeated_paths <- plot_df %>%
    filter(is_nom_sig) %>%
    count(direction_f, pathway, name = "n_sig_cells") %>%
    filter(n_sig_cells >= 2) %>%
    select(direction_f, pathway)
  ellipse_df <- plot_df %>%
    filter(is_nom_sig) %>%
    semi_join(repeated_paths, by = c("direction_f", "pathway"))
  # geom_mark_ellipse needs >= 2 points per group; guard the empty case.
  ellipse_layer <- if (nrow(ellipse_df) >= 2) {
    ggforce::geom_mark_ellipse(
      data = ellipse_df, aes(x = jitter_x, y = jitter_y, group = pathway),
      color = "black", fill = NA, linewidth = 0.35,
      expand = unit(2, "mm"), inherit.aes = FALSE, show.legend = FALSE)
  } else NULL

  # Colour map: grey for non-significant, CANONICAL colour per study (by name, not
  # by presence order, so a study keeps its colour across panels).
  studies_present <- sort(unique(plot_df$study_label[plot_df$is_nom_sig]))
  color_vals <- c("Not Significant" = "grey75", imic_study_cols[studies_present])

  p <- ggplot(plot_df, aes(x = jitter_x, y = jitter_y)) +
    ellipse_layer +
    geom_point(aes(color = point_color), size = 1.6, alpha = 0.85) +
    geom_hline(yintercept = -log10(ALPHA), linetype = "dashed", color = "red") +
    geom_vline(xintercept = ER_REF, linetype = "dashed", color = "blue") +
    geom_label_repel(
      data = label_df, aes(label = pathway, color = point_color),
      size = 1.9, label.padding = 0.12, box.padding = 0.4,
      min.segment.length = 0, max.overlaps = 200, show.legend = FALSE) +
    scale_color_manual(values = color_vals, name = "Study") +
    facet_wrap(~ direction_f, nrow = 1, scales = "free") +
    labs(x = "Enrichment Ratio (hits / expected)",
         y = expression(-Log[10]*"(Raw P)"),
         title = title) +
    theme_classic(base_size = 7) +
    theme(
      axis.title   = element_text(size = 7),
      axis.text    = element_text(size = 7),
      plot.title   = element_text(size = 8, face = "bold"),
      strip.text   = element_text(size = 7.5, face = "bold"),
      strip.background = element_rect(fill = "white", color = NA),
      legend.title = element_text(size = 7),
      legend.text  = element_text(size = 7),
      legend.position = "bottom",
      legend.key.size = unit(0.3, "cm"),
      plot.margin  = margin(4, 6, 2, 4))

  ggsave(
    filename = out_png,
    plot   = p,
    width  = 190,   # ~A4 text width (mm); two facets side by side
    height = 115,   # mm
    units  = "mm",
    dpi    = 600,
    device = ragg::agg_png)

  cat(sprintf("wrote %s\n", out_png))
  cat(sprintf("  labelled pathways: up=%d, down=%d\n",
              sum(label_df$direction_f == "Upregulated (ORA)"),
              sum(label_df$direction_f == "Downregulated (ORA)")))
  invisible(plot_df)
}

# ===========================================================================
# PRIMARY
# ===========================================================================
cat("== PRIMARY ==\n")
primary <- read.csv(
  "results/metaboanalyst/primary_combined/primary_combined_supplementary_table.csv",
  check.names = FALSE, stringsAsFactors = FALSE)

split_direction_csvs(primary, "results/metaboanalyst/primary_combined", "primary")
build_direction_figure(primary,
                       "figures/figureS_primary_ora_by_direction.png",
                       "Primary ORA by direction (combined arms)")

# ===========================================================================
# TERTIARY
# ===========================================================================
cat("== TERTIARY ==\n")
tertiary <- read.csv(
  "results/metaboanalyst/tertiary_combined/tertiary_combined_supplementary_table.csv",
  check.names = FALSE, stringsAsFactors = FALSE)

split_direction_csvs(tertiary, "results/metaboanalyst/tertiary_combined", "tertiary")
build_direction_figure(tertiary,
                       "figures/figureS_tertiary_ora_by_direction.png",
                       "Tertiary ORA by direction (combined arms)")

# ---------------------------------------------------------------------------
# Significance summary (for the report / caption).
# ---------------------------------------------------------------------------
sig_summary <- function(d, name) {
  fdr <- tapply(d$fdr_native < 0.05, tolower(d$direction), sum, na.rm = TRUE)
  nom <- tapply(d$raw_p       < 0.05, tolower(d$direction), sum, na.rm = TRUE)
  cat(sprintf("%s  FDR-sig up=%s down=%s | nominal up=%s down=%s\n",
              name, fdr["up"], fdr["down"], nom["up"], nom["down"]))
}
cat("\n== significance summary (fdr_native<0.05) ==\n")
sig_summary(primary,  "primary ")
sig_summary(tertiary, "tertiary")
