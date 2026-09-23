# fig3B-pathway.R
# Fig 3 Panel B = primary KEGG pathway-impact plot.
#
# WEB-INDEPENDENT (2026-08-13): built from a LOCAL MetaboAnalystR KEGG pathway run
# (src/metaboanalyst/run-primary-pathway-local.R -> primary_pathway_local/), NOT the
# downloaded MetaboAnalyst exports. Using the SAME KEGG library the web tool uses
# (SetKEGG.PathLib "metpa", pathlib="kegg") reproduces the submission's IMPACT axis
# EXACTLY (verified 31/31 cells: e.g. NAD/NAM total 15, impact 0.620) -- the earlier
# compression was only because the local run defaulted to the SMPDB library. Raw
# p-values match to a small factor (compound-name DB drift in CrossReferencing);
# pathway identities, hits and impact are identical and significance status is
# preserved (one borderline Thiamine point at impact 0 sits either side of P=0.05).
# The submitted exports are kept at results/metaboanalyst/primary_pathway_trenton/
# for provenance/QC but no longer feed the figure.
#
# Styling follows Trenton's Primary Outcomes (Pathway Analysis).Rmd:
#   x = Impact, y = -log10(Raw p); coloured by study; P<0.05 dashed line;
#   ellipse + labels ONLY for pathways significant in >= 2 study x timepoint cells.
suppressMessages({ library(dplyr); library(ggplot2); library(ggrepel); library(stringr) })
source("figure-scripts/manuscript_figures/study_colors.R")   # canonical study colours (match graphical abstract)
source("figure-scripts/0_figure-functions.R")                # theme_imic() = Science-submission theme

LOCAL_CSV <- "results/metaboanalyst/primary_pathway_local/primary_pathway_all_cells.csv"
OUT_EMBED <- "figures/figure3_panelB_msea.png"
OUT_CMP   <- "Manuscript/figure_comparison/msea_compare/fig3B_replicated_pathway.png"

# Shared 3B plot in the SUBMITTED style so the combined-arm MAIN panel and the
# arm-stratified SUPPLEMENT match exactly. Mirrors the submitted Fig 3B
# (src/trenton-ports/figure-3b-pathway-submitted.R): theme_bw FULL grid, SOLID P<0.05
# (grey) + Q<0.05 (green) reference lines with LEFT-anchored "P-value < ..." text,
# white boxed labels + ellipses, integer y ticks, x breaks at 0.25, and a bottom-LEFT
# legend. Differences from the submission: the study palette (Okabe-Ito since 2026-09-23)
# and a redundant study SHAPE (circle/triangle/square) for colourblind readers.
.plot_3b <- function(tab, ellipse_df, label_df, cols, fdr_p_thr, y_fdr_line) {
  ggplot(tab, aes(impact, logP)) +
    geom_hline(yintercept = -log10(0.05), color = "#BAB0AC", linewidth = 0.4) +   # solid P<0.05
    { if (is.finite(y_fdr_line)) geom_hline(yintercept = y_fdr_line, color = "#59A14F", linewidth = 0.4) } +  # solid Q<0.05
    annotate("text", x = 0, y = -log10(0.05) + 0.18, hjust = 0, size = 2.3, color = "#8A8580",
             parse = TRUE, label = "italic(P)*\"-value\" < 0.05") +
    { if (is.finite(y_fdr_line)) annotate("text", x = 0, y = y_fdr_line + 0.18, hjust = 0, size = 2.3,
             color = "#3C8C3C", parse = TRUE,
             label = paste0("italic(P)*\"-value\" < ", signif(fdr_p_thr, 2))) } +
    { if (nrow(ellipse_df) >= 2) ggforce::geom_mark_ellipse(data = ellipse_df, aes(group = pathway),
             expand = unit(2, "mm"), colour = "black", linewidth = 0.3) } +
    geom_point(aes(color = point_color, shape = point_color), size = 1.6, alpha = 0.9) +
    geom_label_repel(data = label_df, aes(label = lab, color = point_color),
             size = 2.2, fill = "white", label.size = 0.2, box.padding = 0.35,
             point.padding = 0.3, segment.size = 0.35, max.overlaps = Inf,
             show.legend = FALSE, min.segment.length = 0, seed = 123) +
    scale_color_manual(values = cols, name = "Study") +
    scale_shape_manual(values = imic_shapes_for(names(cols)), name = "Study") +   # redundant CVD cue
    scale_x_continuous(breaks = seq(0, 1, 0.25), limits = c(-0.05, 1.05)) +
    scale_y_continuous(breaks = seq(0, ceiling(max(tab$logP, na.rm = TRUE)), 1)) +
    labs(x = "Pathway Impact", y = expression(-log[10](italic(P)*"-value"))) +
    theme_bw(base_size = 8, base_family = "Helvetica") +   # full grid, matching the submission
    # explicit floor overrides (harmonized 2026-08-26): theme_bw()'s sub-elements are
    # RELATIVE to base_size (axis.text/legend.text default to rel(0.8) = 6.4pt at
    # base_size=8), so base_size=8 alone does not actually meet the Reviewer-2 7pt
    # floor the way theme_imic()'s absolute-size overrides do. Kept as theme_bw()
    # rather than switching to theme_imic() because the full grid + border here is a
    # deliberate match to the submitted figure's style, not an oversight.
    theme(panel.grid.minor = element_blank(),
          axis.text        = element_text(size = 8),
          axis.title       = element_text(size = 9),
          legend.text      = element_text(size = 7),
          legend.title     = element_text(size = 8),
          legend.position = "bottom", legend.justification = "left", legend.box = "horizontal") +
    guides(color = guide_legend(override.aes = list(size = 2.5)))   # shape legend merges into this one
}

build <- function() {
  # local KEGG pathway run: already tidy (study, tp, pathway, impact, raw_p, fdr).
  tab <- read.csv(LOCAL_CSV, check.names = FALSE, stringsAsFactors = FALSE) %>%
    transmute(study, tp, pathway, impact = as.numeric(impact),
              raw_p = as.numeric(raw_p), fdr = as.numeric(fdr))

  tab <- tab %>% mutate(
    logP = -log10(raw_p),
    is_sig = raw_p < 0.05,
    point_color = ifelse(is_sig, study, "Not Significant"),
    short = pathway %>%
      str_replace(regex("Nicotinate and nicotinamide metabolism", ignore_case = TRUE), "NAD/NAM Met") %>%
      str_replace(regex("Riboflavin metabolism", ignore_case = TRUE), "Riboflavin Met") %>%
      str_replace(regex("Thiamine metabolism", ignore_case = TRUE), "Thiamine Met") %>%
      str_replace(regex("Vitamin B6 metabolism", ignore_case = TRUE), "Vit B6 Met") %>%
      str_replace(regex("Pantothenate and CoA biosynthesis", ignore_case = TRUE), "Pantothenate Met") %>%
      str_replace("metabolism", "Met"),
    lab = ifelse(is_sig, paste0(short, " (", tp, ")"), NA))

  cols <- c(imic_study_cols, "Not Significant"="grey75")
  rep_df <- tab %>% filter(is_sig) %>% add_count(pathway, name = "pc") %>% filter(pc >= 2)

  # Green FDR reference line = the highest raw p still FDR-significant (pooled across
  # cells), matching the Q line on the MSEA panels (render_msea_panelB.R). FDR is
  # MetaboAnalyst's per-cell BH from the pathway_results.csv exports.
  fdr_p_thr  <- suppressWarnings(max(tab$raw_p[tab$fdr < 0.05], na.rm = TRUE))
  y_fdr_line <- if (is.finite(fdr_p_thr)) -log10(fdr_p_thr) else NA_real_

  # main = combined-arm panel; labels the pathways significant in >= 2 study x time cells.
  p <- .plot_3b(tab, ellipse_df = rep_df, label_df = rep_df, cols = cols,
                fdr_p_thr = fdr_p_thr, y_fdr_line = y_fdr_line)

  # SUBMITTED panel size: 105 x 99 mm (210/2 x 297/3), Trenton's universal panel
  # convention (imicPaperTriglycerides/Proteomics.Rmd) -- near-square (1.06:1). We
  # keep theme_imic fonts but anchor the physical size to the submission.
  for (o in c(OUT_EMBED, OUT_CMP))
    ggsave(o, p, width = 210/2, height = 297/3, units = "mm", dpi = 300, device = ragg::agg_png)
  cat("wrote 3B from LOCAL KEGG pathway run | cells:", n_distinct(paste(tab$study, tab$tp)),
      "| impact range:", round(range(tab$impact, na.rm = TRUE), 3), "| repeated-pathway groups:",
      n_distinct(rep_df$pathway), "\n")
  invisible(p)
}

# ---------------------------------------------------------------------------
# Stratified Panel B (online supplement). Same pathway-impact style as the
# combined panel above, but built from the arm-STRATIFIED scripted pathway run
# (run-primary-pathway-stratified.R -> primary_stratified_pathway_all_cells.csv;
# one cell per study x timepoint x arm-contrast). NOTE: this uses our
# MetaboAnalystR pathway library, so the impact scale is compressed (~0-0.25)
# relative to Trenton's submitted KEGG exports (0-1) used by the combined panel.
# Labels are one per study x pathway (most significant cell) to stay legible
# across the 15 arm-cells; points and ellipses show every significant cell.
# ---------------------------------------------------------------------------
STRAT_CSV <- "results/metaboanalyst/primary_stratified_pathway/primary_stratified_pathway_all_cells.csv"
OUT_STRAT <- "figures/figure3_panelB_stratified.png"

build_stratified <- function() {
  tab <- read.csv(STRAT_CSV, check.names = FALSE, stringsAsFactors = FALSE) %>%
    transmute(
      studyfull = study, contrast = contrast, pathway = pathway,
      impact = as.numeric(impact), raw_p = as.numeric(raw_p), fdr = as.numeric(fdr))

  tab <- tab %>% mutate(
    study = dplyr::case_when(
      grepl("Elicit", studyfull) ~ "ELICIT",
      grepl("Misame", studyfull) ~ "MISAME-III",
      grepl("Vital",  studyfull) ~ "Mumta-LW", TRUE ~ studyfull),
    tp = trimws(gsub("[()]", "", stringr::str_extract(studyfull, "\\(([^)]+)\\)"))),
    logP = -log10(raw_p),
    is_sig = raw_p < 0.05,
    point_color = ifelse(is_sig, study, "Not Significant"),
    short = pathway %>%
      str_replace(regex("Nicotinate and nicotinamide metabolism", ignore_case = TRUE), "NAD/NAM Met") %>%
      str_replace(regex("Riboflavin metabolism", ignore_case = TRUE), "Riboflavin Met") %>%
      str_replace(regex("Thiamine metabolism", ignore_case = TRUE), "Thiamine Met") %>%
      str_replace(regex("Vitamin B6 metabolism", ignore_case = TRUE), "Vit B6 Met") %>%
      str_replace(regex("Pantothenate and CoA biosynthesis", ignore_case = TRUE), "Pantothenate Met") %>%
      str_replace(regex("Biotin metabolism", ignore_case = TRUE), "Biotin Met") %>%
      str_replace(regex("Lactose degradation", ignore_case = TRUE), "Lactose Deg") %>%
      str_replace("metabolism", "Met"))

  cols <- c(imic_study_cols, "Not Significant"="grey75")
  rep_df <- tab %>% filter(is_sig) %>% add_count(pathway, name = "pc") %>% filter(pc >= 2)
  # Ellipses group a pathway's significant cells. In the arm-STRATIFIED view a pathway
  # can land at very different impacts across arms (e.g. Thiamine Metabolism: impact 0.00
  # in one arm-cell, 0.95 in another), which draws a single GIANT ellipse spanning the
  # whole x-axis. Restrict ellipses to pathways whose significant cells actually CLUSTER
  # in impact (span < 0.3) -- keeps the tight NAD/NAM & Vit B6 groups, drops the giant one.
  rep_tight <- rep_df %>% group_by(pathway) %>%
    filter((max(impact) - min(impact)) < 0.3) %>% ungroup()
  # One label per study x pathway (its most significant cell) so 15 arm-cells stay legible.
  lab_df <- tab %>% filter(is_sig) %>% group_by(study, pathway) %>%
    slice_min(raw_p, n = 1, with_ties = FALSE) %>% ungroup() %>% mutate(lab = short)

  fdr_p_thr  <- suppressWarnings(max(tab$raw_p[tab$fdr < 0.05], na.rm = TRUE))
  y_fdr_line <- if (is.finite(fdr_p_thr)) -log10(fdr_p_thr) else NA_real_

  # supplement = arm-stratified; coloured by study (arm/contrast NOT shape-coded -- 7
  # arm x contrast levels across 3 studies would be illegible as shapes). Labels one per
  # study x pathway; ellipses only for impact-clustered pathways (see rep_tight).
  p <- .plot_3b(tab, ellipse_df = rep_tight, label_df = lab_df, cols = cols,
                fdr_p_thr = fdr_p_thr, y_fdr_line = y_fdr_line)

  # SQUARE, matching the main Fig 3B panel (per Andrew).
  ggsave(OUT_STRAT, p, width = 210/2, height = 297/3, units = "mm", dpi = 300, device = ragg::agg_png)
  cat("wrote 3B STRATIFIED | cells:", n_distinct(paste(tab$studyfull, tab$contrast)),
      "| impact range:", round(range(tab$impact, na.rm = TRUE), 3),
      "| repeated-pathway groups:", n_distinct(rep_df$pathway), "\n")
  invisible(p)
}

if (sys.nframe() == 0) { invisible(build()); invisible(build_stratified()) }
