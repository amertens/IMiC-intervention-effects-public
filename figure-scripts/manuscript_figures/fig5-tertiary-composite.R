# fig5-tertiary-composite.R
# =============================================================================
# Tertiary targeted-lipidome figure (manuscript Fig 5), rebuilt from the
# reproducible results/ + results/metaboanalyst/ outputs.
#
# Mirrors the Figure 3 scripts (src/3 visualizations/figure4-*.R) but for the
# TERTIARY outcome group (targeted Biocrates lipids: triglycerides, ceramides,
# diglycerides, phosphatidylcholines, acylcarnitines, ...). It is deliberately
# SELF-CONTAINED: the volcano-panel machinery and the MSEA-panel logic are
# copied here (not sourced) so this script never edits the Fig 3 originals.
#
# Fig 5 has THREE panels (see the placeholder figures/figure5.png for layout):
#   A) Volcano grid of tertiary lipids, coloured by lipid CATEGORY, with the
#      per-category inset bars and top-feature labels. Combined-arm MAIN
#      (one contrast per study x visit) + arm-stratified SUPPLEMENT.
#   B) Tertiary MSEA: signed enrichment-ratio volcano
#      (results/metaboanalyst/tertiary_msea/tertiary_msea_{combined,stratified}.csv).
#   C) Triglyceride -> fatty-acid composition. That analysis
#      (fig5C-triglyceride.R) is blocked on an external
#      Biocrates structure file, so Panel C here embeds the existing placeholder
#      (figures/final-figure-4-c.png) until figures/figure4_panelC_tg_composition.png
#      is produced.
#
# Composite layout: A on top, (B, C) side by side underneath -- matching the
# A-over-(B,C) arrangement in the figure5.png placeholder.
#
# Outputs:
#   figures/figure5_panelB_msea.png        (tertiary MSEA, standalone)
#   figures/figure5_combined_arms.png      (MAIN: combined-arm A + B + C)
#   figures/figure5_stratified_supplement.png (SUPPLEMENT: stratified A + B + C)
#
# Run from the repo root:
#   Rscript "figure-scripts/manuscript_figures/fig5-tertiary-composite.R"
# =============================================================================

suppressMessages({
  library(dplyr); library(tidyr); library(stringr); library(ggplot2)
  library(ggrepel); library(cowplot); library(magick)
})

# A4 page minus 0.5in margins (same geometry as the Fig 3 composite).
PAGE_W <- 8.27 - 2 * 0.5   # 7.27 in
PAGE_H <- 11.69 - 2 * 0.5  # 10.69 in

tableau10 <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728",
               "#9467BD", "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF")

# Assigned once, after the data / category set is known; read by the panel,
# inset and legend builders below.
final_color_palette <- NULL

# X-axis limits for the volcano panels. Assigned once from the tertiary effect
# range so every panel shares one scale (breaks fixed at -1/0/1).
VOLCANO_XLIM <- NULL

# render_msea_panelB.R (Panel B) is sourced HERE, BEFORE this file's local
# plot_imic_volcano_panel definition below. render_msea_panelB.R pulls in
# 0_figure-functions.R (for theme_imic), which also carries a stale
# plot_imic_volcano_panel; sourcing it first lets our local theme_imic version win.
source("figure-scripts/manuscript_figures/render_msea_panelB.R")

# ===========================================================================
# PANEL A -- volcano machinery (copied from figure4-primary-volcano.R, then
# adapted for the tertiary lipidome: labels use the `biomarker` short code
# rather than `label_f`, and the x-scale is widened to the tertiary effect
# range). Behaviour otherwise preserved.
# ===========================================================================

png_to_ggdraw <- function(path) ggdraw() + draw_image(image_read(path))

# Shorten a feature code for the in-panel repel labels: drop any parenthetical
# qualifier and cap the length so long names do not spill past the panel edge.
short_label <- function(x) {
  x <- sub("\\s*\\(.*$", "", x)
  x <- trimws(x)
  ifelse(nchar(x) > 22, paste0(substr(x, 1, 21), "…"), x)
}

# BH critical raw-p: the largest raw p whose BH q is still <= alpha.
get_bh_cutoff <- function(df, p_col = "pval", q_col = "qval", alpha = 0.05,
                          return = c("log10", "raw", "both")) {
  return <- match.arg(return)
  if (!all(c(p_col, q_col) %in% names(df)))
    stop("Specified p_col / q_col not found in the data frame.")
  is_sig <- df[[q_col]] <= alpha & !is.na(df[[q_col]]) & !is.na(df[[p_col]])
  if (!any(is_sig)) { warning("No q-values <= alpha; returning NA."); p_crit <- NA_real_ }
  else               p_crit <- max(df[[p_col]][is_sig])
  switch(return,
         log10 = -log10(p_crit),
         raw   =  p_crit,
         both  =  list(raw = p_crit, log10 = -log10(p_crit)))
}

# Small stacked-bar inset: per chemical category, total features (faint) vs
# FDR-significant features (solid), drawn top-right of each volcano panel.
make_inset_with_labels <- function(df_sum, df_long2) {
  xmax_val <- max(df_sum$n_total) * 1.05
  p <- ggplot(df_long2) +
    geom_rect(aes(xmin = 0, xmax = count, ymin = ymin, ymax = ymax,
                  fill = category, alpha = type), color = "black", size = 0.15) +
    scale_fill_manual(values = final_color_palette) +
    scale_alpha_manual(values = c(n_total = 0.15, n_sig = 1)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, xmax_val)) +
    coord_flip(clip = "off") +
    theme_void(base_size = 6) +
    theme(panel.background = element_rect(fill = "white", color = NA),
          legend.position = "none", plot.margin = margin(0, 0, 0, 0))
  p + annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
               fill = NA, color = "black", size = 0.1)
}

# One volcano panel (scaled ATE vs -log10 p) for a single study/visit/contrast,
# with the category inset and repel-labelled top FDR-significant features.
plot_imic_volcano_panel <- function(res, title = "", n_top_vars = 5, overlap_n = 20) {
  q_cut <- get_bh_cutoff(res, p_col = "pval", q_col = "pval_adj",
                         alpha = 0.05, return = "raw")

  tt_volcano <- res %>%
    arrange(pval_adj) %>%
    mutate(
      ATE        = est,
      logPval    = -log10(pval),
      sig_status = case_when(
        pval <= q_cut ~ "Significant after FDR",
        pval < 0.05   ~ "Significant before FDR",
        TRUE          ~ "Not Significant"),
      color_var  = ifelse(sig_status == "Significant after FDR",
                          as.character(category), sig_status),
      # tertiary panels label with the short biomarker code (e.g. "Tg.18.3_30.0.")
      # to match the placeholder; label_f collapses every TG to "Triacylglyceride".
      lab_src    = ifelse(sig_status == "Significant after FDR", biomarker, ""))

  category_summary <- tt_volcano %>%
    group_by(category) %>%
    summarise(n_total = n(),
              n_sig = sum(sig_status == "Significant after FDR"),
              rate_sig = n_sig / n_total, .groups = "drop") %>%
    mutate(share_sig = n_sig / sum(n_sig))

  df_sum <- category_summary %>%
    mutate(pct_within = n_sig / n_total,
           pct_total  = n_total / sum(n_total),
           cat_num    = as.numeric(factor(category)),
           ymin       = cat_num - 0.5,
           ymax       = cat_num + 0.5,
           non_sig    = n_total - n_sig)

  df_long2 <- df_sum %>%
    select(category, cat_num, ymin, ymax, n_sig, n_total) %>%
    pivot_longer(cols = c(n_sig, n_total), names_to = "type", values_to = "count")

  inset <- make_inset_with_labels(df_sum, df_long2)

  p_line <- -log10(0.05)   # nominal 0.05 reference line
  q_line <- -log10(q_cut)  # BH-FDR reference line

  ymax_data <- max(tt_volcano$logPval[is.finite(tt_volcano$logPval)], na.rm = TRUE)

  p <- ggplot(tt_volcano, aes(x = ATE, y = logPval)) +
    geom_point(aes(colour = color_var, shape = sig_status), size = 1, alpha = 0.75) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = p_line, linetype = "dashed", colour = tableau10[2]) +
    geom_hline(yintercept = q_line, linetype = "dotted", colour = tableau10[3]) +
    xlab("") + ylab("") + ggtitle(title) +
    scale_color_manual(values = final_color_palette, na.value = "#999999") +
    guides(color = guide_legend(title = NULL)) +
    # x fixed to the shared tertiary effect range (breaks at -1/0/1); y gets
    # ~12% headroom so top-feature labels are not jammed against the ceiling.
    scale_x_continuous(breaks = c(-1, 0, 1), limits = VOLCANO_XLIM) +
    scale_y_continuous(breaks = seq(0, ceiling(ymax_data), by = 2),
                       limits = c(0, ymax_data * 1.12)) +
    theme_imic(base_size = 8) +   # Science-submission theme (Helvetica, font floors)
    theme(legend.position = "none",
          # keep a light border so the 8 volcano facets stay delineated (theme_imic drops it)
          panel.border = element_rect(colour = "grey75", fill = NA, linewidth = 0.3),
          # left-justified so the 3-line facet title sits top-LEFT, clear of the
          # top-right inset bar chart (they collided when centered + bold).
          plot.title = element_text(size = 8, face = "bold", hjust = 0, lineheight = 0.9,
                                    margin = margin(b = 1)),
          plot.margin = margin(t = 1, r = 3, b = -4, l = 0),
          panel.spacing = unit(0, "pt"))

  top_vars <- tt_volcano %>% filter(pval_adj < q_cut) %>%
    arrange(-logPval) %>% head(n = n_top_vars) %>%
    mutate(lab = short_label(lab_src))
  # ylim caps repelled labels below the top ~15% of the panel's own plot area -- that
  # band is reserved for the top-right category-histogram inset (drawn separately via
  # draw_plot() below, so ggrepel has no way to know about it / avoid it on its own).
  # Without this, a label attached to a near-ceiling point can get pushed up into the
  # inset's footprint in panels with a smaller y-range (reported 2026-08-26: e.g.
  # "Dg.18.1_18.3." in MISAME-III 1-2m collided with the inset even though the
  # tallest-range panel, MISAME-III 14-21d, looked fine).
  p <- p + geom_text_repel(data = top_vars, aes(label = lab),
                           max.overlaps = getOption("ggrepel.max.overlaps", default = overlap_n),
                           size = 2.5, alpha = 0.5,
                           ylim = c(NA, ymax_data * 1.12 * 0.85))

  ggdraw() + draw_plot(p, 0, 0, 1, 1) +
    # Inset height 0.16 (overlapped plotted points) -> 0.10 (left a gap below the
    # panel's title band, since y=1 in this [0,1] canvas is the top of the whole
    # grob including the 2-line title, not the panel border) -> 0.13 (still left a
    # visible gap BELOW the inset, above the panel's top border -- growing height
    # upward from a fixed bottom never touches that gap) -> 0.20 with bottom lowered
    # to y = 0.79 (closed the gap, but then collided with the top-right data label in
    # the tallest-range facet, e.g. "Tg.18.3_30.0.") -> y = 0.82 / height = 0.17
    # (reported still too tall) -> y = 0.87 / height = 0.12 (reported STILL colliding
    # with a repelled label in a DIFFERENT, smaller-range facet, MISAME-III 1-2m --
    # a fixed local-canvas position can't account for where ggrepel happens to place a
    # label in any given panel). Real fix: geom_text_repel's ylim above now keeps every
    # label out of the inset's footprint regardless of panel, so the inset position
    # itself just needs to look right -- y = 0.88 / height = 0.11, close to the panel border.
    draw_plot(inset, x = 0.72, y = 0.83, width = 0.265, height = 0.16)
}

# Categories rarer than `threshold` of features are pooled into "Other".
collapse_small_categories <- function(x, threshold = 0.04) {
  x <- as.character(x)
  tab <- prop.table(table(x))
  small <- names(tab[tab < threshold])
  factor(ifelse(x %in% small, "Other", x))
}

# Category -> colour map: the two "not/before FDR" states are black; each real
# category gets a tableau10 colour (rainbow fallback beyond 10 categories).
build_palette <- function(categories) {
  categories <- unique(categories[!is.na(categories)])
  # Triglycerides FIRST so it gets tableau10[1] = blue (matches the submitted
  # figure, whose legend lists Triglycerides against the blue dot); the rest are
  # sorted so the mapping is STABLE across re-runs (previously a plain alphabetical
  # sort pushed Triglycerides to brown).
  ordered <- c(intersect("Triglycerides", categories),
               sort(setdiff(categories, "Triglycerides")))
  fixed <- c("Not Significant" = "#000000", "Significant before FDR" = "#000000")
  cat_colors <- tableau10[seq_len(min(length(ordered), length(tableau10)))]
  names(cat_colors) <- ordered[seq_along(cat_colors)]
  if (length(ordered) > length(cat_colors)) {
    extra <- grDevices::rainbow(length(ordered) - length(cat_colors))
    names(extra) <- ordered[(length(cat_colors) + 1):length(ordered)]
    cat_colors <- c(cat_colors, extra)
  }
  c(fixed, cat_colors)
}

# Standalone colour legend (categories only) shown in the last grid cell.
create_category_legend <- function() {
  legend_colors <- final_color_palette[
    !names(final_color_palette) %in% c("Not Significant", "Significant before FDR")]
  legend_data <- data.frame(category = names(legend_colors),
                            y = seq_along(legend_colors), x = 1)
  ggplot(legend_data, aes(x = x, y = y, fill = category)) +
    geom_point(size = 3, shape = 21, color = "black") +
    scale_fill_manual(values = legend_colors) +
    geom_text(aes(label = category), hjust = 0, nudge_x = 0.2, size = 2.5) +
    xlim(0.8, 3) + ylim(0.5, length(legend_colors)) +
    theme_void() + theme(legend.position = "none") +
    ggtitle("Significant Categories") +
    theme(plot.title = element_text(size = 8, hjust = 0.5))
}

# ===========================================================================
# PANEL B -- tertiary MSEA signed-enrichment-ratio volcano (adapted from
# figure4-panelB-msea.R). Input columns differ from the primary table:
#   study, timepoint, contrast, direction, pathway, total, expected, hits,
#   raw_p, fdr_native, enrichment_ratio
# where `enrichment_ratio` is ALREADY signed by direction (negative = down),
# so no hits/expected recomputation is needed. No `significant`/FDR-flag column
# exists and (as of this data) NO pathway is FDR-significant, so -- exactly as
# the Fig 3 panel-B template does -- we label the NOMINALLY significant
# (raw_p < 0.05) pathways and colour points by study.
# ===========================================================================
# (render_msea_panelB.R is sourced near the top so our local plot_imic_volcano_panel wins.)
# Pathway-name abbreviations matching the submitted 5B labels.
abbr_5b <- function(x) {
  x <- dplyr::recode(x,
    "Spermidine and Spermine Biosynthesis" = "Spermidine & Spermine Syn",
    "Arginine and Proline Metabolism" = "Arg & Pro Met", "Glycine and Serine Metabolism" = "Gly & Ser Met",
    "Methionine Metabolism" = "Met Met", "Glutamate Metabolism" = "Glu Met",
    "Aspartate Metabolism" = "Asp Met", "Alanine Metabolism" = "Ala Met",
    "Pyruvate Metabolism" = "Pyruvate Met", "Glutathione Metabolism" = "GSH Met",
    "Ammonia Recycling" = "NH3 Rec", "Warburg Effect" = "Warburg",
    "Carnitine Synthesis" = "Carn Synthesis", "Homocysteine Degradation" = "Hcy Deg",
    "Oxidation of Branched Chain Fatty Acids" = "Ox of BCFAs",
    "Beta Oxidation of Very Long Chain Fatty Acids" = "β-Ox VLCFAs",
    "Bile Acid Biosynthesis" = "Bile Acid Biosyn", "Ketone Body Metabolism" = "Ketone Body Met",
    "Mitochondrial Electron Transport Chain" = "Mito ETC",
    "Phytanic Acid Peroxisomal Oxidation" = "Phytanic Ox", "Butyrate Metabolism" = "Butyrate Met",
    "Purine Metabolism" = "Purine Met", .default = x)
  gsub(" Metabolism", " Met", x)
}
build_panelB <- function(msea_csv, out_png, show_legend = TRUE) {
  # Delegates to the shared renderer so Fig 5B and Fig 6A are drawn by ONE code
  # path with identical conventions (size floor >= 3, FDR labelling, P/Q lines).
  # The original in-script implementation below is now UNREACHABLE (kept only for
  # reference); render_msea_panelB.R is the single source of truth.
  # Submitted 5B style: colour ALL nominally-significant pathways by study (no
  # "Sig before FDR" open tier), P<0.05 grey + Q<0.05 green lines, black vline at 0,
  # label every nominally-significant pathway.
  # min_size = 1: the submitted 5B (from Trenton's downloaded tables) applied no set-
  # size floor, so 1-2 member lipid pathways (Ketone Body, Phytanic Ox, ...) are shown.
  # theme_imic (Helvetica, base_size 9); Fig 5 panel B, sized as a half-page sub-panel
  # of the full-page (7.25 in) Fig 5 composite. label_max trimmed for legibility at the
  # larger submission font sizes.
  # Restrict the repel labels to ONLY the pathways discussed in the Results (every
  # point is still plotted; only the TEXT LABELS are limited). Exact data-column names
  # (note "Mitochondrial Electron Transport Chain" for the 3-4-month up-reversal);
  # this drops the previously-labelled-but-undiscussed Spermidine & Spermine Biosynthesis,
  # Valine/Leucine/Isoleucine Degradation, Propanoate Metabolism, and Tryptophan Metabolism.
  label_allow_5b <- c(
    "Methionine Metabolism", "Glutamate Metabolism", "Glycine and Serine Metabolism",
    "Arginine and Proline Metabolism", "Fatty Acid Biosynthesis", "Bile Acid Biosynthesis",
    "Ammonia Recycling", "Glutathione Metabolism", "Carnitine Synthesis",
    "Mitochondrial Electron Transport Chain", "Ketone Body Metabolism",
    "Phytanic Acid Peroxisomal Oxidation")
  return(render_msea_panelB(msea_csv, out_png, title = NULL, min_size = 1,
                            submitted_style = TRUE, upper_line = "fdr",
                            label_which = "nominal", label_allow = label_allow_5b, vline0 = TRUE,
                            abbr_fun = abbr_5b, xlab = "Enrichment ratio",
                            width_in = 210/25.4/2, height_in = 297/25.4/3,
                            show_legend = show_legend))
  supp <- read.csv(msea_csv, check.names = FALSE, stringsAsFactors = FALSE)

  tableau20 <- c(
    "#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
    "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC",
    "#86BCD6", "#FFBE7D", "#FF5850", "#A0CBE8", "#8CD17D",
    "#B6992D", "#499894", "#FABFD2", "#D37295", "#B7B7B7")

  ALPHA  <- 0.05   # nominal significance threshold
  ER_REF <- 2      # enrichment-ratio reference (blue dashed line)

  set.seed(123)    # reproducible jitter for overlapping non-significant points

  plot_df <- supp %>%
    mutate(
      # study is already recoded in these CSVs; recode() leaves matches unchanged.
      study_label = recode(study,
                           Elicit = "ELICIT",
                           Misame = "MISAME-III",
                           Vital  = "Mumta-LW",
                           .default = study),
      enrichment_signed = enrichment_ratio,   # already signed by direction
      logP        = -log10(raw_p),
      is_sig      = raw_p < ALPHA,
      point_color = if_else(is_sig, study_label, "Not Significant"),
      jitter_x = if_else(is_sig, enrichment_signed,
                         enrichment_signed + runif(n(), -0.6, 0.6)),
      jitter_y = if_else(is_sig, logP,
                         pmax(0, logP + runif(n(), -0.15, 0.15)))
    )

  studies_present <- sort(unique(plot_df$study_label[plot_df$is_sig]))
  color_vals <- setNames(tableau20[seq_along(studies_present)], studies_present)
  color_vals <- c("Not Significant" = "grey75", color_vals)

  x_hi <- max(c(plot_df$jitter_x, ER_REF + 1), na.rm = TRUE) * 1.05
  x_lo <- min(c(plot_df$jitter_x, 0), na.rm = TRUE)
  x_lo <- if (x_lo < 0) x_lo * 1.1 - 0.5 else -0.5
  y_hi <- max(plot_df$logP, na.rm = TRUE) * 1.12

  # One label per pathway x study (a pathway can recur across visits). Label only
  # FDR-significant pathways (matches the submitted panel's key set); nominally-
  # significant-but-not-FDR points stay coloured but unlabelled for legibility.
  label_df <- plot_df %>%
    filter(fdr_native < ALPHA, hits >= 1) %>%
    group_by(study_label, pathway) %>%
    slice_max(logP, n = 1, with_ties = FALSE) %>%
    ungroup()

  panelB <- ggplot(plot_df, aes(x = jitter_x, y = jitter_y)) +
    geom_point(aes(color = point_color), size = 1.6, alpha = 0.85) +
    geom_hline(yintercept = -log10(ALPHA), linetype = "dashed", color = "red") +
    geom_vline(xintercept = ER_REF, linetype = "dashed", color = "blue") +
    geom_label_repel(
      data = label_df, aes(label = pathway, color = point_color),
      size = 2, label.padding = 0.12, box.padding = 0.4,
      min.segment.length = 0, max.overlaps = 200, show.legend = FALSE) +
    scale_color_manual(values = color_vals, name = "Study") +
    scale_x_continuous(limits = c(x_lo, x_hi)) +
    scale_y_continuous(limits = c(0, y_hi)) +
    labs(x = "Enrichment Ratio (signed by direction)",
         y = expression(-Log[10]*"(Raw P)")) +
    theme_classic(base_size = 7) +
    theme(
      axis.title   = element_text(size = 7),
      axis.text    = element_text(size = 7),
      legend.title = element_text(size = 7),
      legend.text  = element_text(size = 7),
      legend.position = "bottom",
      legend.key.size = unit(0.3, "cm"),
      plot.margin  = margin(4, 6, 2, 4))

  ggsave(filename = out_png, plot = panelB,
         width = 210 / 2, height = 297 / 3, units = "mm",
         dpi = 600, device = ragg::agg_png)

  cat("wrote", out_png, "\n")
  cat("  pathways plotted:", nrow(plot_df),
      "| nominally significant (raw p <", ALPHA, "):", sum(plot_df$is_sig),
      "| FDR-significant:", sum(plot_df$fdr_native < ALPHA, na.rm = TRUE), "\n")
  cat("  studies with significant pathways:",
      paste(studies_present, collapse = ", "), "\n")
  invisible(panelB)
}

# One shared "Study" legend for Panels B+C (both coloured by the same canonical
# imic_study_cols map, plus "Not Significant"), used in place of each panel's
# own per-panel legend so the composite carries just one.
create_study_legend <- function() {
  study_cols_all <- c(imic_study_cols, "Not Significant" = "grey75")
  legend_df <- data.frame(x = 1, y = seq_along(study_cols_all),
                          grp = factor(names(study_cols_all), levels = names(study_cols_all)))
  p <- ggplot(legend_df, aes(x, y, color = grp)) +
    geom_point(size = 2) +
    scale_color_manual(values = study_cols_all, name = "Study") +
    guides(colour = guide_legend(nrow = 1)) +
    theme_void(base_size = 9) +
    theme(legend.position = "bottom", legend.key.size = unit(0.35, "cm"))
  # cowplot::get_legend() can grab an empty guide-box when a theme_void() plot
  # has several (mostly-empty) guide-box slots -- pull every guide-box and keep
  # the one that actually has content instead of trusting "the first one".
  comps <- cowplot::get_plot_component(p, "guide-box", return_all = TRUE)
  good  <- Filter(function(g) !inherits(g, "zeroGrob"), comps)
  if (length(good) == 0) stop("create_study_legend(): no non-empty guide-box found")
  good[[1]]
}

# ===========================================================================
# COMPOSE -- Panel A grid on top, Panels B & C side by side underneath, one
# shared Study legend beneath B+C.
# ===========================================================================
compose_figure <- function(volcano_grid, panel_b_png, panel_c_png,
                            rel_a = 2.1) {
  grid_labeled <- ggdraw(volcano_grid) +
    draw_label("Scaled Average Treatment Effect", x = 0.5, y = 0.0025,
               hjust = 0.5, vjust = 0, size = 8) +
    draw_label("-log10(P-value)", x = 0.0025, y = 0.5, angle = 90,
               hjust = 0.5, vjust = 1, size = 8)

  bottom_row <- plot_grid(
    png_to_ggdraw(panel_b_png), png_to_ggdraw(panel_c_png),
    labels = c("B", "C"), ncol = 2, nrow = 1, rel_widths = c(1, 1))

  bottom_with_legend <- plot_grid(bottom_row, create_study_legend(),
                                  ncol = 1, rel_heights = c(1, 0.06))

  plot_grid(grid_labeled, bottom_with_legend,
            labels = c("A", ""), ncol = 1, nrow = 2,
            rel_heights = c(rel_a, 1))
}

save_figure <- function(fig, path) {
  ggsave(filename = path, plot = fig, width = PAGE_W, height = PAGE_H,
         units = "in", dpi = 300, device = ragg::agg_png, bg = "white")
  cat("wrote", path, "\n")
}

# ===========================================================================
# DATA (tertiary targeted lipidome)
# ===========================================================================
combined_arms   <- readRDS("results/combined_intervention_effects_results_combined_arms.RDS") %>%
  filter(outcome_group == "tertiary", measure == "ATE")
stratified_arms <- readRDS("results/combined_intervention_effects_results_stratified_arms.RDS") %>%
  filter(outcome_group == "tertiary", measure == "ATE")

# MAIN: every study in its combined-arm framing (one contrast per study x visit).
res_combined <- combined_arms
# SUPPLEMENT: Misame/Vital stratified per arm; Elicit uses its combined-arm rows.
res_stratified <- bind_rows(stratified_arms %>% filter(study != "Elicit"),
                            combined_arms  %>% filter(study == "Elicit"))

res_combined$category   <- collapse_small_categories(res_combined$category)
res_stratified$category <- collapse_small_categories(res_stratified$category)

# One shared palette across both figures so a category is the same colour in the
# main and supplement versions.
final_color_palette <- build_palette(c(as.character(res_combined$category),
                                       as.character(res_stratified$category)))

# One shared x-scale for every volcano panel, from the tertiary effect range.
.est_all <- c(res_combined$est, res_stratified$est)
VOLCANO_XLIM <- c(floor(min(.est_all, na.rm = TRUE) * 10) / 10 - 0.1,
                  ceiling(max(.est_all, na.rm = TRUE) * 10) / 10 + 0.15)

category_legend <- create_category_legend()
blank_plot <- ggplot() + theme_void()

panel <- function(res, study, visit, contrast, title)
  plot_imic_volcano_panel(res %>% filter(study == !!study, visit == !!visit,
                                         contrast == !!contrast), title = title)

# ===========================================================================
# PANEL B (standalone) -- rendered once, embedded in both composites.
# ===========================================================================
# Fig 5B source (2026-08-13, web-independence): the LOCAL dual runner
# (src/metaboanalyst/run-tertiary-msea-dual.R) reproduces the submitted panel from
# raw ATE with no web tool -- metabolite + LIPID-MAPS-converted lipid ORA passes
# against the QER reference metabolome, combined+stratified arms. Recovers 16/20
# submitted pathways (+4 borderline). Was: tertiary_msea_fromTables.csv (built from
# hand-downloaded MetaboAnalyst web ORA tables in data/msea/).
build_panelB("results/metaboanalyst/tertiary_msea/tertiary_msea_dual.csv",
             "figures/figure5_panelB_msea.png", show_legend = FALSE)

# Panel C: prefer the no-legend twin of the real triglyceride->fatty-acid
# composition volcano (compose_figure() below draws ONE shared Study legend for
# B+C instead); fall back to the legend-bearing version, then the committed
# placeholder (the Fig 5C analysis is blocked on the external Biocrates
# structure file).
panel_c_png <- if (file.exists("figures/figure4_panelC_tg_composition_nolegend.png")) {
  "figures/figure4_panelC_tg_composition_nolegend.png"
} else if (file.exists("figures/figure4_panelC_tg_composition.png")) {
  "figures/figure4_panelC_tg_composition.png"
} else {
  "figures/final-figure-4-c.png"
}

# ===========================================================================
# MAIN figure -- combined arms (7 volcano panels + legend, 4x2)
# ===========================================================================
# Facet titles: study + visit on ONE visible row, no arm label (all combined-arm here).
# Two LEADING blank rows ("\n\n<title>") put the facet title on the BOTTOM row of the
# 3-line title band, so the white space is ABOVE the title (under the top-of-panel
# histogram) rather than below it; the band height is unchanged.
c_e1 <- panel(res_combined, "Elicit", "1 mo.",      "Nico", "\n\nELICIT 1m")
c_m1 <- panel(res_combined, "Misame", "14-21 days", "BEP",  "\n\nMISAME-III 14-21d")
c_v1 <- panel(res_combined, "Vital",  "1.5 mo.",    "BEP",  "\n\nMumta-LW 1.5m")
c_m2 <- panel(res_combined, "Misame", "1-2 mo.",    "BEP",  "\n\nMISAME-III 1-2m")
c_v2 <- panel(res_combined, "Vital",  "2 mo.",      "BEP",  "\n\nMumta-LW 2m")
c_e2 <- panel(res_combined, "Elicit", "5 mo.",      "Nico", "\n\nELICIT 5m")
c_m3 <- panel(res_combined, "Misame", "3-4 mo.",    "BEP",  "\n\nMISAME-III 3-4m")

# Study-major layout: one study per row, ordered MISAME-III -> Mumta-LW -> ELICIT.
# Legend on the 2nd row (Mumta has only 2 visits, so its 3rd cell holds the legend);
# the empty cell falls to the 3rd row instead.
combined_grid <- plot_grid(
  c_m1, c_m2, c_m3,                 # row 1: MISAME-III: 14-21d, 1-2m, 3-4m
  c_v1, c_v2, category_legend,      # row 2: Mumta-LW: 1.5m, 2m + legend
  c_e1, c_e2, blank_plot,           # row 3: ELICIT: 1m, 5m
  ncol = 3, nrow = 3, align = "hv", axis = "lrtb")

fig5_main <- compose_figure(combined_grid, "figures/figure5_panelB_msea.png",
                            panel_c_png, rel_a = 2.1)
# 3-way export (PDF + EPS + PNG) via the shared save_figure_3way() (0_figure-functions.R):
# name="figure5" writes figures/figure5.{pdf,eps,png} directly. Panel A is real ggplot
# vector; panels B/C are embedded raster (draw_image), which cairo keeps as raster while
# leaving A vector -> reasonable file size. The .qmd embeds figure5.png.
save_figure_3way(fig5_main, "figure5", width = PAGE_W, height = PAGE_H)
# keep the _combined_arms alias in sync (some docs reference it).
file.copy("figures/figure5.png", "figures/figure5_combined_arms.png", overwrite = TRUE)
cat("wrote figures/figure5.{png,pdf,eps} (+ figure5_combined_arms.png alias)\n")

# ===========================================================================
# SUPPLEMENT figure -- arm-stratified (6x3 layout, 15 volcano panels + legend)
# ===========================================================================
s1  <- panel(res_stratified, "Elicit", "1 mo.",      "Nico",         "ELICIT\n1m\nNico.")
s2  <- panel(res_stratified, "Misame", "14-21 days", "BEP/BEP",      "MISAME-III\n14-21d\nBEP/BEP")
s3  <- panel(res_stratified, "Misame", "14-21 days", "BEP/IFA",      "MISAME-III\n14-21d\nBEP/IFA")
s4  <- panel(res_stratified, "Misame", "14-21 days", "IFA/BEP",      "MISAME-III\n14-21d\nIFA/BEP")
s5  <- panel(res_stratified, "Vital",  "1.5 mo.",    "BEP+ExBf",     "Mumta-LW\n1.5m\nBEP")
s6  <- panel(res_stratified, "Vital",  "1.5 mo.",    "BEP+ExBf+AZT", "Mumta-LW\n1.5m\nBEP+AZT")
s7  <- panel(res_stratified, "Misame", "1-2 mo.",    "BEP/BEP",      "MISAME-III\n1-2m\nBEP/BEP")
s8  <- panel(res_stratified, "Misame", "1-2 mo.",    "BEP/IFA",      "MISAME-III\n1-2m\nBEP/IFA")
s9  <- panel(res_stratified, "Misame", "1-2 mo.",    "IFA/BEP",      "MISAME-III\n1-2m\nIFA/BEP")
s10 <- panel(res_stratified, "Vital",  "2 mo.",      "BEP+ExBf",     "Mumta-LW\n2m\nBEP")
s11 <- panel(res_stratified, "Vital",  "2 mo.",      "BEP+ExBf+AZT", "Mumta-LW\n2m\nBEP+AZT")
s12 <- panel(res_stratified, "Elicit", "5 mo.",      "Nico",         "ELICIT\n5m\nNico.")
s13 <- panel(res_stratified, "Misame", "3-4 mo.",    "BEP/BEP",      "MISAME-III\n3-4m\nBEP/BEP")
s14 <- panel(res_stratified, "Misame", "3-4 mo.",    "BEP/IFA",      "MISAME-III\n3-4m\nBEP/IFA")
s15 <- panel(res_stratified, "Misame", "3-4 mo.",    "IFA/BEP",      "MISAME-III\n3-4m\nIFA/BEP")

stratified_grid <- plot_grid(
  s1,  s2,  s3,  s4,  s5,  s6,
  blank_plot, s7, s8, s9, s10, s11,
  s12, s13, s14, s15, blank_plot, category_legend,
  ncol = 6, nrow = 3, align = "hv", axis = "lrtb")

save_figure(
  compose_figure(stratified_grid, "figures/figure5_panelB_msea.png", panel_c_png,
                 rel_a = 2.1),
  "figures/figure5_stratified_supplement.png")
