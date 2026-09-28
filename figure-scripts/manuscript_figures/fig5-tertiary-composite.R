# =============================================================================
# fig5-tertiary-composite.R
#
# Builds Fig 5 (tertiary outcomes, the targeted Biocrates metabolites and lipids).
# Panel A: volcano grid of the combined-arm intervention effects, one panel per
# study x visit; FDR-significant features are coloured and shaped by chemical
# category, the top 5 are labelled, and an inset bar chart shows the per-category
# share of features tested (light) and FDR-significant (dark). Panel B: MSEA
# signed enrichment-ratio volcano, drawn here by render_msea_panelB(). Panel C:
# triglyceride fatty-acid composition, the PNG written by fig5C-triglyceride.R.
# Layout: A on top, B and C side by side underneath with one shared study legend.
#
# Inputs:  results/combined_intervention_effects_results_{combined,stratified}_arms.RDS
#            (src/2 analysis/clean_results.R; the stratified file is used only for the
#            shared category palette and x-range, so the panels match the submitted figure)
#          results/metaboanalyst/tertiary_msea/tertiary_msea_dual.csv
#            (src/metaboanalyst/run-tertiary-msea-dual.R; also Table S3)
#          figures/figure5_panelC_tg_composition_nolegend.png (fig5C-triglyceride.R)
# Outputs: figures/figure5_panelB_msea.png, figures/figure5.{png,pdf,eps}
# [needs an on-request file] the combined-results RDS files are feature-level result
# files that are not shipped (size). Panel C also needs the Biocrates structure file
# (see fig5C-triglyceride.R).
# Run from the repo root:
#   Rscript "figure-scripts/manuscript_figures/fig5-tertiary-composite.R"
# =============================================================================

suppressMessages({
  library(dplyr); library(tidyr); library(stringr); library(ggplot2)
  library(ggrepel); library(cowplot); library(magick)
})

# A4 page minus 0.5 in margins.
PAGE_W <- 8.27 - 2 * 0.5   # 7.27 in
PAGE_H <- 11.69 - 2 * 0.5  # 10.69 in

# Assigned once, after the data / category set is known; read by the panel,
# inset and legend builders below.
final_color_palette <- NULL

# X-axis limits for the volcano panels. Assigned once from the tertiary effect
# range so every panel shares one scale (breaks fixed at -1/0/1).
VOLCANO_XLIM <- NULL

# render_msea_panelB.R (Panel B) also pulls in 0_figure-functions.R (theme_imic(),
# get_bh_cutoff(), imic_logp_title) and study_colors.R.
source("figure-scripts/manuscript_figures/render_msea_panelB.R")

# ===========================================================================
# PANEL A -- volcano panels, as in Fig 3A but for the tertiary features: labels
# use the `biomarker` code rather than `label_f`, and the x-scale is widened to
# the tertiary effect range.
# ===========================================================================

png_to_ggdraw <- function(path) ggdraw() + draw_image(image_read(path))

# Shorten a feature code for the in-panel repel labels: drop any parenthetical
# qualifier and cap the length so long names do not spill past the panel edge.
short_label <- function(x) {
  x <- sub("\\s*\\(.*$", "", x)
  x <- trimws(x)
  ifelse(nchar(x) > 22, paste0(substr(x, 1, 21), "…"), x)
}

# Readable names for the Biocrates column codes the panel labels: "Tg.18.3_32.1." ->
# "TG 18:3_32:1", "Dg.18.3_18.3." -> "DG 18:3_18:3", "Hipacid" -> "Hippuric acid".
pretty_feature <- function(x) {
  x <- sub("^(Tg|Dg)[.](\\d+)[.](\\d+)_(\\d+)[.](\\d+)[.]$", "\\U\\1\\E \\2:\\3_\\4:\\5", x, perl = TRUE)
  dplyr::recode(x, "Hipacid" = "Hippuric acid")
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
plot_imic_volcano_panel <- function(res, title = "", n_top_vars = 5, overlap_n = 20,
                                    title_size = 8) {
  # BH critical raw p (0_figure-functions.R): the largest raw p whose q is <= 0.05
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
      pt_size    = ifelse(sig_status == "Significant after FDR", 1.4, 1),
      # tertiary panels label with the biomarker code (e.g. "Tg.18.3_30.0.", printed as
      # "TG 18:3_30:0" by pretty_feature()); label_f collapses every TG to "Triacylglyceride".
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
    # FDR-significant points: category colour and category symbol (colour-independent
    # cue), drawn slightly larger; the other tiers keep circle (n.s.) / square
    # (nominal only) in black. fill carries the colour for the filled shape 25.
    geom_point(aes(colour = color_var, fill = color_var, shape = color_var, size = pt_size),
               alpha = 0.75) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    # P < 0.05 (grey) and Q < 0.05 (green) solid, as in Fig 3
    geom_hline(yintercept = p_line, colour = "#BAB0AC", linewidth = 0.4) +
    geom_hline(yintercept = q_line, colour = "#59A14F", linewidth = 0.4) +
    xlab("") + ylab("") + ggtitle(title) +
    scale_color_manual(values = final_color_palette, na.value = "#999999") +
    scale_fill_manual(values = final_color_palette, na.value = "#999999", guide = "none") +
    scale_shape_manual(values = final_shape_palette, guide = "none") +
    scale_size_identity() +
    guides(color = guide_legend(title = NULL)) +
    # x fixed to the shared tertiary effect range (breaks at -1/0/1); y gets
    # ~12% headroom so top-feature labels are not jammed against the ceiling.
    scale_x_continuous(breaks = c(-1, 0, 1), limits = VOLCANO_XLIM) +
    scale_y_continuous(breaks = seq(0, ceiling(ymax_data), by = 2),
                       limits = c(0, ymax_data * 1.12)) +
    theme_imic(base_size = 8) +   # shared theme (Helvetica, font-size floors)
    theme(legend.position = "none",
          # keep a light border so the 8 volcano facets stay delineated (theme_imic drops it)
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3),
          # left-justified so the 3-line facet title sits top-LEFT, clear of the
          # top-right inset bar chart (they collided when centered + bold).
          plot.title = element_text(size = title_size, face = "plain", hjust = 0, lineheight = 0.9,
                                    margin = margin(b = 1)),
          plot.margin = margin(t = 1, r = 3, b = -4, l = 0),
          panel.spacing = unit(0, "pt"))

  top_vars <- tt_volcano %>% filter(pval_adj < q_cut) %>%
    arrange(-logPval) %>% head(n = n_top_vars) %>%
    mutate(lab = short_label(pretty_feature(lab_src)))
  # No ylim cap on the labels: a cap pins the label of any point above it to the cap,
  # where ggrepel cannot separate it from its neighbours. The inset sits above the
  # panel border, and ggrepel keeps labels inside the panel.
  p <- p + geom_text_repel(data = top_vars, aes(label = lab),
                           max.overlaps = getOption("ggrepel.max.overlaps", default = overlap_n),
                           size = 2.5, alpha = 0.5,
                           seed = 123,   # draw-time placement; required for byte-reproducible Fig 5A
                           # max.time = Inf stops on the iteration count alone, so placement
                           # cannot depend on machine speed (as in Fig 3B)
                           box.padding = 0.3, max.time = Inf, max.iter = 1e5)

  # Inset in the top-right corner of the panel's canvas. y = 1 in this [0, 1] canvas
  # is the top of the whole grob including the title band, not the panel border.
  ggdraw() + draw_plot(p, 0, 0, 1, 1) +
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
# category gets an imic_cat_cols colour (study_colors.R; rainbow fallback beyond 7).
build_palette <- function(categories) {
  categories <- unique(categories[!is.na(categories)])
  # Triglycerides first so it gets imic_cat_cols[1] = blue, as in the submitted
  # figure; the rest are sorted so the mapping is stable across re-runs.
  ordered <- c(intersect("Triglycerides", categories),
               sort(setdiff(categories, "Triglycerides")))
  fixed <- c("Not Significant" = "#000000", "Significant before FDR" = "#000000")
  # colourblind-safe set, paired 1:1 with imic_cat_shapes (build_shape_palette())
  cat_colors <- imic_cat_cols[seq_len(min(length(ordered), length(imic_cat_cols)))]
  names(cat_colors) <- ordered[seq_along(cat_colors)]
  if (length(ordered) > length(cat_colors)) {
    extra <- grDevices::rainbow(length(ordered) - length(cat_colors))
    names(extra) <- ordered[(length(cat_colors) + 1):length(ordered)]
    cat_colors <- c(cat_colors, extra)
  }
  c(fixed, cat_colors)
}

# Symbol for each category, in the same order build_palette() hands out colours, so
# colour i and shape i always travel together (non-significant tiers: circle / square).
build_shape_palette <- function(palette) {
  cats <- setdiff(names(palette), c("Not Significant", "Significant before FDR"))
  shp  <- rep_len(imic_cat_shapes, length(cats))
  c("Not Significant" = 16, "Significant before FDR" = 15, stats::setNames(shp, cats))
}

# Standalone legend (categories only) shown in an empty grid cell: colour and symbol,
# listed top-to-bottom in the inset bars' left-to-right (alphabetical) order.
create_category_legend <- function() {
  legend_colors <- final_color_palette[
    !names(final_color_palette) %in% c("Not Significant", "Significant before FDR")]
  cats <- sort(names(legend_colors))
  # fixed 1-unit key pitch in a fixed 10-unit band, title just above the first key
  # (same geometry as Fig 3A's legend)
  legend_data <- data.frame(category = cats, y = -(seq_along(cats) - 1), x = 1)
  ggplot(legend_data, aes(x = x, y = y, colour = category, fill = category, shape = category)) +
    geom_point(size = 3) +
    scale_colour_manual(values = legend_colors) +
    scale_fill_manual(values = legend_colors) +
    scale_shape_manual(values = final_shape_palette) +
    geom_text(aes(label = category), colour = "black", hjust = 0, nudge_x = 0.2, size = 2.5) +
    xlim(0.8, 3) + ylim(-9, 0.35) +
    theme_void() + theme(legend.position = "none") +
    ggtitle("Significant Categories") +
    theme(plot.title = element_text(size = 8, hjust = 0.5, margin = margin(b = 0)),
          plot.margin = margin(t = 6, b = 0))
}

# ===========================================================================
# PANEL B -- tertiary MSEA signed-enrichment-ratio volcano. Input columns:
#   study, timepoint, contrast, direction, pathway, total, expected, hits,
#   raw_p, fdr_native, enrichment_ratio
# where `enrichment_ratio` is already signed by direction (negative = down).
# ===========================================================================
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
  # Drawn by the shared renderer, so Fig 5B and Fig 6A follow the same conventions.
  # Submitted 5B style: colour all nominally significant pathways by study (no
  # "Sig before FDR" open tier), P<0.05 grey + Q<0.05 green lines, black vline at 0.
  # min_size = 1: the submitted 5B applied no set-size floor, so 1-2 member lipid
  # pathways (Ketone Body, Phytanic Ox, ...) are shown. Sized as a half-page
  # sub-panel (105 x 99 mm) of the full-page Fig 5 composite.
  # Labels are restricted to the pathways discussed in the Results (every point is
  # still plotted; only the text labels are limited). Names are exact data-column
  # names (note "Mitochondrial Electron Transport Chain" for the 3-4-month
  # up-reversal).
  label_allow_5b <- c(
    "Methionine Metabolism", "Glutamate Metabolism", "Glycine and Serine Metabolism",
    "Arginine and Proline Metabolism", "Fatty Acid Biosynthesis", "Bile Acid Biosynthesis",
    "Ammonia Recycling", "Glutathione Metabolism", "Carnitine Synthesis",
    "Mitochondrial Electron Transport Chain", "Ketone Body Metabolism",
    "Phytanic Acid Peroxisomal Oxidation")
  return(render_msea_panelB(msea_csv, out_png, title = NULL, min_size = 1,
                            submitted_style = TRUE, upper_line = "fdr",
                            label_which = "nominal", label_allow = label_allow_5b, vline0 = TRUE,
                            abbr_fun = abbr_5b, xlab = "Enrichment Ratio",
                            width_in = 210/25.4/2, height_in = 297/25.4/3,
                            show_legend = show_legend))
}

# One shared "Study" legend for Panels B+C (both coloured by the same canonical
# imic_study_cols map, plus "Not Significant"), used in place of each panel's
# own per-panel legend so the composite carries just one.
create_study_legend <- function() {
  study_cols_all <- c(imic_study_cols, "Not Significant" = "grey75")
  legend_df <- data.frame(x = 1, y = seq_along(study_cols_all),
                          grp = factor(names(study_cols_all), levels = names(study_cols_all)))
  p <- ggplot(legend_df, aes(x, y, color = grp, shape = grp)) +
    geom_point(size = 2) +
    scale_color_manual(values = study_cols_all, name = "Study") +
    scale_shape_manual(values = imic_shapes_for(names(study_cols_all)), name = "Study") +  # study symbols (CVD cue)
    guides(colour = guide_legend(nrow = 1), shape = guide_legend(nrow = 1)) +
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
  # x title in its own strip below the grid so it cannot sit on the bottom row's ticks
  grid_labeled <- plot_grid(
    ggdraw(volcano_grid) +
      draw_label(imic_logp_title, x = 0.0025, y = 0.5, angle = 90,
                 hjust = 0.5, vjust = 1, size = 8),
    ggdraw() + draw_label("Scaled Average Treatment Effect", x = 0.5, y = 0.5, size = 8),
    ncol = 1, rel_heights = c(1, 0.028))

  bottom_row <- plot_grid(
    png_to_ggdraw(panel_b_png), png_to_ggdraw(panel_c_png),
    labels = c("B", "C"), ncol = 2, nrow = 1, rel_widths = c(1, 1))

  bottom_with_legend <- plot_grid(bottom_row, create_study_legend(),
                                  ncol = 1, rel_heights = c(1, 0.06))

  plot_grid(grid_labeled, bottom_with_legend,
            labels = c("A", ""), ncol = 1, nrow = 2,
            rel_heights = c(rel_a, 1))
}

# ===========================================================================
# DATA (tertiary targeted metabolites and lipids)
# ===========================================================================
combined_arms   <- readRDS("results/combined_intervention_effects_results_combined_arms.RDS") %>%
  filter(outcome_group == "tertiary", measure == "ATE")
stratified_arms <- readRDS("results/combined_intervention_effects_results_stratified_arms.RDS") %>%
  filter(outcome_group == "tertiary", measure == "ATE")

# Plotted: every study in its combined-arm framing (one contrast per study x visit).
res_combined <- combined_arms
# Arm-stratified results (Misame/Vital per arm; Elicit uses its combined-arm rows).
# Not plotted; they enter only the shared category palette and x-range below, which
# span both framings.
res_stratified <- bind_rows(stratified_arms %>% filter(study != "Elicit"),
                            combined_arms  %>% filter(study == "Elicit"))

res_combined$category   <- collapse_small_categories(res_combined$category)
res_stratified$category <- collapse_small_categories(res_stratified$category)

final_color_palette <- build_palette(c(as.character(res_combined$category),
                                       as.character(res_stratified$category)))
final_shape_palette <- build_shape_palette(final_color_palette)   # symbol per category (CVD cue)

# One shared x-scale for every volcano panel, from the tertiary effect range.
.est_all <- c(res_combined$est, res_stratified$est)
VOLCANO_XLIM <- c(floor(min(.est_all, na.rm = TRUE) * 10) / 10 - 0.1,
                  ceiling(max(.est_all, na.rm = TRUE) * 10) / 10 + 0.15)

category_legend <- create_category_legend()
blank_plot <- ggplot() + theme_void()

panel <- function(res, study, visit, contrast, title, title_size = 8)
  plot_imic_volcano_panel(res %>% filter(study == !!study, visit == !!visit,
                                         contrast == !!contrast), title = title,
                          title_size = title_size)

# ===========================================================================
# PANEL B (standalone PNG, embedded in the composite)
# ===========================================================================
# tertiary_msea_dual.csv comes from the local MetaboAnalystR runner
# (src/metaboanalyst/run-tertiary-msea-dual.R): a metabolite pass plus a
# LIPID MAPS-converted lipid pass of ORA against the reference metabolome, for
# combined and stratified arms, computed from the ATEs without the web tool.
build_panelB("results/metaboanalyst/tertiary_msea/tertiary_msea_dual.csv",
             "figures/figure5_panelB_msea.png", show_legend = FALSE)

# Panel C: prefer the no-legend version of the triglyceride fatty-acid composition
# panel (compose_figure() draws one shared Study legend for B+C instead); fall back
# to the legend-bearing version.
panel_c_png <- if (file.exists("figures/figure5_panelC_tg_composition_nolegend.png")) {
  "figures/figure5_panelC_tg_composition_nolegend.png"
} else if (file.exists("figures/figure5_panelC_tg_composition.png")) {
  "figures/figure5_panelC_tg_composition.png"
} else {
  stop("Panel C PNG not found: run figure-scripts/manuscript_figures/fig5C-triglyceride.R first")
}

# ===========================================================================
# Fig 5 -- combined arms (7 volcano panels + legend, 3x3)
# ===========================================================================
# Facet titles: study + visit on one visible row, no arm label (all combined-arm here).
# Two leading blank rows ("\n\n<title>") put the facet title on the bottom row of the
# 3-line title band, so the white space is above the title (under the top-of-panel
# histogram) rather than below it; the band height is unchanged.
c_e1 <- panel(res_combined, "Elicit", "1 mo.",      "Nico", "\n\nELICIT (1 mo.)")
c_m1 <- panel(res_combined, "Misame", "14-21 days", "BEP",  "\n\nMISAME-III (14-21 days)")
c_v1 <- panel(res_combined, "Vital",  "1.5 mo.",    "BEP",  "\n\nMumta-LW (1.5 mo.)")
c_m2 <- panel(res_combined, "Misame", "1-2 mo.",    "BEP",  "\n\nMISAME-III (1-2 mo.)")
c_v2 <- panel(res_combined, "Vital",  "2 mo.",      "BEP",  "\n\nMumta-LW (2 mo.)")
c_e2 <- panel(res_combined, "Elicit", "5 mo.",      "Nico", "\n\nELICIT (5 mo.)")
c_m3 <- panel(res_combined, "Misame", "3-4 mo.",    "BEP",  "\n\nMISAME-III (3-4 mo.)")

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
# 3-way export (PDF + EPS + PNG) via save_figure_3way() (0_figure-functions.R). Panel A
# is ggplot vector; panels B/C are embedded raster (draw_image), which cairo keeps as
# raster while leaving A vector.
save_figure_3way(fig5_main, "figure5", width = PAGE_W, height = PAGE_H)
cat("wrote figures/figure5.{png,pdf,eps}\n")
