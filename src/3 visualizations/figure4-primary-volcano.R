# figure4-primary-volcano.R
# =============================================================================
# Primary untargeted-metabolome volcano composite (manuscript Fig 3).
#
# Adapted from trenton scripts/old/figure_scripts/figure4.Rmd into a standalone,
# push-button R script, with the 2026-07-30 Andrew x Trenton figure change:
#   * MAIN figure now uses the COMBINED arms (one contrast per study x visit),
#     for consistency with the rest of the main text -> figures/figure4.png
#   * the ARM-STRATIFIED layout (Misame BEP/BEP, BEP/IFA, IFA/BEP; Mumta two
#     arms) is kept for the online supplement
#     -> figures/figure4_stratified_supplement.png
#
# Structure of the composite (unchanged from Trenton's figure): panel A is the
# grid of per-study/visit volcano panels built here; panels B and C are Trenton's
# pre-rendered MSEA network PNGs (figures/final-figure-4-b.png / -c.png), embedded
# as static images. B/C are the same in both versions (study-level networks); if
# the combined-arm main needs combined-arm B/C, those are regenerated separately.
#
# Run from the repo root:
#   Rscript "src/3 visualizations/figure4-primary-volcano.R"
# =============================================================================

suppressMessages({
  library(dplyr); library(tidyr); library(ggplot2)
  library(ggrepel); library(cowplot); library(magick)
})

# A4 page minus 0.5in margins (same as Trenton's figure geometry)
PAGE_W <- 8.27 - 2 * 0.5   # 7.27 in
PAGE_H <- 11.69 - 2 * 0.5  # 10.69 in

tableau10 <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728",
               "#9467BD", "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF")

# final_color_palette is a category -> colour map referenced by the panel/inset/
# legend builders below. It is assigned once (after the data is loaded and the
# category set is known) and read at call time, so defining the helpers first is
# fine.
final_color_palette <- NULL

# ---------------------------------------------------------------------------
# Helpers ported verbatim (behaviour-preserving) from figure4.Rmd
# ---------------------------------------------------------------------------

# Embed a PNG as a ggdraw canvas so it can sit in a cowplot grid.
png_to_ggdraw <- function(path) ggdraw() + draw_image(image_read(path))

# Shorten a metabolite name for the in-panel repel labels: drop the
# "(expressed as ...)" / "(niacin catabolite)" parenthetical qualifiers and cap
# the length, so long names don't spill past the narrow panel edges.
short_label <- function(x) {
  x <- sub("\\s*\\(.*$", "", x)          # remove from the first "(" onward
  x <- trimws(x)
  ifelse(nchar(x) > 22, paste0(substr(x, 1, 21), "…"), x)
}

# BH critical raw-p: the largest raw p whose BH q is still <= alpha. Used to draw
# the FDR reference line and to flag "significant after FDR".
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
# FDR-significant features (solid), drawn in the top-right of each volcano panel.
make_inset_with_labels <- function(df_sum, df_long2) {
  xmax_val <- max(df_sum$n_total) * 1.05  # pad x-range slightly
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
plot_imic_volcano_panel <- function(res, title = "", label_type = "label",
                                    n_top_vars = 5, overlap_n = 20) {
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
      label_f    = ifelse(sig_status == "Significant after FDR", label_f, ""))

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

  p <- ggplot(tt_volcano, aes(x = ATE, y = logPval)) +
    geom_point(aes(colour = color_var, shape = sig_status), size = 1, alpha = 0.75) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = p_line, linetype = "dashed", colour = tableau10[2]) +
    geom_hline(yintercept = q_line, linetype = "dotted", colour = tableau10[3]) +
    xlab("") + ylab("") + ggtitle(title) +
    scale_color_manual(values = final_color_palette, na.value = "#999999") +
    guides(color = guide_legend(title = NULL)) +
    # x tightened to the actual effect-size range (-0.63..1.54): drops the empty
    # left third and stops clipping the largest positive effect. y gets ~12%
    # headroom so top-feature labels are not jammed against the panel ceiling.
    scale_x_continuous(breaks = c(0, 1), limits = c(-0.7, 1.6)) +
    scale_y_continuous(breaks = seq(0, ceiling(max(tt_volcano$logPval)), by = 2),
                       limits = c(0, ceiling(max(tt_volcano$logPval)) * 1.12)) +
    theme_bw() +
    theme(legend.position = "none",
          axis.text = element_text(size = 7), axis.title = element_text(size = 7),
          plot.title = element_text(size = 6.8, lineheight = 0.9,
                                    margin = margin(b = 1)),
          plot.margin = margin(t = 1, r = 3, b = -4, l = 0),
          panel.spacing = unit(0, "pt"))

  top_vars <- tt_volcano %>% filter(pval_adj < q_cut) %>%
    arrange(-logPval) %>% head(n = n_top_vars) %>%
    mutate(lab = short_label(label_f))
  # label with the (shortened) full metabolite name, not the short biomarker code,
  # to match the published figure (e.g. "riboflavin" rather than "Ribo").
  p <- p + geom_text_repel(data = top_vars, aes(label = lab),
                           max.overlaps = getOption("ggrepel.max.overlaps", default = overlap_n),
                           size = 2.5, alpha = 0.5)

  # overlay the category inset in the top-right corner (the title occupies the
  # top-left). Kept compact so the up-regulated feature labels retain the
  # center-right space.
  ggdraw() + draw_plot(p, 0, 0, 1, 1) +
    draw_plot(inset, x = 0.72, y = 0.835, width = 0.265, height = 0.15)
}

# Categories rarer than `threshold` of features are pooled into "Other" so the
# palette and inset stay legible.
collapse_small_categories <- function(x, threshold = 0.04) {
  x <- as.character(x)
  tab <- prop.table(table(x))
  small <- names(tab[tab < threshold])
  factor(ifelse(x %in% small, "Other", x))
}

# Build the category -> colour map (Trenton's construction): the two "not/before
# FDR" states are black; each real category gets a tableau10 colour, with rainbow
# fallback if there are more categories than tableau colours.
build_palette <- function(categories) {
  categories <- sort(unique(categories[!is.na(categories)]))
  fixed <- c("Not Significant" = "#000000", "Significant before FDR" = "#000000")
  cat_colors <- tableau10[seq_len(min(length(categories), length(tableau10)))]
  names(cat_colors) <- categories[seq_along(cat_colors)]
  if (length(categories) > length(cat_colors)) {
    extra <- grDevices::rainbow(length(categories) - length(cat_colors))
    names(extra) <- categories[(length(cat_colors) + 1):length(categories)]
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

# Wrap a volcano grid with shared axis titles + the two pre-rendered MSEA panels.
compose_figure <- function(volcano_grid) {
  grid_labeled <- ggdraw(volcano_grid) +
    draw_label("Scaled Average Treatment Effect", x = 0.5, y = 0.0025,
               hjust = 0.5, vjust = 0, size = 7) +
    draw_label("-log10(P-value)", x = 0.0025, y = 0.5, angle = 90,
               hjust = 0.5, vjust = 1, size = 7)
  # Manuscript Fig 3 has TWO panels: A = primary volcano (above), B = the
  # reproducible B-vitamin MSEA (figure4-panelB-msea.R). There is NO panel C here
  # -- the triglyceride -> fatty-acid panel belongs to Fig 5 (tertiary outcomes).
  pick <- function(real, placeholder) if (file.exists(real)) real else placeholder
  panel_b <- png_to_ggdraw(pick("figures/figure4_panelB_msea.png", "figures/final-figure-4-b.png"))
  plot_grid(grid_labeled, panel_b,
            labels = c("A", "B"), align = "h", axis = "l",
            ncol = 1, nrow = 2, rel_heights = c(2 / 3, 1 / 3))
}

save_figure <- function(fig, path) {
  ggsave(filename = path, plot = fig, width = PAGE_W, height = PAGE_H,
         units = "in", dpi = 300, device = ragg::agg_png)
  cat("wrote", path, "\n")
}

# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------
# NOTE: this is the PRIMARY (targeted HM-nutrient / B-vitamin) volcano = manuscript
# Fig 3. Trenton's old figure4.Rmd filtered outcome_group=="tertiary", but on the
# current data that is the untargeted lipidome (Triglycerides/Ceramides/...), which
# does NOT match the published B-vitamin figure. The nutrient data (Micronutrient,
# B1/B2/B3/B6, Macronutrient) lives in outcome_group=="primary" — corrected here.
combined_arms   <- readRDS("results/combined_intervention_effects_results_combined_arms.RDS") %>%
  filter(outcome_group == "primary", measure == "ATE")
stratified_arms <- readRDS("results/combined_intervention_effects_results_stratified_arms.RDS") %>%
  filter(outcome_group == "primary", measure == "ATE")

# MAIN: every study in its combined-arm framing (one contrast per study x visit).
res_combined <- combined_arms
# SUPPLEMENT: Misame/Mumta stratified per arm; Elicit has no strata so it uses its
# combined-arm rows (matches Trenton's original figure).
res_stratified <- bind_rows(stratified_arms %>% filter(study != "Elicit"),
                            combined_arms  %>% filter(study == "Elicit"))

res_combined$category   <- collapse_small_categories(res_combined$category)
res_stratified$category <- collapse_small_categories(res_stratified$category)

# One shared palette across both figures so a category is the same colour in the
# main and supplement versions.
final_color_palette <- build_palette(c(as.character(res_combined$category),
                                       as.character(res_stratified$category)))
category_legend <- create_category_legend()
blank_plot <- ggplot() + theme_void()

panel <- function(res, study, visit, contrast, title)
  plot_imic_volcano_panel(res %>% filter(study == !!study, visit == !!visit,
                                         contrast == !!contrast), title = title)

# ---------------------------------------------------------------------------
# MAIN figure -- combined arms (7 panels, ordered by collection time)
# ---------------------------------------------------------------------------
c_e1 <- panel(res_combined, "Elicit", "1 mo.",      "Nico", "ELICIT\n1m\nNico.")
c_m1 <- panel(res_combined, "Misame", "14-21 days", "BEP",  "MISAME-III\n14-21d\nBEP")
c_v1 <- panel(res_combined, "Vital",  "1.5 mo.",    "BEP",  "Mumta-LW\n1.5m\nBEP")
c_m2 <- panel(res_combined, "Misame", "1-2 mo.",    "BEP",  "MISAME-III\n1-2m\nBEP")
c_v2 <- panel(res_combined, "Vital",  "2 mo.",      "BEP",  "Mumta-LW\n2m\nBEP")
c_e2 <- panel(res_combined, "Elicit", "5 mo.",      "Nico", "ELICIT\n5m\nNico.")
c_m3 <- panel(res_combined, "Misame", "3-4 mo.",    "BEP",  "MISAME-III\n3-4m\nBEP")

combined_grid <- plot_grid(
  c_e1, c_m1, c_v1, c_m2,
  c_v2, c_e2, c_m3, category_legend,
  ncol = 4, nrow = 2, align = "hv", axis = "lrtb")

# Written to a distinct name (not figure4.png) so the published figure is not
# overwritten until the combined version is reviewed and the MSEA B/C panels are
# regenerated for combined arms.
save_figure(compose_figure(combined_grid), "figures/figure4_combined_arms.png")

# ---------------------------------------------------------------------------
# SUPPLEMENT figure -- arm-stratified (Trenton's original 6x3 layout, 15 panels)
# ---------------------------------------------------------------------------
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

save_figure(compose_figure(stratified_grid), "figures/figure4_stratified_supplement.png")
