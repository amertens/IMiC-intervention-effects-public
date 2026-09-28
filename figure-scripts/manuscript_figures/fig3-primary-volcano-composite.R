# =============================================================================
# fig3-primary-volcano-composite.R
#
# Builds Fig 3. Panel A: volcano plots of the combined-arm intervention effects on
# the primary milk outcomes, one panel per study and visit; FDR-significant points
# are coloured and shaped by outcome category and the top 3 are labelled, and an
# inset bar chart shows, per category, the share of outcomes tested (light) and
# FDR-significant (dark). Panel B is the KEGG pathway-impact plot written by
# fig3B-pathway.R, embedded as a PNG.
#
# Inputs:  results/combined_intervention_effects_results_combined_arms.RDS
#            (written by src/2 analysis/clean_results.R)
#          figures/figure3_panelB_msea.png (fig3B-pathway.R; run it first)
# Outputs: figures/figure3.{png,pdf,eps}
# [needs an on-request file] the combined-results RDS is a feature-level result
# file that is not shipped (size).
# =============================================================================
suppressPackageStartupMessages(library(here))
# the relative '../../' paths below resolve from this folder
setwd(file.path(here::here(), "figure-scripts/manuscript_figures"))

library(tidyverse)
library(cowplot)
library(magick)     # image_read() for the embedded Panel B PNG
library(ggrepel)

# theme_imic(), get_bh_cutoff(), label helpers and save_figure_3way()
source(file.path(here::here(), "figure-scripts/0_figure-functions.R"))
source(file.path(here::here(), "figure-scripts/manuscript_figures/study_colors.R"))  # CVD-safe palettes/shapes

# ---- Page dimensions -----------------------------------------------------------
# Portrait, W/H = 0.682: a full-width volcano facet grid (~2/3 height) over a
# left-justified, roughly square pathway Panel B (~1/3 height).
page_w <- 7.25                                # 2-column width
page_h <- round(page_w / 0.682, 2)            # 10.63 in

blank_plot <- ggplot() + theme_void()

# wrap a PNG as a ggdraw() canvas
png_to_ggdraw <- function(path){
  ggdraw() + draw_image(image_read(path))
}

# Inset bar chart for one volcano panel: per category, all outcomes tested (light)
# overlaid with the FDR-significant ones (dark).
  make_inset_with_labels <- function(df_sum, df_long2) {

    xmax_val <- max(df_sum$n_total) * 1.05  # Pad x-range slightly for labels

   p = ggplot(df_long2) +
      # Draw bars (no spacing)
      geom_rect(
        aes(xmin = 0, xmax = count,
            ymin = ymin, ymax = ymax,
            fill = category, alpha = type),
        color = "black", size = 0.15
      ) +
      scale_fill_manual(values = final_color_palette) +
      scale_alpha_manual(values = c(n_total = 0.15, n_sig = 1)) +
      scale_y_continuous(expand = c(0,0)) +
      scale_x_continuous(expand = c(0,0), limits = c(0, xmax_val)) +
      coord_flip(clip = "off") +

      theme_void(base_size = 6) +
      theme(
        panel.background = element_rect(fill = "white", color = NA),
        legend.position = "none",
        plot.margin = margin(0,0,0,0)
      )

   # add a thin border
  inset_with_border <- p + annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = NA, color = "black", size = 0.1)
  inset_with_border

  }


# One volcano panel with its category inset: category colour plus symbol for
# FDR-significant points, P and Q reference lines, and the top n_top_vars labelled
# with short names.
plot_imic_volcano_panel <- function(res,
                                    title              = "",
                                    label_type         = "label",
                                    n_top_vars         = 3,
                                    overlap_n          = 8,
                                    title_size         = 8) {




  # Q < 0.05 expressed as the largest raw P that is FDR-significant
  q_cut= get_bh_cutoff(df=res,
                p_col  = "pval",
                q_col  = "pval_adj",
                alpha  = 0.05,
                return = c("raw"))

  tt_volcano <- res %>%
    arrange(pval_adj) %>%
    mutate(
      ATE        = est,
      logPval    = -log10(pval),
      sig_status = case_when(
        pval <= q_cut        ~ "Significant after FDR",
        pval     < 0.05        ~ "Significant before FDR",
        TRUE                         ~ "Not Significant"
      ),
      color_var  = ifelse(sig_status == "Significant after FDR",
                          category, sig_status),
      pt_size    = ifelse(sig_status == "Significant after FDR", 1.4, 1),
      label_f    = ifelse(sig_status == "Significant after FDR", label_f, "")
    )


  category_summary <- tt_volcano %>%
    group_by(category) %>%
    summarise(
      n_total = n(),
      n_sig = sum(sig_status == "Significant after FDR"),
      rate_sig = n_sig / n_total
    ) %>%
    mutate(
      share_sig = n_sig / sum(n_sig)
    )


    df_sum <- category_summary %>%
    mutate(
      pct_within = n_sig / n_total,
      pct_total  = n_total / sum(n_total),
      cat_num    = as.numeric(factor(category)),
      ymin       = cat_num - 0.5,
      ymax       = cat_num + 0.5,
      non_sig    = n_total - n_sig
    )

  df_long2 <- df_sum %>%
    select(category, cat_num, ymin, ymax, n_sig, n_total) %>%
    pivot_longer(
      cols = c(n_sig, n_total),
      names_to = "type",
      values_to = "count"
    )


  inset <- make_inset_with_labels(df_sum, df_long2)

  # reference lines
  p_line <- -log10(0.05)
  q_line <- -log10(q_cut)

  p <- ggplot(tt_volcano, aes(x = ATE, y = logPval)) +
    # FDR-significant points: category colour and category symbol (colour-independent
    # cue), drawn slightly larger; the other tiers keep circle (n.s.) / square
    # (nominal only) in black. fill carries the colour for the filled shape 25.
    geom_point(aes(colour = color_var, fill = color_var, shape = color_var, size = pt_size),
               alpha = 0.75) +
    geom_vline(xintercept = 0,            linetype = "dashed") +
    # P < 0.05 (grey) and Q < 0.05 (green) lines solid, matching Panel B
    geom_hline(yintercept = p_line,       colour = "#BAB0AC", linewidth = 0.4) +
    geom_hline(yintercept = q_line,       colour = "#59A14F", linewidth = 0.4) +
    xlab("") + ylab("") + ggtitle(title) +
    scale_color_manual(values = final_color_palette, na.value = "#999999") +
    scale_fill_manual(values = final_color_palette, na.value = "#999999", guide = "none") +
    scale_shape_manual(values = final_shape_palette, guide = "none") +
    scale_size_identity() +
    guides(color = guide_legend(title = NULL)) +

    scale_x_continuous(breaks = c(-1,0,1),
                       limits = c(-1,1.5)) +
    scale_y_continuous(
      breaks = scales::breaks_pretty(n = 4),
      # ~8% headroom so the top point's label sits just below the panel border,
      # clear of the category histogram that caps the top-right corner.
      limits = c(0, ceiling(max(tt_volcano$logPval)) * 1.08)
    ) +

    theme_imic(base_size = 8) +   # shared theme (Helvetica, font-size floors)
    theme(
      legend.position = "none",
      # keep a light border so the 8 volcano facets stay delineated (theme_imic drops it)
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3),
      # The facet title sits on the bottom row of the title band (call sites pass
      # "\n\n<title>"), so the white space is above the title, just under the top-of-panel
      # histogram. The top plot.margin adds that space; the title's small bottom margin
      # keeps it snug above the panel border.
      plot.title = element_text(size = title_size, face = "plain", hjust = 0, lineheight = 0.9, margin = margin(t = 0, b = 1)),
      plot.margin = margin(t = 10, r = 0, b = -4, l = 0),
      panel.spacing = unit(0,"pt")
    )

  # label top variables after FDR
  top_vars <- tt_volcano %>%
    filter(pval_adj < q_cut) %>%
    arrange(-logPval) %>%
    head(n = n_top_vars)

  p <- p + geom_text_repel(
    data            = top_vars,
    aes(label       = plotmath_label(abbr_label(biomarker))),   # short forms so labels fit the narrow panels (Fig 3A)
    parse           = TRUE,   # B-vitamin subscripts via plotmath (Arial has no subscript-digit glyphs for PDF/EPS)
    max.overlaps    = getOption("ggrepel.max.overlaps", default = overlap_n),
    size            = 2.5,
    alpha           = 0.5,
    seed            = 123)   # draw-time label placement; without it Fig 3 is not byte-reproducible

    # combine volcano + inset. The inset's BOTTOM (y) is placed flush on the plot's
    # top border so the histogram caps the panel with no white gap; with the taller
    # top plot.margin above, the panel top sits lower, so y is lowered to match.
  ggdraw() +
    draw_plot(p, 0, 0, 1, 1) +
    draw_plot(inset,
              x = 0.675,   # adjust horizontal location
              y = 0.782,   # histogram baseline flush on the plot's top border
              width = 0.321,
              height = 0.1725)

}


# Combined-arm results: one contrast per study x visit.
res_comb_all <- readRDS(file="../../results/combined_intervention_effects_results_combined_arms.RDS") %>%
  filter(outcome_group=="primary", measure=="ATE")
res <- res_comb_all


# Category -> colour map, keyed by name so a category keeps its colour across
# re-runs (assigning colours in first-appearance order reshuffled them whenever
# the category order changed). Colours are the colourblind-safe imic_cat_cols
# (study_colors.R), and every category also gets its own symbol (category_shapes)
# so it is identifiable without colour. Non-significant and nominal-only points
# are black.
category_color_palette <- c(
  "Not Significant"        = "#000000",
  "Significant before FDR" = "#000000"
)
category_colors <- c(
  "Other B vitamins" = imic_cat_cols[1],   # blue
  "B2"               = imic_cat_cols[2],   # orange
  "B3 or related"    = imic_cat_cols[3],   # bluish green
  "B1"               = imic_cat_cols[4],   # vermillion
  "B6"               = imic_cat_cols[5],   # reddish purple
  "Micronutrient"    = imic_cat_cols[6],   # sky blue
  "Macronutrient"    = imic_cat_cols[7]    # wine
)
category_shapes <- setNames(imic_cat_shapes, names(category_colors))
# defensive fallback for any category not in the fixed map
.unmapped <- setdiff(unique(res$category[!is.na(res$category)]), names(category_colors))
if (length(.unmapped)) category_colors <- c(category_colors, setNames(rainbow(length(.unmapped)), .unmapped))
final_color_palette <- c(category_color_palette, category_colors)
if (length(.unmapped)) category_shapes <- c(category_shapes, setNames(rep(1, length(.unmapped)), .unmapped))
final_shape_palette <- c("Not Significant" = 16, "Significant before FDR" = 15, category_shapes)


# ---- Panel A: one volcano per study x visit ------------------------------------
# Two leading blank rows ("\n\n<title>") place the facet title on the bottom row of the
# 3-line title band, so the white space is above the title (under the top-of-panel
# histogram) rather than below it. The 3-line band height is unchanged, so the plot area
# and the inset histogram stay put.
c_m1 <- plot_imic_volcano_panel(res_comb_all %>% filter(study=="Misame", visit=="14-21 days", contrast=="BEP"),  title="\n\nMISAME-III (14-21 days)")
c_m2 <- plot_imic_volcano_panel(res_comb_all %>% filter(study=="Misame", visit=="1-2 mo.",    contrast=="BEP"),  title="\n\nMISAME-III (1-2 mo.)")
c_m3 <- plot_imic_volcano_panel(res_comb_all %>% filter(study=="Misame", visit=="3-4 mo.",    contrast=="BEP"),  title="\n\nMISAME-III (3-4 mo.)")
c_v1 <- plot_imic_volcano_panel(res_comb_all %>% filter(study=="Vital",  visit=="1.5 mo.",    contrast=="BEP"),  title="\n\nMumta-LW (1.5 mo.)")
c_v2 <- plot_imic_volcano_panel(res_comb_all %>% filter(study=="Vital",  visit=="2 mo.",      contrast=="BEP"),  title="\n\nMumta-LW (2 mo.)")
c_e1 <- plot_imic_volcano_panel(res_comb_all %>% filter(study=="Elicit", visit=="1 mo.",      contrast=="Nico"), title="\n\nELICIT (1 mo.)")
c_e2 <- plot_imic_volcano_panel(res_comb_all %>% filter(study=="Elicit", visit=="5 mo.",      contrast=="Nico"), title="\n\nELICIT (5 mo.)")


# Colour + symbol key for the significant categories (fills one empty grid cell).
create_category_legend <- function() {
  # Get only the category colors (excluding "Not Significant" and "Significant before FDR")
  legend_colors <- final_color_palette[!names(final_color_palette) %in% c("Not Significant", "Significant before FDR")]
  cats <- sort(names(legend_colors))   # = factor(category) order = inset bar order (left to right)

  # Top-to-bottom in inset order on a fixed 1-unit pitch; the y range below is a fixed
  # 10-unit band, so the title sits just above the first key
  legend_data <- data.frame(
    category = cats,
    y = -(seq_along(cats) - 1),
    x = 1
  )

  # Create the legend plot: each category's colour and symbol
  legend_plot <- ggplot(legend_data, aes(x = x, y = y, colour = category, fill = category, shape = category)) +
    geom_point(size = 3) +
    scale_colour_manual(values = legend_colors) +
    scale_fill_manual(values = legend_colors) +
    scale_shape_manual(values = final_shape_palette) +
    geom_text(aes(label = category), colour = "black", hjust = 0, nudge_x = 0.2, size = 2.5) +
    xlim(0.5, 3) +
    ylim(-9, 0.35) +
    theme_void() +
    theme(legend.position = "none") +
    ggtitle("Significant Categories") +
    theme(plot.title = element_text(size = 8, hjust = 0.5, margin = margin(b = 0)),
          plot.margin = margin(t = 6, b = 0))

  return(legend_plot)
}

category_legend <- create_category_legend()


# ---- Panel B: KEGG pathway impact (PNG written by fig3B-pathway.R) --------------
p3b_comb_path  <- "../../figures/figure3_panelB_msea.png"
p3b_comb  <- if (file.exists(p3b_comb_path))  png_to_ggdraw(p3b_comb_path)  else blank_plot

# Shared axis titles for the volcano grid. 8 pt = the size Panel B's axis titles
# print at (fig3B-pathway.R EMBED_SCALE). The x title gets its own strip below the
# grid so it cannot sit on the bottom row's ticks.
label_grid <- function(g) plot_grid(
  ggdraw(g) +
    draw_label(imic_logp_title, x = 0.0025, y = 0.5, angle = 90, hjust = 0.5, vjust = 1, size = 8),
  ggdraw() + draw_label("Scaled Average Treatment Effect", x = 0.5, y = 0.5, size = 8),
  ncol = 1, rel_heights = c(1, 0.028))

# Panel A grid, study-major (MISAME-III / Mumta-LW / ELICIT rows), in the same order
# as Fig 5A. The legend fills Mumta-LW's empty 3rd cell; the blank cell falls to
# the 3rd row.
comb_grid <- plot_grid(
  c_m1, c_m2, c_m3,
  c_v1, c_v2, category_legend,
  c_e1, c_e2, blank_plot,
  ncol = 3, nrow = 3, align = "hv", axis = "lrtb")

# Panel B sits bottom-left as a roughly square plot with blank space to its right
# (the pathway plot does not span the full page width).
b_row_comb  <- plot_grid(p3b_comb,  blank_plot, ncol = 2, rel_widths = c(0.52, 0.48))

fig3 <- plot_grid(label_grid(comb_grid), b_row_comb, labels = "AUTO",
                  ncol = 1, nrow = 2, rel_heights = c(0.68, 0.32))

save_figure_3way(fig3, name = "figure3",
                 width = page_w, height = page_h,
                 dir = file.path(here::here(), "figures"))
