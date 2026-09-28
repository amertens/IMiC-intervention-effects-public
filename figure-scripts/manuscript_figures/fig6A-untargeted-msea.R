# =============================================================================
# fig6A-untargeted-msea.R
#
# Builds Fig 6A, the MSEA of the untargeted milk metabolome, with the shared renderer
# render_msea_panelB() in the style of Figs 3B and 5B: every pathway with P < 0.05
# filled in its study colour and symbol, the rest grey; reference lines at P < 0.05
# (grey) and the Q < 0.05 FDR frontier (green); only FDR-significant pathways are
# labelled (at most 8). Also writes the key shared by Panels A-C, which
# fig6-composite.R places once under Panel C.
#
# Inputs:  results/metaboanalyst/untargeted_msea/untargeted_msea_combined.csv
#            (src/metaboanalyst/run-untargeted-msea.R; its FDR-significant rows are Table S5)
# Outputs: figures/figure6_panelA_untargeted_msea.png, figures/figure6_legend_ABC.png
# Run from the repo root:  Rscript "figure-scripts/manuscript_figures/fig6A-untargeted-msea.R"
# =============================================================================
source("figure-scripts/manuscript_figures/render_msea_panelB.R")
source("figure-scripts/manuscript_figures/fig6_layout.R")   # FIG6_* geometry + print-scale text sizes

MSEA_CSV <- "results/metaboanalyst/untargeted_msea/untargeted_msea_combined.csv"
OUT_PNG  <- "figures/figure6_panelA_untargeted_msea.png"
# FIG6_Y_MAX (shared A-C -log10(P) height; 9 covers Panel C's max ~7.72 with repel
# headroom) and FIG6_PANEL_H_MM come from fig6_layout.R. Fig 6 A/B/C are drawn without
# their own legends at one shared size, and a single key (this panel's legend: 3 studies
# + Not Significant) is written separately and placed once under Panel C by
# fig6-composite.R.
LEGEND_PNG <- "figures/figure6_legend_ABC.png"

if (sys.nframe() == 0) {
  # submitted_style = TRUE: P < 0.05 filled in study colour, non-sig grey (as Figs 3B/5B);
  # upper_line = "fdr" draws the Q<0.05 green frontier; label_which = "fdr" labels only
  # FDR-significant pathways.
  p <- render_msea_panelB(MSEA_CSV, OUT_PNG, title = NULL, min_size = 1,
                     submitted_style = TRUE, upper_line = "fdr",
                     label_which = "fdr", label_max = 8, vline0 = TRUE,
                     xlab = "Fold Enrichment (signed by direction)",   # same title as Panels B-C
                     study_cols = imic_study_cols,  # shared study colours (ELICIT blue / MISAME orange / Mumta purple)
                     y_max = FIG6_Y_MAX,  # shared Fig 6 A/B/C y-axis height (matches submitted style)
                     x_margin_hi = 1.30,  # extra right-edge headroom so the "Glucose-Alanine Cycle" label
                     # (rightmost point) doesn't get clipped by the panel border
                     # base_size = 8 and label_size = 2.2 match Fig 3A/5A's text sizing
                     # (theme_imic(base_size = 8), repel size 2.5): the boxed labels read larger
                     # than Fig 3A/5A's plain, semi-transparent labels at the same nominal size,
                     # so the size is trimmed a little. Fig 5B keeps the renderer defaults.
                     base_size = 8, label_size = 2.2,
                     show_legend = FALSE,   # shared key drawn once under Panel C (see LEGEND_PNG)
                     axis_text_size  = fig6_pt(FIG6_TICK_PT,  FIG6_ABC_SCALE),   # print at Figs 3/5 sizes
                     axis_title_size = fig6_pt(FIG6_TITLE_PT, FIG6_ABC_SCALE),
                     width_in = FIG6_PANEL_W_MM/25.4, height_in = FIG6_PANEL_H_MM/25.4)

  # The shared Fig 6 A-C key = this panel's legend, drawn on its own.
  g <- ggplotGrob(p + theme(legend.position = "bottom",
                            legend.text  = element_text(size = fig6_pt(FIG6_LEGEND_PT, FIG6_ABC_SCALE)),
                            legend.title = element_text(size = fig6_pt(FIG6_TITLE_PT, FIG6_ABC_SCALE)),
                            legend.key.size = unit(0.35 / FIG6_ABC_SCALE, "cm")))
  leg <- g$grobs[[which(g$layout$name %in% c("guide-box-bottom", "guide-box"))[1]]]
  ragg::agg_png(LEGEND_PNG, width = FIG6_PANEL_W_MM, height = 25, units = "mm", res = 300,
                background = "white")
  grid::grid.newpage(); grid::grid.draw(leg); invisible(dev.off())
  cat("wrote", LEGEND_PNG, "\n")
}
