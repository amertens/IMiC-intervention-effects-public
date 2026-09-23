# fig6A-untargeted-msea.R
# =============================================================================
# Manuscript Figure 6, Panel A: MSEA of the untargeted milk metabolome.
#
# Thin wrapper around the shared renderer render_msea_panelB(). EXPLORATORY-FDR style
# (2026-08-14): 6A is second-platform untargeted corroboration, so the panel is drawn to
# read as exploratory, not as strong stand-alone evidence:
#   - points coloured by FDR (Q<0.05 -> study colour, solid); nominally-significant-but-
#     not-FDR shown as OPEN circles ("Sig before FDR"); non-significant grey,
#   - reference lines at P < 0.05 (grey) and the Q < 0.05 FDR frontier (green),
#   - only FDR-significant pathways are labelled (capped).
# This replaces the earlier submitted-mimic style (colour by nominal p, all filled).
#
# Input: results/metaboanalyst/untargeted_msea/untargeted_msea_combined.csv, produced by
# src/metaboanalyst/run-untargeted-msea.R. Background = MetaboAnalyst WHOLE-LIBRARY default
# (the un-inflated, submission-faithful universe; see that script's 2026-08-14 reconciliation
# header). The earlier 615-name reference over-stated significance (annotation/reference
# sensitivity, e.g. Bile Acid Biosynthesis n.s. in the submission but p=7e-5 under 615).
#
# Run from repo root:  Rscript "figure-scripts/manuscript_figures/fig6A-untargeted-msea.R"
# =============================================================================
source("figure-scripts/manuscript_figures/render_msea_panelB.R")

MSEA_CSV <- "results/metaboanalyst/untargeted_msea/untargeted_msea_combined.csv"
OUT_PNG  <- "figures/figure6_panelA_untargeted_msea.png"
# Shared y-axis ceiling across Fig 6 Panels A/B/C, matching the submitted figure's style
# (all three panels share one -log10(P) height). 9 comfortably covers the tallest point
# across all three panels' source data (Panel C's max ~7.72) with label-repel headroom;
# see figure-scripts/manuscript_figures/fig6B-mummichog.R and fig6C-proteomics.R for the
# matching constant (2026-08-26).
FIG6_Y_MAX <- 9

if (sys.nframe() == 0) {
  # EXPLORATORY-FDR style: submitted_style = FALSE colours by FDR (Q<0.05 study-colour
  # solid, nominal-only OPEN "Sig before FDR", non-sig grey); upper_line = "fdr" draws the
  # Q<0.05 green frontier; label_which = "fdr" labels only FDR-significant pathways. Under
  # the whole-library background there are far fewer FDR-sig pathways, so the panel reads
  # as modest/exploratory -- matching the submitted magnitude and the manuscript framing.
  render_msea_panelB(MSEA_CSV, OUT_PNG, title = NULL, min_size = 1,
                     submitted_style = FALSE, upper_line = "fdr",
                     label_which = "fdr", label_max = 8, vline0 = TRUE,
                     nominal_shape = 1,   # OPEN circles for the nominal-only tier: its light blue sat too close to
                     # ELICIT blue for colourblind readers when both were filled circles (2026-09-23)
                     study_cols = imic_study_cols,  # canonical study colours (ELICIT blue / MISAME orange / Mumta green)
                     y_max = FIG6_Y_MAX,  # shared Fig 6 A/B/C y-axis height (matches submitted style)
                     x_margin_hi = 1.30,  # extra right-edge headroom so the "Glucose-Alanine Cycle" label
                     # (rightmost point) doesn't get clipped by the panel border
                     # base_size = 8 (was 9) + label_size = 2.2 (was default 2.5): matches Fig 3A/5A's
                     # text sizing (theme_imic(base_size=8), repel size 2.5) more closely -- Fig 6's
                     # boxed geom_label_repel labels read visually larger than Fig 3A/5A's plain,
                     # semi-transparent geom_text_repel at the identical nominal size because of the
                     # label box itself, so the size is trimmed a bit further to compensate (reported
                     # 2026-08-26: "text seems large" vs the other figures). Fig 5B (this same shared
                     # renderer) is intentionally left at its own default -- this scoped to Fig 6 only.
                     base_size = 8, label_size = 2.2,
                     # half-page panel of the full-page Fig 6 (7.25 in), 105x99 mm (submitted panel size).
                     width_in = 210/25.4/2, height_in = 297/25.4/3)
}
