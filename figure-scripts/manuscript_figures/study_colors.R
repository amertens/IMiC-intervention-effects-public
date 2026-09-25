# study_colors.R -- SINGLE SOURCE OF TRUTH for IMiC study colours AND shapes.
#
# Every study-coloured manuscript panel (Fig 3B pathway, Fig 5B/6A MSEA, Fig 5C
# fatty-acid, Fig 6B mummichog, Fig 6C proteome GO) must use THESE named maps so a
# study keeps the same colour and symbol in every figure.
#
# COLOURBLIND-SAFE since 2026-09-23 (Science editor, round 2: "make sure [figures] are
# colorblind-accessible, and either avoid red-green color coding or at least make sure
# there is a backup feature ... such as different symbol shapes"). The previous map
# (ELICIT #1F77B5 / MISAME-III #FE7F01 / Mumta-LW #2CA02C, set 2026-08-13) put
# MISAME-III and Mumta-LW at CIEDE2000 dE 1.5 under simulated protanopia -- i.e. the same
# colour -- with no second cue. Now:
#   * colours = Okabe-Ito blue / orange / reddish purple. Mumta-LW was bluish green
#     until 2026-09-24; blue vs bluish green read as two dark cool colours in small points
#     (and merge under tritanopia, dE 12.1), so it moved to a separate hue family. Worst
#     pairwise CIEDE2000 now: normal 41.1, protan 12.3 (blue/purple), deutan 24.9,
#     tritan 11.1 (orange/purple), grayscale 8.0 -- all backed by the shapes below;
#   * shapes  = a REDUNDANT study cue (circle / triangle / square), so a reader who
#     cannot separate the colours can still tell the studies apart.
# NOTE: the graphical abstract at src/3 visualizations/graphical_abstract_boxplot.R
# `study_colors` still uses an older mapping.
#
# IMPORTANT: always index these maps BY STUDY NAME (imic_study_cols[studies_present]),
# never assign colours by study-presence order -- order-based assignment is what
# made a study's colour change between panels (e.g. Fig 6C has only 2 studies).
imic_study_cols <- c(
  "ELICIT"     = "#0072B2",  # Okabe-Ito blue
  "MISAME-III" = "#E69F00",  # Okabe-Ito orange
  "Mumta-LW"   = "#CC79A7"   # Okabe-Ito reddish purple
)
imic_study_shapes <- c(
  "ELICIT"     = 16,         # filled circle
  "MISAME-III" = 17,         # filled triangle
  "Mumta-LW"   = 15          # filled square
)

# Categorical palette for NON-study encodings (the milk-component categories that
# colour the FDR-significant points in the Fig 3A / Fig 5A volcanoes). Okabe-Ito minus
# black (used for non-significant points) and yellow (too faint on white), plus Tol's wine
# (#882255; replaced a brown on 2026-09-24 that matched vermillion under protanopia, dE 5.4).
# Worst pairwise CIEDE2000 across normal/protan/deutan/tritan vision is 11.1 for all seven.
# That is only just above the ~10 "distinguishable" line, so these are ALWAYS paired
# with imic_cat_shapes, one symbol per category, as a colour-independent cue; the non-significant
# (circle 16) and nominal-only (square 15) tiers keep their own shapes. Shape 25 takes
# its colour from `fill`, so volcano panels map fill as well as colour.
imic_cat_cols   <- c("#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9", "#882255")
imic_cat_shapes <- c(17, 18, 25, 8, 4, 3, 10)   # triangle, diamond, down-triangle, asterisk, x, plus, circle-plus

# Shape map aligned to a colour map's keys: each study key gets its study shape and
# any other key (e.g. "Not Significant") a filled circle. Pass the SAME keys to
# scale_colour_manual() and scale_shape_manual() (same name, same breaks) so ggplot
# merges the two into one legend.
imic_shapes_for <- function(keys) {
  s <- unname(imic_study_shapes[keys])
  s[is.na(s)] <- 16
  stats::setNames(s, keys)
}
