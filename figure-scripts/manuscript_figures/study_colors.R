# =============================================================================
# study_colors.R
#
# Shared colour and shape maps, sourced by the figure scripts. Every study-coloured
# panel (Fig 3B pathway impact, Fig 5B and 6A MSEA, Fig 5C fatty acids, Fig 6B
# mummichog, Fig 6C proteome GO) uses these named maps so a study keeps the same
# colour and symbol in every figure; the category palette colours the
# FDR-significant points of the Fig 3A and Fig 5A volcanoes.
#
# Inputs: none. Outputs: none (defines objects only).
# =============================================================================

# Study colours are Okabe-Ito blue / orange / reddish purple, chosen to stay
# distinguishable under simulated colour-vision deficiency (worst pairwise CIEDE2000:
# normal 41.1, protan 12.3 (blue/purple), deutan 24.9, tritan 11.1 (orange/purple),
# grayscale 8.0). Shapes (circle / triangle / square) are a redundant study cue, so
# the studies can be told apart without colour.
#
# Always index these maps by study name (imic_study_cols[studies_present]), never
# by the order in which studies appear: some panels have only two studies (e.g.
# Fig 6C), and order-based assignment would change a study's colour between panels.
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

# Categorical palette for non-study encodings (the milk-component categories that
# colour the FDR-significant points in the Fig 3A / Fig 5A volcanoes). Okabe-Ito minus
# black (used for non-significant points) and yellow (too faint on white), plus Tol's
# wine (#882255). Worst pairwise CIEDE2000 across normal/protan/deutan/tritan vision is
# 11.1 for all seven, only just above the ~10 "distinguishable" line, so these are
# always paired with imic_cat_shapes, one symbol per category, as a colour-independent
# cue; the non-significant (circle 16) and nominal-only (square 15) tiers keep their
# own shapes. Shape 25 takes its colour from `fill`, so volcano panels map fill as
# well as colour.
imic_cat_cols   <- c("#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9", "#882255")
imic_cat_shapes <- c(17, 18, 25, 8, 4, 3, 10)   # triangle, diamond, down-triangle, asterisk, x, plus, circle-plus

# Shape map aligned to a colour map's keys: each study key gets its study shape and
# any other key (e.g. "Not Significant") a filled circle. Pass the same keys to
# scale_colour_manual() and scale_shape_manual() (same name, same breaks) so ggplot
# merges the two into one legend.
imic_shapes_for <- function(keys) {
  s <- unname(imic_study_shapes[keys])
  s[is.na(s)] <- 16
  stats::setNames(s, keys)
}
