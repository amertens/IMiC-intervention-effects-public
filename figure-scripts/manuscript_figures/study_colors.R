# study_colors.R -- SINGLE SOURCE OF TRUTH for IMiC study colours.
#
# Every study-coloured manuscript panel (Fig 3B pathway, Fig 5B/6A MSEA, Fig 5C
# fatty-acid, Fig 6B mummichog, Fig 6C proteome GO) must use THIS named map so a
# study keeps the same colour in every figure.
#   MISAME-III = orange, Mumta-LW = green, ELICIT = blue.
# (Palette set 2026-08-13 per author decision: MISAME-III orange, Mumta-LW green.
# NOTE: the graphical abstract at src/3 visualizations/graphical_abstract_boxplot.R
# `study_colors` still uses the older MISAME-green / Mumta-orange mapping and should
# be updated to match if visual consistency with the graphical abstract is wanted.)
#
# IMPORTANT: always index this map BY STUDY NAME (imic_study_cols[studies_present]),
# never assign colours by study-presence order -- order-based assignment is what
# made a study's colour change between panels (e.g. Fig 6C has only 2 studies).
imic_study_cols <- c(
  "ELICIT"     = "#1F77B5",  # blue
  "MISAME-III" = "#FE7F01",  # orange
  "Mumta-LW"   = "#2CA02C"   # green
)
