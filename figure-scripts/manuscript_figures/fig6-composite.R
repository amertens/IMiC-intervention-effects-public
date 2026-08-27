# fig6-composite.R
# Assemble Fig 6 in TWO COLUMNS so it fits a Science print page:
#   LEFT column (stacked, half-width): A = untargeted MSEA (fig6A-untargeted-msea.R)
#                                      B = Mummichog volcano (fig6B-mummichog.R)
#                                      C = milk proteome GO ORA (fig6C-proteomics.R)
#   RIGHT column (full height):        D = cross-compartment transfer (fig6D-crosscompartment.R)
# A single-column stack was legible but ~51 cm tall (over a print page); two columns
# keep it page-sized. Panel B/A/C labels were reduced/abbreviated for the narrow
# half-width. magick trims whitespace, labels panels, and assembles.
#
# Output: figures/figure6.png (300 dpi) + figures/figure6.pdf (300-dpi raster in a
# PDF container for the typesetter). NOTE: the PDF is raster, not true vector, because
# the panels are assembled from PNGs; a true-vector Fig 6 would need the panels
# composited as ggplot objects (patchwork/cowplot) rather than images.
suppressMessages(library(magick))
FIG <- file.path(here::here(), "figures")

prep <- function(file) image_trim(image_read(file.path(FIG, file)), fuzz = 1)

# ---- LEFT column: A, B, C stacked -------------------------------------------
A <- prep("figure6_panelA_untargeted_msea.png")
B <- prep("figure6_panelB_mummichog.png")
haveC <- file.exists(file.path(FIG, "figure6_panelC_proteomics_go.png"))
left_panels <- list(A, B)
if (haveC) left_panels <- c(left_panels, list(prep("figure6_panelC_proteomics_go.png")))

WL <- max(vapply(left_panels, function(x) image_info(x)$width, numeric(1)))
left_panels <- lapply(left_panels, function(x) image_resize(x, paste0(WL, "x")))

strip <- round(WL * 0.020)
addlab <- function(img, L, w) {
  img <- image_border(img, "white", paste0("0x", strip))
  image_annotate(img, L, size = round(w * 0.030), weight = 700,
                 location = paste0("+", round(w * 0.006), "+4"))
}
left_panels <- Map(function(img, L) addlab(img, L, WL), left_panels,
                   c("A", "B", "C")[seq_along(left_panels)])
left_col <- image_append(image_join(left_panels), stack = TRUE)
# A/B/C-only stack (used by the comparison doc's "new A/B/C" slot).
image_write(left_col, file.path(FIG, "figure6_leftcol_ABC.png"))

# ---- RIGHT column: D ---------------------------------------------------------
haveD <- file.exists(file.path(FIG, "figure6_panelD_crosscompartment.png"))
if (haveD) {
  D <- prep("figure6_panelD_crosscompartment.png")
  D <- image_border(D, "white", "45x6")
  D <- addlab(D, "D", image_info(D)$width)
  # Panel D fills the ENTIRE right column: scale it to the exact height of the stacked
  # A/B/C left column (preserving D's aspect), so there is no whitespace beside the
  # tall slope plot. The left column keeps its native height.
  H <- image_info(left_col)$height
  D <- image_resize(D, paste0("x", H))
  fig <- image_append(image_join(list(left_col, D)), stack = FALSE)
} else {
  fig <- left_col
}

image_write(fig, file.path(FIG, "figure6.png"))
# Typesetter copy. Fig 6 is assembled from panel PNGs (magick), so the PDF is a high-res
# (300 dpi) RASTER container -- print-quality but not true vector. We deliberately do NOT
# write EPS: magick's raster EPS is uncompressed (~47 MB) and unusable. A true-vector Fig 6
# (and a compact EPS) would require compositing the panels as ggplot grobs (cowplot),
# the way Figs 3 and 5 are built -- a larger refactor, deferred.
image_write(fig, file.path(FIG, "figure6.pdf"), format = "pdf", density = "300")
info <- image_info(fig)
cat(sprintf("wrote figures/figure6.png + .pdf  %dx%d px  (left: A+B+%s | right: %s)\n",
            info$width, info$height,
            if (haveC) "C" else "C=PLACEHOLDER",
            if (haveD) "D" else "D=MISSING"))
