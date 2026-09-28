# =============================================================================
# fig6-composite.R
#
# Assembles Fig 6 from the four panel PNGs, in two columns so it fits a print page:
#   left column (stacked, half-width): A = untargeted MSEA (fig6A-untargeted-msea.R)
#                                      B = mummichog volcano (fig6B-mummichog.R)
#                                      C = milk proteome GO ORA (fig6C-proteomics.R)
#                                      + one shared A-C key under C
#   right column (full height):        D = cross-compartment comparison (fig6D-crosscompartment.R)
# magick trims whitespace, adds the panel letters and assembles. A missing Panel C
# or D is left out (the log says so).
#
# Inputs:  figures/figure6_panel{A_untargeted_msea,B_mummichog,C_proteomics_go,
#          D_crosscompartment}.png and figures/figure6_legend_ABC.png (panel scripts)
# Outputs: figures/figure6.png (300 dpi) and figures/figure6.pdf (the same 300-dpi
#          raster in a PDF container: the panels are PNGs, so the PDF is not vector)
# Panel D needs an on-request file (see fig6D-crosscompartment.R).
# =============================================================================
suppressMessages(library(magick))
FIG <- file.path(here::here(), "figures")
source(file.path(here::here(), "figure-scripts/manuscript_figures/fig6_layout.R"))

# Panel letters print at the size used in Figs 1-5 (cowplot's 14 pt bold) once the
# figure is printed PRINT_W_IN wide.
PRINT_W_IN <- 7.25
LETTER_PT  <- 14

prep <- function(file) image_trim(image_read(file.path(FIG, file)), fuzz = 1)

# ---- LEFT column: A, B, C stacked -------------------------------------------
A <- prep("figure6_panelA_untargeted_msea.png")
B <- prep("figure6_panelB_mummichog.png")
haveC <- file.exists(file.path(FIG, "figure6_panelC_proteomics_go.png"))
left_panels <- list(A, B)
if (haveC) left_panels <- c(left_panels, list(prep("figure6_panelC_proteomics_go.png")))

WL <- max(vapply(left_panels, function(x) image_info(x)$width, numeric(1)))
last_w <- image_info(left_panels[[length(left_panels)]])$width   # bottom panel, pre-resize
left_panels <- lapply(left_panels, function(x) image_resize(x, paste0(WL, "x")))

strip <- round(WL * 0.020)
# Each panel gets a white band on top that carries its letter; the letters themselves are
# drawn after assembly, once the final width (and so the print scale) is known.
band <- round(WL * 0.065)
addband <- function(img) image_append(c(image_blank(image_info(img)$width, band, "white"), img),
                                      stack = TRUE)
abc_scale_px <- WL / vapply(left_panels, function(x) image_info(x)$width, numeric(1))[1]
left_panels <- lapply(left_panels, addband)
panel_tops <- cumsum(c(0, head(vapply(left_panels, function(x) image_info(x)$height, numeric(1)), -1)))
# One shared key for A-C (Panel A's legend: 3 studies + Not
# Significant; written by fig6A-untargeted-msea.R) under the bottom panel. The panels
# carry no legends of their own, so all three plot areas are the same size. The key is
# scaled by the same factor as the panel above it (so its text matches the panels'),
# then centred on the column width -- never stretched to fill it.
leg_file <- "figure6_legend_ABC.png"
if (file.exists(file.path(FIG, leg_file))) {
  key <- prep(leg_file)
  key <- image_resize(key, paste0(round(image_info(key)$width * WL / last_w), "x"))
  key <- image_extent(key, paste0(WL, "x", image_info(key)$height + 2 * strip),
                      gravity = "center", color = "white")
  left_panels <- c(left_panels, list(key))
}
left_col <- image_append(image_join(left_panels), stack = TRUE)

# ---- RIGHT column: D ---------------------------------------------------------
haveD <- file.exists(file.path(FIG, "figure6_panelD_crosscompartment.png"))
if (haveD) {
  D <- prep("figure6_panelD_crosscompartment.png")
  D <- image_border(D, "white", "45x6")
  D_h0 <- image_info(D)$height
  # Panel D fills the entire right column: scale it to the exact height of the stacked
  # A/B/C left column (preserving D's aspect), so there is no whitespace beside the
  # tall plot. The left column keeps its native height.
  H <- image_info(left_col)$height
  D <- image_resize(D, paste0("x", H - band))   # addband() adds the letter band on top
  d_scale_px <- (H - band) / D_h0
  D <- addband(D)
  fig <- image_append(image_join(list(left_col, D)), stack = FALSE)
} else {
  fig <- left_col
}

# Letters, all one size: LETTER_PT at the printed width.
px_per_in <- image_info(fig)$width / PRINT_W_IN
lab_px <- round(LETTER_PT / 72 * px_per_in)
lab_at <- data.frame(L = c("A", "B", "C")[seq_along(panel_tops)], x = 0, y = panel_tops)
if (haveD) lab_at <- rbind(lab_at, data.frame(L = "D", x = image_info(left_col)$width, y = 0))
for (i in seq_len(nrow(lab_at)))
  fig <- image_annotate(fig, lab_at$L[i], size = lab_px, weight = 700, font = "Arial",
                        location = paste0("+", lab_at$x[i] + round(0.15 * lab_px), "+", lab_at$y[i]))

# Print scales of the panels (printed size / rendered size). The panel scripts size their
# text with FIG6_ABC_SCALE / FIG6_D_SCALE (fig6_layout.R) so all four panels print at the
# same sizes; warn when the layout has drifted from those constants.
abc_print <- abc_scale_px * 300 / px_per_in
d_print   <- if (haveD) d_scale_px * 300 / px_per_in else NA_real_
cat(sprintf("print scales: A-C %.3f (FIG6_ABC_SCALE %.2f), D %.3f (FIG6_D_SCALE %.2f); %.2f x %.2f in\n",
            abc_print, FIG6_ABC_SCALE, d_print, FIG6_D_SCALE,
            PRINT_W_IN, image_info(fig)$height / px_per_in))
if (abs(abc_print / FIG6_ABC_SCALE - 1) > 0.03 || (haveD && abs(d_print / FIG6_D_SCALE - 1) > 0.03))
  warning("Fig 6 print scales drifted from fig6_layout.R; update FIG6_ABC_SCALE / FIG6_D_SCALE ",
          "and re-run the panel scripts", call. = FALSE)

image_write(fig, file.path(FIG, "figure6.png"))
# Typesetter copy: a 300 dpi raster in a PDF container. No EPS is written, because
# magick's raster EPS is uncompressed (~47 MB).
image_write(fig, file.path(FIG, "figure6.pdf"), format = "pdf", density = "300")
info <- image_info(fig)
cat(sprintf("wrote figures/figure6.png + .pdf  %dx%d px  (left: A+B+%s | right: %s)\n",
            info$width, info$height,
            if (haveC) "C" else "C=PLACEHOLDER",
            if (haveD) "D" else "D=MISSING"))
