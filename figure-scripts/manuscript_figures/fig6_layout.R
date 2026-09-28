# =============================================================================
# fig6_layout.R
#
# Shared Fig 6 geometry, sourced by fig6A-6D and fig6-composite.R. Fig 6 is assembled
# from panel PNGs (fig6-composite.R) and printed 7.25 in wide, so a panel prints
# smaller or larger than it was rendered: Panels A-C at FIG6_ABC_SCALE of their
# rendered size and Panel D at FIG6_D_SCALE (D is scaled up to the column height,
# A-C are scaled down). fig6-composite.R reports both scales on every run and warns
# when they drift from these values; update them here if the layout changes.
# fig6_pt() turns a printed size into the size to render at, so every panel prints
# its axis text at the same sizes as Figs 3 and 5.
#
# Inputs: none. Outputs: none (defines constants and helpers only).
# =============================================================================
FIG6_Y_MAX      <- 9      # shared -log10(P) ceiling of Panels A-C
FIG6_PANEL_W_MM <- 105    # Panels A-C render size
FIG6_PANEL_H_MM <- 95
FIG6_D_W_IN     <- 4.4    # Panel D render size
FIG6_D_H_IN     <- 9.9
FIG6_ABC_SCALE  <- 0.76   # printed / rendered size (reported by fig6-composite.R)
FIG6_D_SCALE    <- 0.93
FIG6_TICK_PT    <- 7.5    # printed tick-label size
FIG6_TITLE_PT   <- 8      # printed axis-title size (Figs 3 and 5: 8 pt)
FIG6_LEGEND_PT  <- 7      # printed legend-text size
fig6_pt <- function(pt, scale) pt / scale

# Threshold-line captions ("P-value < 0.05", "Q-value < 0.05") in Panels B and C. Points
# reach both x edges there, so a fixed edge position would hide the captions.
# fig6_caption_x() scans left-anchored positions along a line and returns
# the start (hjust = 0) of the stretch whose caption box covers the fewest points. The
# captions sit on opaque white boxes, so crossing the zero line costs only half a point:
# it is used only when every other stretch covers data. fig6_caption_obstacles() puts
# invisible repel obstacles over that box so the pathway labels keep clear of it. A 6.8 pt
# caption spans about 22% of an A-C panel's width and sits 0.05-0.55 y-units above its line.
FIG6_CAPTION_W <- 0.22
.fig6_xlim <- function(x) range(x) + c(-1, 1) * 0.05 * diff(range(x))   # ggplot's default expansion
fig6_caption_x <- function(x, y, y_line, dy = c(0.05, 0.55)) {
  lim <- .fig6_xlim(x); w <- FIG6_CAPTION_W * diff(lim); pad <- 0.02 * diff(lim)
  starts <- seq(lim[1] + pad, lim[2] - w - pad, length.out = 200)
  n_cov <- vapply(starts, function(s) sum(x >= s - pad & x <= s + w + pad &
                                          y >= y_line + dy[1] & y <= y_line + dy[2]), numeric(1))
  n_cov <- n_cov + 0.5 * (starts < 0 & starts + w > 0)   # crossing the zero line
  starts[max(which(n_cov == min(n_cov)))]     # ties: the rightmost, like Panel A's captions
}
fig6_caption_obstacles <- function(x0, y_line, x) {
  data.frame(fe = seq(x0, x0 + FIG6_CAPTION_W * diff(.fig6_xlim(x)), length.out = 12),
             y = y_line + 0.3)
}
