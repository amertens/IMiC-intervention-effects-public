# =============================================================================
# 0_figure-functions.R
#
# Shared figure helpers, sourced by the figure scripts (directly or through
# src/0-config.R): the tableau10 palette, canonical biomarker display labels,
# plotmath helpers for subscripts, the shared -log10(P) axis title, theme_imic()
# (set as the session default), Science figure dimensions, save_figure_3way()
# (PDF + EPS + PNG export) and get_bh_cutoff(). Study colours and shapes live in
# manuscript_figures/study_colors.R.
#
# Inputs: none. Outputs: none (defines objects only; save_figure_3way() writes
# the files it is given).
# =============================================================================

# tableau10: blue (1) and orange (2) mark the nominal and FDR significance tiers in
# Figs 2, S2, S3 and S4; Fig 1B colours the milk modalities with the reversed set.
tableau10 <- c("#1F77B4","#FF7F0E","#2CA02C","#D62728",
               "#9467BD","#8C564B","#E377C2","#7F7F7F","#BCBD22","#17BECF")

# ------------------------------------------------------------------------------
# Canonical biomarker display labels, in one place so every figure uses
# identical labels (Short, Title Case) with correct symbols: Greek tocopherols
# (α/γ), IgA, FGF-21, and B-vitamin subscripts. Keyed by LOWER-CASED biomarker
# code; codes not in the map fall back to the supplied label (HMOs / metabolites keep
# their label_f). Apply via canonical_label(code, fallback).
# ------------------------------------------------------------------------------
biomarker_label_map <- c(
  # macronutrients
  "kcal.l"="Total Energy","fat"="Total Fat","protein"="Total Protein",
  "cho"="Total Carbohydrate","carbohydrate"="Total Carbohydrate",
  # minerals / trace elements
  "as"="Arsenic","ca"="Calcium","ca_bio"="Calcium","cr"="Chromium","cu"="Copper",
  "fe"="Iron","k"="Potassium","mg"="Magnesium","mn"="Manganese","mo"="Molybdenum",
  "na"="Sodium","p"="Phosphorus","se"="Selenium","zn"="Zinc","choline"="Choline",
  # fat-soluble vitamins  (α- / γ-tocopherol)
  "vitamin.a"="Vitamin A","a.tocopherol"="α-Tocopherol","g.tocopherol"="γ-Tocopherol",
  # B-vitamin composites  (subscripts)
  "b1"="Vitamin B₁","b2"="Vitamin B₂","b3"="Vitamin B₃",
  "b6"="Vitamin B₆","b12"="Vitamin B₁₂",
  "pa"="Vitamin B₅","bio"="Vitamin B₇",
  # B-vitamin species
  "t"="Free Thiamin","tmp"="Thiamine Monophosphate","tpp"="Thiamine Pyrophosphate",
  "ribo"="Riboflavin","fmn"="Flavin Mononucleotide","fad"="Flavin Adenine Dinucleotide",
  "nam"="Nicotinamide","nad"="Nicotinamide Adenine Dinucleotide",
  "nmn"="Nicotinamide Mononucleotide","nr"="Nicotinamide Riboside",
  "nufa"="Nudifloramide","trp"="Tryptophan",
  "pl"="Pyridoxal","pm"="Pyridoxamine","pn"="Pyridoxine","plp"="Pyridoxal 5′-Phosphate",
  # bioactive proteins / hormones
  "fgf.21"="FGF-21","iga"="IgA","fsh"="Follicle-Stimulating Hormone",
  "lh"="Luteinizing Hormone","leptin"="Leptin","insulin"="Insulin",
  "calprotectin"="Calprotectin (S100A8/A9)"
)

canonical_label <- function(code, fallback = code) {
  code_lc  <- tolower(as.character(code))
  out      <- unname(biomarker_label_map[code_lc])
  fb       <- as.character(fallback)
  if (length(fb) == 1L) fb <- rep(fb, length(out))
  out[is.na(out)] <- fb[is.na(out)]
  out
}

# abbr_label(): short forms for the crowded primary B-vitamin volcano (Fig 3A),
# where full cofactor names overflow the narrow facet panels and get truncated at
# the panel edge. Falls back to canonical_label. Abbreviations follow the submitted
# panel (NMN, NAD, FAD, FMN, NR, TPP, PLP, ...). Case-insensitive matching.
abbr_label <- function(code) {
  lab <- canonical_label(code)
  rules <- list(
    c("nicotinic acid mononucleotide",     "NaMN"),
    c("nicotinamide mononucleotide",       "NMN"),
    c("nicotinamide adenine dinucleotide", "NAD"),
    c("flavin adenine dinucleotide",       "FAD"),
    c("flavin mononucleotide",             "FMN"),
    c("nicotinamide riboside",             "NR"),
    c("thiamine pyrophosphate",            "TPP"),
    c("thiamine monophosphate",            "TMP"),
    c("pyridoxal.*phosphate",              "PLP"))
  for (r in rules) lab <- ifelse(grepl(r[1], lab, ignore.case = TRUE), r[2], lab)
  lab
}

# plotmath_label(): Arial (what "Helvetica" resolves to on Windows) has no glyphs for
# the Unicode subscript digits U+2080-2089, so "Vitamin B₁₂" drew as "Vitamin B" plus
# an empty box in the cairo PDF/EPS files (ragg PNGs hid this by borrowing the digits
# from a fallback font). This rewrites a label as plotmath source that draws each
# subscript run as a real subscript in the plot's own font:
#   "Vitamin B₁₂" -> "Vitamin B"[12]      "Vitamin B₁ (mg/L)" -> "Vitamin B"[1]*" (mg/L)"
# Labels without subscripts come back as quoted strings, so a whole label column can
# be parsed. Draw with geom_text(_repel)(parse = TRUE), or with plotmath_expr() as a
# scale's `labels` for axis text.
plotmath_label <- function(x) {
  x <- as.character(x)
  vapply(x, function(s) {
    if (is.na(s)) return(NA_character_)
    if (!nzchar(s)) return('""')
    runs <- regmatches(s, gregexpr("[₀-₉]+|[^₀-₉]+", s, perl = TRUE))[[1]]
    out <- ""
    for (r in runs) {
      if (grepl("^[₀-₉]+$", r, perl = TRUE)) {
        digits <- chartr("₀₁₂₃₄₅₆₇₈₉",
                         "0123456789", r)
        out <- paste0(if (nzchar(out)) out else '""', "[", digits, "]")
      } else {
        q <- paste0('"', gsub('(["\\\\])', "\\\\\\1", r), '"')
        out <- if (nzchar(out)) paste0(out, "*", q) else q
      }
    }
    out
  }, character(1), USE.NAMES = FALSE)
}
plotmath_expr <- function(x) parse(text = plotmath_label(x), keep.source = FALSE)

# Shared title for every -log10(P) axis: en dash (no space), capital L, subscript 10,
# italic P, i.e. "–Log₁₀(P-value)". plotmath, so it works in labs(), ylab() and
# cowplot::draw_label() and draws in PDF/EPS too.
imic_logp_title <- expression("–"*Log[10]*"("*italic(P)*"-value)")

#-------------------------------------------------------------------------------
# theme_imic(): shared theme with minimum font sizes for print legibility:
#   axis-tick      ≥ 7 pt
#   axis-title     ≥ 8 pt
#   plot title     ≥ 10 pt, bold
#   strip text     ≥ 10 pt, bold, on a white background
#   legend text    ≥ 7 pt
# No panel border, minor gridlines, horizontal gridlines or y ticks; legend on top.
# Calling theme_imic() with a larger base_size scales everything up; the floors
# above are enforced by `max(...)`.
#-------------------------------------------------------------------------------
theme_imic <- function(base_size = 9, base_family = "Helvetica") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.border       = element_blank(),
      axis.ticks.y       = element_blank(),
      axis.text.y        = element_text(size = max(7, base_size), hjust = 0),
      axis.text.x        = element_text(size = max(7, base_size)),
      axis.title         = element_text(size = max(8, base_size + 1)),
      plot.title         = element_text(face = "bold", size = max(10, base_size + 2)),
      strip.background   = element_rect(fill = "white", colour = NA),
      strip.text         = element_text(face = "bold", size = max(10, base_size + 1)),
      legend.position    = "top",
      legend.title       = element_text(size = max(7, base_size)),
      legend.text        = element_text(size = max(7, base_size - 1))
    )
}

theme_set(theme_imic())

#-------------------------------------------------------------------------------
# Science journal figure dimensions. Recommended widths (final published size):
#   1-column:        2.24 in (57 mm)
#   1.5-column:      4.76 in (121 mm)
#   2-column / full: 7.25 in (184 mm)
# Maximum height: 9.5 in (~24 cm) per page.
# All numbered manuscript figures in this paper are 2-column (full-page).
#-------------------------------------------------------------------------------
science_dims <- list(
  one_col      = list(width = 2.24, height_max = 9.5),
  one_and_half = list(width = 4.76, height_max = 9.5),
  full_page    = list(width = 7.25, height_max = 9.5)
)

#-------------------------------------------------------------------------------
# Multi-format figure export: saves the same plot as PDF (vector, for print), EPS
# (vector, for typesetters) and PNG (raster), each into `dir`.
#-------------------------------------------------------------------------------
save_figure_3way <- function(plot, name, width = 7.25, height = 9.5, dpi = 300,
                              dir = "figures") {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  # UTF-8-capable devices so α/γ-tocopherol render (the default Windows png/postscript
  # devices fail with an mbcsToSbcs conversion error). Cairo does not fall back to another
  # font for glyphs Arial lacks, so B-vitamin subscripts must go through plotmath_label().
  devices <- list(pdf = grDevices::cairo_pdf,
                  eps = grDevices::cairo_ps,
                  png = ragg::agg_png)
  for (ext in names(devices)) {
    f <- file.path(dir, paste0(name, ".", ext))
    ggplot2::ggsave(filename = f, plot = plot,
                    width = width, height = height,
                    units = "in", dpi = dpi, bg = "white",
                    device = devices[[ext]])
  }
  invisible(file.path(dir, paste0(name, c(".pdf", ".eps", ".png"))))
}

# ------------------------------------------------------------
# get_bh_cutoff(): the BH critical raw p, i.e. the largest raw p whose BH q is
# still <= alpha; used for the Q < 0.05 line on the volcano panels.
# ------------------------------------------------------------
# df      : data-frame that already contains raw-p and BH-q columns
# p_col   : name of the raw-p column
# q_col   : name of the BH-adjusted q column
# alpha   : desired FDR level (default 0.05)
# return  : "log10"  → –log10(p_crit)   (ready for ggplot y-axis)
#           "raw"    → p_crit           (raw p threshold)
#           "both"   → list(raw = p_crit,
#                            log10 = –log10(p_crit))
# ------------------------------------------------------------
get_bh_cutoff <- function(df,
                          p_col  = "pval",
                          q_col  = "qval",
                          alpha  = 0.05,
                          return = c("log10", "raw", "both")) {

  return <- match.arg(return)

  # basic checks
  if (!all(c(p_col, q_col) %in% names(df)))
    stop("Specified p_col / q_col not found in the data frame.")

  # which tests are FDR-significant?
  sig <- df[[q_col]] <= alpha & !is.na(df[[q_col]]) & !is.na(df[[p_col]])

  if (!any(sig)) {
    warning("No q-values ≤ alpha; returning NA.")
    p_crit <- NA_real_
  } else {
    p_crit <- max(df[[p_col]][sig])
  }

  out <- switch(return,
                log10 = -log10(p_crit),
                raw   =  p_crit,
                both  =  list(raw   = p_crit,
                              log10 = -log10(p_crit)))
  out
}


