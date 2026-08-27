# render_msea_panelB.R -- reusable signed-enrichment-ratio MSEA volcano renderer.
#
# One function for the tertiary Fig 5B and untargeted Fig 6A. Draws the full
# pathway table: non-significant pathways form the grey cloud; significant ones
# are coloured by study and labelled. The submitted panels colour by NOMINAL
# significance (raw p < 0.05), all points filled, with a P<0.05 (grey) line and an
# upper green line (Q<0.05 for 5B, P<0.01 for 6A), and a black vertical line at 0.
#
# Input CSV columns (the full, unfiltered pathway table from run-*-msea):
#   study, timepoint, contrast, direction, pathway, total, expected, hits,
#   raw_p, fdr_native, enrichment_ratio   (enrichment_ratio already signed)
suppressMessages({ library(dplyr); library(ggplot2); library(ggrepel) })
source("figure-scripts/manuscript_figures/study_colors.R")   # canonical study colours (match graphical abstract)
source("figure-scripts/0_figure-functions.R")                # theme_imic() = the Science-submission theme
                                                             # (Helvetica, base_size 9, Reviewer-2 font floors)

render_msea_panelB <- function(msea_csv, out_png,
                               title      = NULL,
                               min_size   = 1,     # drop pathways with < min_size members
                               alpha      = 0.05,  # nominal significance line
                               er_ref     = 2,     # enrichment-ratio reference line (non-submitted)
                               submitted_style = FALSE, # colour by NOMINAL significance, all filled,
                               #   no "Sig before FDR" open tier (matches submitted 5B/6A)
                               p_hi       = 0.01,  # upper green P line when upper_line = "phi"
                               upper_line = "fdr", # "fdr" = Q-frontier green line; "phi" = P < p_hi
                               label_which= "fdr", # "fdr" = label FDR-sig; "nominal" = label raw p<alpha
                               label_max  = Inf,   # cap on labelled pathways (FDR-sig prioritised, then top -logP)
                               label_allow = NULL, # optional explicit allow-list of pathway names to label
                               #   (case-insensitive match on the pathway column). When set, ONLY these
                               #   pathways are labelled and label_max/label_which selection are bypassed;
                               #   a listed pathway with no significant cell still gets its single most-
                               #   significant cell labelled, so the allow-list is always represented.
                               vline0     = FALSE, # TRUE = black vline at 0; FALSE = blue dashed at er_ref
                               # Text/point sizes. Defaults sized for a half-page (~3.6 in) sub-panel
                               # at theme_imic base_size 9, so geom text reads at ~7-8 pt in print.
                               point_size = 2, label_size = 2.5,
                               nominal_shape = 1,  # ggplot shape for the nominal-only ("Sig before FDR")
                               #   tier: 1 = open circle (default), 16 = closed/filled circle
                               label_padding = 0.12, box_padding = 0.4, base_size = 9,
                               study_cols = imic_study_cols,  # canonical map by default; a caller may override
                               abbr_fun   = function(x) x,     # pathway-name abbreviator (submitted 5B used short names)
                               show_legend = TRUE,  # FALSE = suppress this panel's own legend (e.g. when the
                               #   caller draws one shared legend for a multi-panel composite instead)
                               xlab = "Enrichment Ratio (signed by direction)",  # caller override for the x-axis title
                               x_margin_hi = 1.05,  # right-edge headroom multiplier on x_hi; bump this if a
                               #   labelled point near the right edge gets its label box clipped by the panel
                               y_max = NULL,  # caller override for the y-axis ceiling (-log10 P scale), with
                               #   integer breaks 0:floor(y_max); NULL (default) keeps the old per-panel
                               #   auto-scaled range (max(logP)*1.08). Set for Fig 6A so it shares a common
                               #   y-axis height with Panels B/C, matching the submitted figure (2026-08-26).
                               # Science figure sizing (0_figure-functions.R): numbered figures are
                               # 2-column/full-page (7.25 in); these MSEA volcanoes are half-page
                               # sub-panels, so default ~3.6 x 3.5 in. Callers override per composite slot.
                               width_in = 3.6, height_in = 3.5) {
  supp <- read.csv(msea_csv, check.names = FALSE, stringsAsFactors = FALSE)
  supp <- supp[!is.na(supp$total) & supp$total >= min_size, , drop = FALSE]

  tableau20 <- c(
    "#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
    "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC")

  set.seed(123)  # reproducible jitter for overlapping non-significant points
  plot_df <- supp %>%
    mutate(
      study_label = recode(study, Elicit = "ELICIT", Misame = "MISAME-III",
                           Vital = "Mumta-LW", .default = study),
      logP        = -log10(raw_p),
      fdr_sig     = !is.na(fdr_native) & fdr_native < alpha,       # survives FDR
      is_sig      = raw_p < alpha,                                 # nominally significant
      nom_before  = is_sig & !fdr_sig,
      point_color = dplyr::case_when(fdr_sig ~ study_label,
                                     nom_before ~ "Sig before FDR",
                                     TRUE ~ "Not Significant"),
      pt_shape    = if_else(nom_before, "Sig before FDR", "other"),
      jitter_x = if_else(is_sig, enrichment_ratio,
                         enrichment_ratio + runif(n(), -0.6, 0.6)),
      jitter_y = if_else(is_sig, logP, pmax(0, logP + runif(n(), -0.15, 0.15))))

  # Submitted style: colour by NOMINAL significance, all points filled, drop the
  # "Sig before FDR" open tier.
  if (isTRUE(submitted_style)) {
    plot_df$point_color <- ifelse(plot_df$is_sig, plot_df$study_label, "Not Significant")
    plot_df$pt_shape    <- "other"
  }

  sig_mask <- if (isTRUE(submitted_style)) plot_df$is_sig else plot_df$fdr_sig
  studies_present <- sort(unique(plot_df$study_label[sig_mask]))
  # Colour each study by the CANONICAL named map (never by presence order), so a
  # study keeps the same colour across panels / regardless of which studies appear.
  color_vals <- c(study_cols[studies_present],
                  "Sig before FDR" = "#6BAED6", "Not Significant" = "grey75")

  x_hi <- max(c(plot_df$jitter_x, er_ref + 1), na.rm = TRUE) * x_margin_hi
  x_lo <- min(c(plot_df$jitter_x, 0), na.rm = TRUE); x_lo <- if (x_lo < 0) x_lo * 1.1 - 0.5 else -0.5
  y_hi <- if (!is.null(y_max)) y_max else max(plot_df$logP, na.rm = TRUE) * 1.08

  # Labels. If an explicit allow-list is supplied it wins: label ONLY those pathways
  # (case-insensitive), one row per study x pathway at its most-significant cell.
  # Otherwise fall back to FDR-significant, or ALL nominally-significant (submitted 5B).
  if (!is.null(label_allow)) {
    allow_lc   <- tolower(trimws(label_allow))
    is_allowed <- function(p) tolower(trimws(p)) %in% allow_lc
    # ONE label per allowed pathway, at its single most-significant cell across all
    # studies/timepoints (group_by pathway only, not study x pathway). This avoids the
    # 3x over-labelling when a pathway is significant in several study cells; the point
    # colour + "(timepoint)" tag still show which study/when. The chosen cell is the
    # strongest one, which for these pathways is also the discussed cell (the amino-acid/
    # energy pathways -> their ELICIT/MISAME down cell; the ETC/Ketone-body/Phytanic
    # reversal pathways -> MISAME-III 3-4 mo. UP). Fall back to any-cell if a listed
    # pathway is never nominally significant, so the allow-list is always represented.
    lab_sig <- plot_df %>% filter(is_sig, hits >= 1, is_allowed(pathway)) %>%
      group_by(pathway) %>% slice_max(logP, n = 1, with_ties = FALSE) %>% ungroup()
    seen      <- tolower(trimws(unique(lab_sig$pathway)))
    lab_extra <- plot_df %>% filter(hits >= 1, is_allowed(pathway),
                                    !tolower(trimws(pathway)) %in% seen) %>%
      group_by(pathway) %>% slice_max(logP, n = 1, with_ties = FALSE) %>% ungroup()
    label_df <- bind_rows(lab_sig, lab_extra)
  } else {
    label_src <- if (identical(label_which, "nominal"))
                   plot_df %>% filter(is_sig, hits >= 1)
                 else if (identical(label_which, "phi"))
                   plot_df %>% filter(raw_p < p_hi, hits >= 1)   # submitted 6A: label P < 0.01
                 else
                   plot_df %>% filter(fdr_native < alpha, hits >= 1)
    label_df <- label_src %>%
      group_by(study_label, pathway) %>% slice_max(logP, n = 1, with_ties = FALSE) %>% ungroup()
    # Thin to label_max labels: keep FDR-significant first, then the most significant
    # nominal ones, so a dense nominal set stays legible.
    if (is.finite(label_max) && nrow(label_df) > label_max)
      label_df <- label_df %>% arrange(desc(fdr_sig), desc(logP)) %>% slice_head(n = label_max)
  }
  # Submitted panels append an abbreviated timepoint, e.g. "Warburg Effect (1.5 m)".
  .abbr_tp <- function(tp) trimws(gsub("\\s*days?", " d", gsub("\\s*mo\\.?", " m", tp)))
  label_df$plab <- if (isTRUE(submitted_style))
    paste0(abbr_fun(label_df$pathway), " (", .abbr_tp(label_df$timepoint), ")") else abbr_fun(label_df$pathway)

  # Upper (green) reference line: FDR frontier (Q<alpha) or a fixed P<p_hi line.
  if (identical(upper_line, "phi")) {
    y_hi_line <- -log10(p_hi)
    upper_lab <- paste0("P < ", p_hi)
  } else {
    fdr_p_thr <- suppressWarnings(max(supp$raw_p[supp$fdr_native < alpha], na.rm = TRUE))
    y_hi_line <- if (is.finite(fdr_p_thr)) -log10(fdr_p_thr) else NA_real_
    upper_lab <- "Q < 0.05"
  }
  col_breaks <- c(studies_present,
                  if (!isTRUE(submitted_style)) "Sig before FDR",
                  "Not Significant")

  p <- ggplot(plot_df, aes(x = jitter_x, y = jitter_y)) +
    geom_point(aes(color = point_color, shape = pt_shape), size = point_size, alpha = 0.85) +
    geom_hline(yintercept = -log10(alpha), linetype = "dashed", color = "grey55") +  # P < 0.05
    { if (is.finite(y_hi_line)) geom_hline(yintercept = y_hi_line, linetype = "dashed", color = "#59A14F") } +
    { if (isTRUE(submitted_style)) annotate("text", x = -Inf, y = -log10(alpha), label = "P < 0.05", hjust = -0.05, vjust = -0.4, size = 2.4, color = "grey45") } +
    { if (isTRUE(submitted_style) && is.finite(y_hi_line)) {
        # Fig 6A (upper_line = "phi") clusters its labels at the LEFT on this line,
        # so anchor the threshold label at the RIGHT edge to keep it clear; Fig 5B
        # ("fdr") keeps the original left placement.
        .up_x  <- if (identical(upper_line, "phi")) Inf  else -Inf
        .up_hj <- if (identical(upper_line, "phi")) 1.05 else -0.05
        annotate("text", x = .up_x, y = y_hi_line, label = upper_lab,
                 hjust = .up_hj, vjust = -0.4, size = 2.4, color = "#3C8C3C")
      } } +
    { if (isTRUE(vline0)) geom_vline(xintercept = 0, color = "black", linewidth = 0.4)
      else geom_vline(xintercept = er_ref, linetype = "dashed", color = "blue") } +
    geom_label_repel(data = label_df, aes(label = plab, color = point_color),
                     size = label_size, label.padding = label_padding, box.padding = box_padding,
                     min.segment.length = 0, max.overlaps = 200, show.legend = FALSE) +
    scale_color_manual(values = color_vals, name = "Study",
                       breaks = col_breaks) +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +   # wrap so it fits the narrow half-page panel
    scale_shape_manual(values = c("Sig before FDR" = nominal_shape, "other" = 16), guide = "none") +
    scale_x_continuous(limits = c(x_lo, x_hi)) +
    scale_y_continuous(limits = c(0, y_hi),
                       breaks = if (!is.null(y_max)) seq(0, floor(y_max), 1) else waiver()) +
    labs(x = xlab,
         y = expression(-Log[10]*"(Raw P)"), title = title) +
    theme_imic(base_size = base_size) +   # Science-submission theme (Helvetica, font floors)
    theme(legend.position = if (isTRUE(show_legend)) "bottom" else "none",
          legend.key.size = unit(0.35, "cm"),
          # the volcano needs its L-shaped axes back (theme_imic drops the theme_bw border)
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3),
          plot.margin = margin(4, 6, 2, 4))

  ggsave(filename = out_png, plot = p, width = width_in, height = height_in,
         units = "in", dpi = 300, device = ragg::agg_png)
  cat("wrote", out_png, "| pathways:", nrow(plot_df),
      "| nominally sig:", sum(plot_df$is_sig),
      "| FDR-sig:", sum(plot_df$fdr_native < alpha, na.rm = TRUE), "\n")
  invisible(p)
}
