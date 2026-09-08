# ---------------------------------------------------------------------------
# Auto-converted from the corresponding .Rmd (knitr::purl). Run via Rscript.
# Sets the working directory to figure-scripts/manuscript_figures/ so the
# script's '../../' relative paths resolve to the repo root.
# ---------------------------------------------------------------------------
suppressPackageStartupMessages(library(here))
setwd(file.path(here::here(), "figure-scripts/manuscript_figures"))


library(tidyverse)
library(cowplot)
library(magick)   
library(tinytex)
library(tidygraph)
library(ggraph)
library(ggrepel)

# Science aesthetic: load shared palette (tableau10), theme_imic() floors, and
# save_figure_3way() (PDF + EPS + PNG export for the typesetter).
source(file.path(here::here(), "figure-scripts/0_figure-functions.R"))

knitr::opts_chunk$set(
  echo       = FALSE,
  warning    = FALSE,
  message    = FALSE,
  fig.width  = 8.27,         # should be exactly A4
  fig.height = 11.69,
  out.width  = "\\paperwidth",
  out.height = "\\paperheight",
  fig.align  = "center",
  fig.pos    = "!h"          
)

# ── automatic page break after any chunk with 'pagebreak=TRUE' ───────────
knitr::knit_hooks$set(
  pagebreak = function(before, options, envir) {
    if (!before) return("\\newpage")  # insert after the chunk runs
  }
)

# ── Page dimensions ──────────────────────────────────────────────────────
# The SUBMITTED Fig 3 is portrait (W/H = 0.682): a full-width volcano facet grid
# (~2/3 height) over a left-justified ~square pathway Panel B (~1/3 height). We
# reproduce that PROPORTION (theme_imic fonts kept) rather than the Science 7.25x9.5.
page_w <- 7.25                                # keep 2-col width
page_h <- round(page_w / 0.682, 2)            # 10.63 in -> submitted portrait aspect

#blank plots
blank_plot <- ggplot() + theme_void()

# ── helper: wrap a PNG as a ggdraw() canvas ─────────────────────────────
png_to_ggdraw <- function(path){
  ggdraw() + draw_image(image_read(path))
}

tableau10 <- c("#1F77B4","#FF7F0E","#2CA02C","#D62728",
               "#9467BD","#8C564B","#E377C2","#7F7F7F","#BCBD22","#17BECF")






# ------------------------------------------------------------
# get_bh_cutoff()
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


  make_inset_with_labels <- function(df_sum, df_long2) {
    
    xmax_val <- max(df_sum$n_total) * 1.05  # Pad x-range slightly for labels
    
   p = ggplot(df_long2) +
      # Draw bars (no spacing)
      geom_rect(
        aes(xmin = 0, xmax = count,
            ymin = ymin, ymax = ymax,
            fill = category, alpha = type),
        color = "black", size = 0.15
      ) +
      scale_fill_manual(values = final_color_palette) +
      scale_alpha_manual(values = c(n_total = 0.15, n_sig = 1)) +
      
      # # Left label: % within category
      # geom_text(
      #   data = df_sum,
      #   aes(
      #     x = 2, 
      #     y = cat_num,
      #     #label = paste0(round(100 * pct_within), "%")
      #     label = round(100 * pct_within,1)
      #   ),
      #   hjust = 0.5, size = 2
      # ) +
      # 
      # # Right label: % of total significant
      # geom_text(
      #   data = df_sum,
      #   aes(
      #     x = xmax_val-2,
      #     y = cat_num,
      #     #label = paste0(round(100 * pct_total), "%")
      #     label = round(100 * pct_total,1)
      #   ),
      #   hjust = 0.5, size = 2
      # ) +
      
      scale_y_continuous(expand = c(0,0)) +
      scale_x_continuous(expand = c(0,0), limits = c(0, xmax_val)) +
      coord_flip(clip = "off") +
      
      theme_void(base_size = 6) +
      theme(
        panel.background = element_rect(fill = "white", color = NA),
        legend.position = "none",
        plot.margin = margin(0,0,0,0)
      )
   
   # Wrap in ggdraw to add a border
  inset_with_border <- p + annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = NA, color = "black", size = 0.1)
  inset_with_border
   
  }
  

plot_imic_volcano_panel <- function(res,
                                    title              = "",
                                    label_type         = "label",
                                    n_top_vars         = 3,
                                    overlap_n          = 8) {
  
  

  
  q_cut= get_bh_cutoff(df=res,
                p_col  = "pval",
                q_col  = "pval_adj",
                alpha  = 0.05,
                return = c("raw"))

  tt_volcano <- res %>% 
    arrange(pval_adj) %>% 
    mutate(
      ATE        = est,
      logPval    = -log10(pval),
      sig_status = case_when(
        pval <= q_cut        ~ "Significant after FDR",
        pval     < 0.05        ~ "Significant before FDR",
        TRUE                         ~ "Not Significant"
      ),
      color_var  = ifelse(sig_status == "Significant after FDR",
                          category, sig_status),
      label_f    = ifelse(sig_status == "Significant after FDR", label_f, "")
    )

  
  category_summary <- tt_volcano %>%
    group_by(category) %>%
    summarise(
      n_total = n(),
      n_sig = sum(sig_status == "Significant after FDR"),
      rate_sig = n_sig / n_total
    ) %>%
    mutate(
      share_sig = n_sig / sum(n_sig)
    )
  
  
    df_sum <- category_summary %>%
    mutate(
      pct_within = n_sig / n_total,
      pct_total  = n_total / sum(n_total),
      cat_num    = as.numeric(factor(category)),
      ymin       = cat_num - 0.5,
      ymax       = cat_num + 0.5,
      non_sig    = n_total - n_sig
    )
  
  df_long2 <- df_sum %>%
    select(category, cat_num, ymin, ymax, n_sig, n_total) %>%
    pivot_longer(
      cols = c(n_sig, n_total), 
      names_to = "type", 
      values_to = "count"
    )

  
  inset <- make_inset_with_labels(df_sum, df_long2) 

  # reference lines
  p_line <- -log10(0.05)
  q_line <- -log10(q_cut)
  
  p <- ggplot(tt_volcano, aes(x = ATE, y = logPval)) +
    geom_point(aes(colour = color_var, shape=sig_status), size = 1, alpha = 0.75) +
    geom_vline(xintercept = 0,            linetype = "dashed") +
    geom_hline(yintercept = p_line,       linetype = "dashed", colour = tableau10[2]) +
    geom_hline(yintercept = q_line,       linetype = "dotted", colour = tableau10[3]) +
    xlab("") + ylab("") + ggtitle(title) +
    scale_color_manual(values = final_color_palette, na.value = "#999999") +
    guides(color = guide_legend(title = NULL)) +
    
    # scale_x_continuous(
    #   breaks = seq(floor(min(tt_volcano$ATE)),
    #                ceiling(max(tt_volcano$ATE)), by = 2)
    # ) +
    scale_x_continuous(breaks = c(-1,0,1),
                       limits = c(-1,1.5)) +
    scale_y_continuous(
      breaks = scales::breaks_pretty(n = 4),
      # ~8% headroom so the top point's label sits just below the panel border,
      # clear of the category histogram that now caps the top-right corner.
      limits = c(0, ceiling(max(tt_volcano$logPval)) * 1.08)
    ) +
    
    theme_imic(base_size = 8) +   # Science-submission theme (Helvetica, font floors)
    theme(
      legend.position = "none",
      # keep a light border so the 8 volcano facets stay delineated (theme_imic drops it)
      panel.border = element_rect(colour = "grey75", fill = NA, linewidth = 0.3),
      # Facet title now sits on the BOTTOM row of the title band (call sites pass
      # "\n\n<title>"), so the white space is ABOVE the title, just under the top-of-panel
      # histogram. A generous top plot.margin adds that breathing room (Fig 3A "more white
      # space"); the title's small bottom margin keeps it snug above the panel border.
      plot.title = element_text(size = 8, face = "bold", hjust = 0, lineheight = 0.9, margin = margin(t = 0, b = 1)),
      plot.margin = margin(t = 10, r = 0, b = -4, l = 0),
      panel.spacing = unit(0,"pt")
    )
  
  # label top variables after FDR
  top_vars <- tt_volcano %>% 
    filter(pval_adj < q_cut) %>% 
    arrange(-logPval) %>% 
    head(n = n_top_vars)
  
  p <- p + geom_text_repel(
    data            = top_vars,
    aes(label       = abbr_label(biomarker)),   # short forms so labels fit the narrow panels (Fig 3A)
    max.overlaps    = getOption("ggrepel.max.overlaps", default = overlap_n),
    size            = 2.5,
    alpha           = 0.5,
    seed            = 123)   # draw-time label placement; without it Fig 3 is not byte-reproducible
  
    # combine volcano + inset. The inset's BOTTOM (y) is placed flush on the plot's
    # top border so the histogram caps the panel with no white gap; with the taller
    # top plot.margin above, the panel top sits lower, so y is lowered to match.
  ggdraw() +
    draw_plot(p, 0, 0, 1, 1) +
    draw_plot(inset,
              x = 0.675,   # adjust horizontal location
              y = 0.782,   # histogram baseline flush on the plot's top border
              width = 0.321,
              height = 0.1725)
  
}





# TWO framings of the primary volcano grid:
#   MAIN Fig 3A  = COMBINED arms (one contrast per study x visit).
#   SUPPLEMENT   = arm-STRATIFIED (Misame/Vital per arm; Elicit has no arm split,
#                  so it uses its combined rows) -- mirrors fig5-tertiary-composite.R.
res_comb_all <- readRDS(file="../../results/combined_intervention_effects_results_combined_arms.RDS") %>%
  filter(outcome_group=="primary", measure=="ATE")
res_strat <- readRDS(file="../../results/combined_intervention_effects_results_stratified_arms.RDS") %>%
  filter(study!="Elicit", outcome_group=="primary", measure=="ATE")
res_strat_all <- bind_rows(res_strat, res_comb_all %>% filter(study=="Elicit"))

# Palette + shared x-range are built from the combined framing (all biomarkers present).
res <- res_comb_all
table(res$study,  res$contrast)


# EXPLICIT, STABLE category -> colour map. The previous version assigned
# tableau10[1:7] onto unique(res$category) in first-appearance ORDER, so any
# re-run that reordered categories reshuffled the colours (Macronutrient and
# Micronutrient swapped hues between the submission and the re-run, recolouring
# ~250 points even though every point sat in the same place). Mapping by NAME
# keeps colours identical across re-runs and matches the submitted figure
# (Macronutrient = brown, Micronutrient = pink).
category_color_palette <- c(
  "Not Significant"        = "#000000",
  "Significant before FDR" = "#000000"
)
category_colors <- c(
  "Other B vitamins" = tableau10[1],
  "B2"               = tableau10[2],
  "B3 or related"    = tableau10[3],
  "B1"               = tableau10[4],
  "B6"               = tableau10[5],
  "Macronutrient"    = tableau10[6],
  "Micronutrient"    = tableau10[7]
)
# defensive fallback for any category not in the fixed map
.unmapped <- setdiff(unique(res$category[!is.na(res$category)]), names(category_colors))
if (length(.unmapped)) category_colors <- c(category_colors, setNames(rainbow(length(.unmapped)), .unmapped))
final_color_palette <- c(category_color_palette, category_colors)


table(res$study, res$contrast)
table(res$study, res$visit)


p1 <- plot_imic_volcano_panel(res_strat_all %>% filter(study=="Elicit", visit=="1 mo.", contrast=="Nico"), title="ELICIT\n1m\nNico.")

p2 <- plot_imic_volcano_panel(res_strat_all %>% filter(study=="Misame", visit=="14-21 days", contrast=="BEP/BEP"), title="MISAME-III\n14-21d\nBEP/BEP")
p3 <- plot_imic_volcano_panel(res_strat_all %>% filter(study=="Misame", visit=="14-21 days", contrast=="BEP/IFA"), title="MISAME-III\n14-21d\nBEP/IFA")
p4 <- plot_imic_volcano_panel(res_strat_all %>% filter(study=="Misame", visit=="14-21 days", contrast=="IFA/BEP"), title="MISAME-III\n14-21d\nIFA/BEP")

p5 <- plot_imic_volcano_panel(res_strat_all %>% filter(study=="Vital", visit=="1.5 mo.", contrast=="BEP+ExBf"), title="Mumta-LW\n1.5m\nBEP")
p6 <- plot_imic_volcano_panel(res_strat_all %>% filter(study=="Vital", visit=="1.5 mo.", contrast=="BEP+ExBf+AZT"), title="Mumta-LW\n1.5m\nBEP+AZT")


p7 <- plot_imic_volcano_panel(res_strat_all %>% filter(study=="Misame", visit=="1-2 mo.", contrast=="BEP/BEP"), title="MISAME-III\n1-2m\nBEP/BEP")
p8 <- plot_imic_volcano_panel(res_strat_all %>% filter(study=="Misame", visit=="1-2 mo.", contrast=="BEP/IFA"), title="MISAME-III\n1-2m\nBEP/IFA")
p9 <- plot_imic_volcano_panel(res_strat_all %>% filter(study=="Misame", visit=="1-2 mo.", contrast=="IFA/BEP"), title="MISAME-III\n1-2m\nIFA/BEP")

p10 <- plot_imic_volcano_panel(res_strat_all %>% filter(study=="Vital", visit=="2 mo.", contrast=="BEP+ExBf"), title="Mumta-LW\n2m\nBEP")
p11 <- plot_imic_volcano_panel(res_strat_all %>% filter(study=="Vital", visit=="2 mo.", contrast=="BEP+ExBf+AZT"), title="Mumta-LW\n2m\nBEP+AZT")


p12 <- plot_imic_volcano_panel(res_strat_all %>% filter(study=="Elicit", visit=="5 mo.", contrast=="Nico"), title="ELICIT\n5m\nNico.")

p13 <- plot_imic_volcano_panel(res_strat_all %>% filter(study=="Misame", visit=="3-4 mo.", contrast=="BEP/BEP"), title="MISAME-III\n3-4m\nBEP/BEP")
p14 <- plot_imic_volcano_panel(res_strat_all %>% filter(study=="Misame", visit=="3-4 mo.", contrast=="BEP/IFA"), title="MISAME-III\n3-4m\nBEP/IFA")
p15 <- plot_imic_volcano_panel(res_strat_all %>% filter(study=="Misame", visit=="3-4 mo.", contrast=="IFA/BEP"), title="MISAME-III\n3-4m\nIFA/BEP")


# ---- MAIN (combined-arm) panels: one contrast per study x visit -------------
# Two LEADING blank rows ("\n\n<title>") place the facet title on the BOTTOM row of the
# 3-line title band, so the white space is ABOVE the title (under the top-of-panel
# histogram) rather than below it. The 3-line band height is unchanged, so the plot area
# and the inset histogram stay put.
c_m1 <- plot_imic_volcano_panel(res_comb_all %>% filter(study=="Misame", visit=="14-21 days", contrast=="BEP"),  title="\n\nMISAME-III 14-21d")
c_m2 <- plot_imic_volcano_panel(res_comb_all %>% filter(study=="Misame", visit=="1-2 mo.",    contrast=="BEP"),  title="\n\nMISAME-III 1-2m")
c_m3 <- plot_imic_volcano_panel(res_comb_all %>% filter(study=="Misame", visit=="3-4 mo.",    contrast=="BEP"),  title="\n\nMISAME-III 3-4m")
c_v1 <- plot_imic_volcano_panel(res_comb_all %>% filter(study=="Vital",  visit=="1.5 mo.",    contrast=="BEP"),  title="\n\nMumta-LW 1.5m")
c_v2 <- plot_imic_volcano_panel(res_comb_all %>% filter(study=="Vital",  visit=="2 mo.",      contrast=="BEP"),  title="\n\nMumta-LW 2m")
c_e1 <- plot_imic_volcano_panel(res_comb_all %>% filter(study=="Elicit", visit=="1 mo.",      contrast=="Nico"), title="\n\nELICIT 1m")
c_e2 <- plot_imic_volcano_panel(res_comb_all %>% filter(study=="Elicit", visit=="5 mo.",      contrast=="Nico"), title="\n\nELICIT 5m")




# Create blank plots as placeholders
blank_plot <- ggplot() + theme_void()



# Create a color legend subplot
create_category_legend <- function() {
  # Get only the category colors (excluding "Not Significant" and "Significant before FDR")
  legend_colors <- final_color_palette[!names(final_color_palette) %in% c("Not Significant", "Significant before FDR")]
  
  # Create a data frame for the legend
  legend_data <- data.frame(
    category = names(legend_colors),
    y = c(2,2.66,3.33,4,2.33,3,3.66),
    x = c(1,1,1,1,1,1,1)
  )
  
  # Create the legend plot
  legend_plot <- ggplot(legend_data, aes(x = x, y = y, fill = category)) +
    geom_point(size = 3, shape = 21, color = "black") +
    scale_fill_manual(values = legend_colors) +
    geom_text(aes(label = category), hjust = 0, nudge_x = 0.2, size = 2.5) +
    xlim(0.5, 3) +
    ylim(1.5, 4) +
    theme_void() +
    theme(legend.position = "none") +
    ggtitle("Significant Categories") +
    theme(plot.title = element_text(size = 8, hjust = 0.5),
          plot.margin = margin(b = -10)) 

  return(legend_plot)
}


# Create the legend
category_legend <- create_category_legend()


# library(cowplot)
# library(grid)
# 
# legend_gt <- ggplotGrob(category_legend)
# 
# 
# # Disable clipping in every grob that has a viewport
# for (i in seq_along(legend_gt$grobs)) {
#   g <- legend_gt$grobs[[i]]
#   if ("vp" %in% names(g)) {
#     g$vp$clip <- "off"
#     legend_gt$grobs[[i]] <- g
#   }
# }
# 
# # Wrap as a proper ggdraw object
# legend_noclip <- ggdraw() + draw_grob(legend_gt)


# Arrange plots in a 6x3 grid (18 total positions)
# Top 2/3 of page will contain plots, with last 3 positions blank for legend
plot_grid_top <- plot_grid(
  p1, p2, p3, p4, p5, p6,
  blank_plot, p7, p8, p9, p10, p11, 
  p12, p13, p14, p15, category_legend, #blank_plot,
  ncol = 6, 
  nrow = 3,
  rel_widths  = rep(1, 6),
  rel_heights = rep(1, 3),
  align = "h",
  axis = "l"
)

plot_grid_top_labeled <- ggdraw(plot_grid_top) +
  draw_label(
    "Scaled Average Treatment Effect",
    x = 0.5, y = 0.0025,
    hjust = 0.5, vjust = 0,
    size = 7
  ) +
  draw_label(
    "-log10(P-value)",
    x = 0.0025, y = 0.5,
    angle = 90,
    hjust = 0.5, vjust = 1,
    size = 7
  )


# Panel B = reproduced primary/B-vitamin MSEA (figure6-untargeted style, primary
# outcomes) from figure4-panelB-msea.R. This replaces the former grey placeholder
# (final-figure-3-b.png) — the original artwork had a broken "$"-prefixed filename.
# Panel B (KEGG pathway-impact, landscape, from fig3B-pathway.R):
#   MAIN       = COMBINED-arm pathway impact   (figures/figure3_panelB_msea.png)
#   SUPPLEMENT = arm-STRATIFIED pathway impact (figures/figure3_panelB_stratified.png)
p3b_comb_path  <- "../../figures/figure3_panelB_msea.png"
p3b_strat_path <- "../../figures/figure3_panelB_stratified.png"
p3b_comb  <- if (file.exists(p3b_comb_path))  png_to_ggdraw(p3b_comb_path)  else blank_plot
p3b_strat <- if (file.exists(p3b_strat_path)) png_to_ggdraw(p3b_strat_path) else blank_plot

# axis labels applied to whichever volcano grid is passed in
label_grid <- function(g) ggdraw(g) +
  draw_label("Scaled Average Treatment Effect", x = 0.5, y = 0.0025, hjust = 0.5, vjust = 0, size = 7) +
  draw_label("-log10(P-value)", x = 0.0025, y = 0.5, angle = 90, hjust = 0.5, vjust = 1, size = 7)

# MAIN Fig 3A: combined arms, study-major (MISAME / Mumta / ELICIT rows), matching
# the reordered Fig 5A. Legend on the 2nd row (Mumta's empty 3rd cell); the blank
# cell falls to the 3rd row.
comb_grid <- plot_grid(
  c_m1, c_m2, c_m3,
  c_v1, c_v2, category_legend,
  c_e1, c_e2, blank_plot,
  ncol = 3, nrow = 3, align = "hv", axis = "lrtb")

# Panel B sits bottom-LEFT as a ~square with blank space to its right, exactly as in
# the submitted Fig 3 (the pathway plot does not span the full page width).
b_row_comb  <- plot_grid(p3b_comb,  blank_plot, ncol = 2, rel_widths = c(0.52, 0.48))
b_row_strat <- plot_grid(p3b_strat, blank_plot, ncol = 2, rel_widths = c(0.52, 0.48))

fig3 <- plot_grid(label_grid(comb_grid), b_row_comb, labels = "AUTO",
                  ncol = 1, nrow = 2, rel_heights = c(0.68, 0.32))

# SUPPLEMENT: arm-stratified 15-panel grid (interleaved, mirrors Fig 5 supplement)
# + stratified pathway-impact Panel B (also bottom-left).
fig3_supp <- plot_grid(plot_grid_top_labeled, b_row_strat, labels = "AUTO",
                       ncol = 1, nrow = 2, rel_heights = c(0.68, 0.32))





#save MAIN as png (legacy) + 3-way (PDF/EPS/PNG) for Science submission
ggsave(
  filename = "../../figures/figure3.png",
  plot = fig3,
  width = page_w, height = page_h,
  units = "in", dpi = 300
)
save_figure_3way(fig3, name = "figure3",
                 width = page_w, height = page_h,
                 dir = file.path(here::here(), "figures"))

#save arm-stratified SUPPLEMENT (PNG + 3-way)
ggsave(
  filename = "../../figures/figure3_stratified_supplement.png",
  plot = fig3_supp,
  width = page_w, height = page_h,
  units = "in", dpi = 300
)
save_figure_3way(fig3_supp, name = "figure3_stratified_supplement",
                 width = page_w, height = page_h,
                 dir = file.path(here::here(), "figures"))


