# =============================================================================
# figure-scripts/manuscript_figures/figS3-hmo-bioactives.R
#
# Reads:  results/subsetted results/secondary_bioactives.csv
#         results/subsetted results/secondary_hmo.csv
#
# Paths above were recovered from this script's syntax tree and are
# repo-relative; they resolve from the repo root via here::here().
#
# Header generated from the code itself; it makes no claim about method.
# See README.md for run order and results/ARTIFACT_MANIFEST.csv for the
# exhibit each script feeds.
# =============================================================================



rm(list=ls())
source(paste0(here::here(),"/src/0-config.R"))
library(cowplot)

# Standalone "Significance (FDR-corrected)" legend, hand-built from a tiny synthetic
# data frame rather than extracted from either forest plot's own ggplot legend.
# Necessary because of a ggplot2 4.0.x rendering bug (confirmed 2026-08-26 via an
# isolated repro): a discrete scale level with ZERO underlying data rows -- as "Sig"
# (FDR-significant) is here, since no biomarker in this figure clears FDR -- does not
# draw its point glyph in the legend even with scale_*_manual(..., drop = FALSE); the
# text label appears but the coloured/shaped key next to it is blank. Same pattern as
# create_study_legend() in fig5-tertiary-composite.R.
create_significance_legend <- function() {
  lev <- c("Not Significant", "Sig before FDR", "Sig")
  legend_df <- data.frame(x = 1, y = 1:3, sigcat = factor(lev, levels = lev))
  p <- ggplot(legend_df, aes(x, y, colour = sigcat, shape = sigcat)) +
    geom_point(size = 2) +
    scale_shape_manual(values = c("Not Significant" = 1, "Sig before FDR" = 1, "Sig" = 19),
                       name = "Significance\n(FDR-corrected)") +
    scale_colour_manual(values = c("Not Significant" = "grey60", "Sig before FDR" = tableau10[1],
                                   "Sig" = tableau10[2]),
                        name = "Significance\n(FDR-corrected)") +
    theme_void(base_size = 8) +
    theme(legend.position = "right", legend.title = element_text(size = 8),
          legend.text = element_text(size = 7))
  comps <- cowplot::get_plot_component(p, "guide-box", return_all = TRUE)
  good  <- Filter(function(g) !inherits(g, "zeroGrob"), comps)
  if (length(good) == 0) stop("create_significance_legend(): no non-empty guide-box found")
  good[[1]]
}


tableau10 <- c("#1F77B4","#FF7F0E","#2CA02C","#D62728",
               "#9467BD","#8C564B","#E377C2","#7F7F7F","#BCBD22","#17BECF")


custom_labels <- function(labels) {
  # This assumes you're using vars(gender, education) in facet_wrap
  lapply(labels, function(x) paste(x, collapse = " | "))
}

forest_plot <- function(tab, arm_strat=F, tertiary=F, pos="none", y_limits=NULL){

  # canonical biomarker labels (single source of truth: 0_figure-functions.R)
  tab$label_f <- canonical_label(tab$biomarker, tab$label_f)
  tab$sigFDR <- factor(tab$sigFDR, levels=c(0,1))
  tab$sig <- factor(tab$sig, levels=c(0,1))
  tab$study <- dplyr::recode(as.character(tab$study),
                             "Vital"="Mumta-LW", "Elicit"="ELICIT", "Misame"="MISAME-III")
  tab$study <- factor(tab$study, levels=c("Mumta-LW","ELICIT","MISAME-III"))
  tab$visit <- factor(tab$visit, levels=c(  "1.5 mo.","2 mo.","1 mo.","5 mo.","14-21 days", "1-2 mo.","3-4 mo."))
  tab=tab %>% mutate(sigcat=case_when(
    sigFDR == 1 & sig == 1 ~ "Sig",
    sigFDR == 0 & sig == 1 ~ "Sig before FDR",
    sigFDR == 0 & sig == 0 ~ "Not Significant"
  ), sigcat=factor(sigcat, levels=c("Not Significant", "Sig before FDR", "Sig")))
  unique(tab$visit)
  unique(tab$study)
  
  if(arm_strat){
    ggplot(tab, aes(x = reorder(label_f, -est), y = est, color=contrast, shape=sigcat, group=contrast)) +
      geom_point(position = position_dodge(width = 0.5), size = 2) +
      geom_errorbar(aes(ymin = cil, ymax = ciu), width = 0.2, position = position_dodge(width = 0.5)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
      coord_flip() +
      scale_shape_manual(values = c(1, 19, 17), 
                         labels = c("Not Significant","Sig before FDR", "Sig")) +
      scale_color_manual(values = c(tableau10)) +
      facet_wrap(ncol=4, study~visit, labeller = labeller(.multi_line = FALSE)) +
      labs(title = "Forest Plot of Estimates",
           x = "Variable",
           y = "Estimate") +
      theme_minimal() + theme(legend.position = "right")
  }else{
    
    if(tertiary){
      tab=tab %>% filter(sigFDR==1) %>% droplevels()
      ggplot(tab,
             aes(x = reorder(label_f, est), y = est, color=category, group=contrast)) +
        geom_point(position = position_dodge(width = 0.5), size = 2, alpha=0.5) +
        geom_errorbar(aes(ymin = cil, ymax = ciu), width = 0.2, position = position_dodge(width = 0.5), alpha=0.5) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
        coord_flip() +
        # scale_color_manual(values = c("grey60", tableau10[2]), 
        #                    labels = c("Not Significant", "Significant")) +
        facet_wrap(ncol=4, study~visit, labeller = labeller(.multi_line = FALSE)) +
        labs(#title = "Forest Plot of Estimates",
          x = "",
          y = "") +
        theme_minimal() + theme(legend.position = "right", axis.text.y=element_blank())
    }else{
      # 3-tier scheme consistent with Fig 2 / Fig S3: grey = not significant,
      # tableau BLUE = significant before FDR, tableau ORANGE = significant after FDR
      # (open circle before FDR, filled after). theme_bw() draws a border around each
      # facet, matching the submitted forests.
      ggplot(tab, aes(x = reorder(label_f, -est), y = est, color=sigcat, shape=sigcat, group=contrast)) +
        geom_point(position = position_dodge(width = 0.5), size = 2) +
        geom_errorbar(aes(ymin = cil, ymax = ciu), width = 0.2, position = position_dodge(width = 0.5)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
        coord_flip() +
        scale_shape_manual(values = c("Not Significant"=1, "Sig before FDR"=1, "Sig"=19),
                           name="Significance\n(FDR-corrected)", drop = FALSE,
                           limits = c("Not Significant","Sig before FDR","Sig")) +  # force orange "Sig" key even with no sig data
        scale_color_manual(values = c("Not Significant"="grey60",
                                      "Sig before FDR"=tableau10[1], "Sig"=tableau10[2]),
                           name="Significance\n(FDR-corrected)", drop = FALSE,
                           limits = c("Not Significant","Sig before FDR","Sig")) +
        facet_wrap(ncol=4, study~visit, labeller = labeller(.multi_line = FALSE)) +
        { if (!is.null(y_limits)) scale_y_continuous(limits = y_limits) } +  # shared x-axis (post-flip) across A/B, 2026-08-26
        labs(x = "", y = "") +
        theme_bw() + theme(legend.position = pos, axis.text.y=element_text(size=7),
                           panel.grid.minor = element_blank(),
                           strip.background = element_blank(),   # white strips (match Fig 2 / submitted)
                           strip.text = element_text(size = 7),
                           legend.title = element_text(size = 8),
                           legend.text = element_text(size = 7))
    }
  }
}






# This script produces supplement Fig S4 = SECONDARY outcomes (HMOs + bioactive
# proteins), COMBINED-arm forest (arm_strat=FALSE, coloured by per-arm sigFDR).
# NOTE: the primary-outcome combined forest is main-text Fig 3 (fig2-primary-forest.R),
# and Fig S3 is the deficiency-RR forest (figS2-milq-deficiency.R), neither
# is produced here.

hmo_dat <- read.csv(here("results/subsetted results/secondary_hmo.csv")) %>%
  dplyr::filter(biomarker != "secretor")
bioact_dat <- read.csv(here("results/subsetted results/secondary_bioactives.csv"))
# Shared x-axis (post-coord_flip) range across Panels A and B, so both read on the
# same scale instead of each auto-scaling to its own data (reported 2026-08-26: the
# panels previously showed different tick ranges, e.g. A: -0.8/0.4, B: -0.4/0.8).
shared_ylim <- range(c(hmo_dat$cil, hmo_dat$ciu, bioact_dat$cil, bioact_dat$ciu), na.rm = TRUE)
# no manual padding here -- scale_y_continuous()'s own default expansion adds it,
# same as each panel's previous auto-scaled range did

#### HMOs - combined arms
# Drop "Maternal Secretor Status" (biomarker "secretor") to match the submitted
# Fig S4, which did not include this row.
# pos="none": both panels' own ggplot legends are suppressed -- the submitted figure
# shows only ONE legend total, and ggplot's own legend for the "Sig" (FDR-significant)
# level is broken here anyway (a ggplot2 4.0.x bug: a zero-row discrete level draws no
# key glyph even with drop=FALSE, confirmed 2026-08-26). A hand-built standalone legend
# (create_significance_legend(), below) is composited in instead.
p4= forest_plot(hmo_dat, pos="none", y_limits = shared_ylim)

#### Bioactives - combined arms
# pos="none": see p4 above -- previously this panel's legend floated at fixed coords
# (0.9, 0.55), landing inside the real ELICIT-5mo facet content and getting clipped by
# its border, hiding the orange "Sig" key (reported 2026-08-26).
p5= forest_plot(bioact_dat, pos="none", y_limits = shared_ylim)

# ---- Fig S4: secondary outcomes (HMOs + bioactive proteins), combined arms ----
# p4/p5 alignment is exactly the original call (align="v"/axis="lr" between two RAW
# ggplots) -- a first attempt that nested p4 with the legend into its own plot_grid
# BEFORE this outer align step broke p4/p5's facet/axis alignment badly (cowplot's
# align machinery expects raw ggplot cells, not a pre-composed grid, 2026-08-26).
# Instead, the legend column is built as its OWN separate nested grid (legend over a
# blank spacer, same 1:0.4 height split as stack_ab) and placed beside stack_ab with
# NO alignment requested between the two columns -- this keeps the legend vertically
# scoped to Panel A's row (matching the submitted figure) without ever asking cowplot
# to axis-align a non-ggplot legend grob against a real panel.
stack_ab <- plot_grid(
  p4+ theme(plot.margin = margin(0,0,0,0)),
  p5+ theme(plot.margin = margin(0,0,0,0)),
  ncol = 1,
  labels = "AUTO",
  rel_heights = c(1,0.4),
  align = "v",           # Align vertically
  axis = "lr"            # Align left and right axes
)
# Overlay the legend directly on the empty facet cell in Panel A's own 4-col grid
# (row 2 only has 3 of 4 columns -- MISAME-III's three visits -- leaving column 4
# empty, same slot the submitted figure's own auto-placed legend used). Placed as an
# ggdraw() overlay AFTER stack_ab is fully assembled/aligned (not nested beforehand),
# so it never touches the p4/p5 alignment step. x/y tuned to land in that gap: lower
# (row 2 of Panel A's facet grid, not vertically centered on the whole A+B stack) and
# further left (inside Panel A's own boundary, not a separate strip past its right
# edge) than the first two attempts (reported 2026-08-26).
fig_s4 <- ggdraw(stack_ab) +
  draw_plot(create_significance_legend(), x = 0.78, y = 0.30, width = 0.20, height = 0.22)

fig_s4

save_figure_3way(fig_s4, name = "figureS3_hmo_bioactives", width = 7.25, height = 10,
                 dir = paste0(here::here(), "/figures"))

