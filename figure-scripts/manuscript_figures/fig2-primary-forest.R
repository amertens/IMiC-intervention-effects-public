# =============================================================================
# figure-scripts/manuscript_figures/fig2-primary-forest.R
#
# Reads:  results/subsetted results/primary_bvit.csv
#         results/subsetted results/primary_macro.csv
#         results/subsetted results/primary_micro.csv
#         results/subsetted results/secondary_bioactives.csv
#         results/subsetted results/secondary_hmo.csv
#         results/subsetted results/tertiary_targeted_metabolomics.csv
# Writes: figure-data/primary_forest_plots.RDS
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


tableau10 <- c("#1F77B4","#FF7F0E","#2CA02C","#D62728",
               "#9467BD","#8C564B","#E377C2","#7F7F7F","#BCBD22","#17BECF")


custom_labels <- function(labels) {
  # This assumes you're using vars(gender, education) in facet_wrap
  lapply(labels, function(x) paste(x, collapse = " | "))
}

forest_plot <- function(tab, arm_strat=F, tertiary=F){

  # canonical biomarker labels (single source of truth: 0_figure-functions.R)
  tab$label_f <- canonical_label(tab$biomarker, tab$label_f)
  tab$sigFDR <- factor(tab$sigFDR, levels=c(0,1))
  tab$sig <- factor(tab$sig, levels=c(0,1))
  # Use the final revision study names (submitted figure used these; the raw
  # keys are Vital/Elicit/Misame).
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
      # 3-tier colouring (matches the paper's volcanoes): grey = not significant,
      # blue = significant before FDR, orange = significant after FDR.
      ggplot(tab, aes(x = reorder(label_f, -est), y = est, color=sigcat, shape=sigcat, group=contrast)) +
        geom_point(position = position_dodge(width = 0.5), size = 2) +
        geom_errorbar(aes(ymin = cil, ymax = ciu), width = 0.2, position = position_dodge(width = 0.5)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
        coord_flip() +
        scale_shape_manual(values = c("Not Significant"=1, "Sig before FDR"=1, "Sig"=19),
                           drop = FALSE) +   # open circle for before-FDR (matches submission)
        scale_color_manual(values = c("Not Significant"="grey60",
                                      "Sig before FDR"=tableau10[1],
                                      "Sig"=tableau10[2]), drop = FALSE) +
        scale_x_discrete(labels = plotmath_expr) +   # B-vitamin subscripts via plotmath (Arial has no subscript-digit glyphs for PDF/EPS)
        facet_wrap(ncol=4, study~visit, labeller = labeller(.multi_line = FALSE)) +
        labs(#title = "Forest Plot of Estimates",
          x = "",
          y = "") +
        theme_bw() + theme(legend.position = "none", axis.text.y=element_text(size=7),
                           panel.grid.minor = element_blank(),
                           strip.background = element_blank(),      # white strips (match submitted)
                           strip.text = element_text(size = 7))     # theme_bw() = border around each facet
    }
  }
}






#### Macronutrients - combined arms

p1= forest_plot(read.csv(here("results/subsetted results/primary_macro.csv")))
p1

#### Micronutrients - combined arms
#Includes fat-soluble vitamins

p2= forest_plot(read.csv(here("results/subsetted results/primary_micro.csv")))

#### B-vitamins - combined arms
p3= forest_plot(read.csv(here("results/subsetted results/primary_bvit.csv")))

#### HMOs - combined arms
p4= forest_plot(read.csv(here("results/subsetted results/secondary_hmo.csv")))


#### Bioactives - combined arms
p5= forest_plot(read.csv(here("results/subsetted results/secondary_bioactives.csv")))

p6= forest_plot(read.csv(here("results/subsetted results/tertiary_targeted_metabolomics.csv")), tertiary=TRUE)
  

# Shared significance legend placed INSIDE the B-vitamin panel's empty bottom-
# right facet cell (the MISAME row has only 3 of 4 columns), matching the
# submitted figure's placement rather than a separate strip.
p3_leg <- p3 +
  theme(legend.position = c(0.995, 0.18),        # low in the empty bottom-right cell, clear of the last panel's plotted points (matches submitted)
        legend.justification = c(1, 0.5),        # right-anchored so it can't spill left
        legend.background = element_rect(fill = "white", colour = NA),
        legend.key = element_blank(),
        legend.key.size = unit(0.35, "cm"),       # keys sized to match the panel's own axis/strip text, not larger
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7)) +
  guides(colour = guide_legend(title = "Significance\n(FDR-corrected)"),
         shape  = guide_legend(title = "Significance\n(FDR-corrected)"))

fig1 <- plot_grid(
  p1 + theme(plot.margin = margin(0,0,-1,0)),
  p2 + theme(plot.margin = margin(0,0,-1,0)),
  p3_leg + theme(plot.margin = margin(0,0,0,0)),
  ncol = 1,
  labels = "AUTO",
  # Panels A/B (macro/micro) given a little more height and Panel C (22-row
  # B-vitamin) trimmed slightly, while keeping C tall enough that its labels
  # don't collide.
  rel_heights = c(0.30, 0.58, 0.80),
  align = "v",           # Align vertically
  axis = "lr")           # Align left and right axes

fig1


p1_alt= forest_plot(bind_rows(
  read.csv(here("results/subsetted results/primary_macro.csv")),
  read.csv(here("results/subsetted results/primary_micro.csv")),
  read.csv(here("results/subsetted results/primary_bvit.csv"))
  
  ))

p2_alt= forest_plot(bind_rows(
  read.csv(here("results/subsetted results/secondary_hmo.csv")),
  read.csv(here("results/subsetted results/secondary_bioactives.csv"))
))


saveRDS(list(p_forest_primary=fig1,
             p_forest_secondary=p5,
             p_forest_tertiary=p6),
        file = here("figure-data/primary_forest_plots.RDS"))

# Manuscript Fig 3 = primary-outcome forest (A: macronutrients, B: micronutrients,
# C: B-vitamins). Export PDF + EPS + PNG (Science format, white background).
save_figure_3way(fig1, name = "figure2", width = 7.25, height = 9.5)

