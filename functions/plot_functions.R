# =============================================================================
# functions/plot_functions.R
#
# plot_SLvim(): the arm-classification AUC plot (CV-AUC with 95% CI by
# predictor group, colored by collection time, faceted by study x arm), used by
# src/3 visualizations/5-SL_VIM_plots.R for the working versions of Fig 1A.
# Sourced by src/0-config.R; the tableau10 palette is defined in
# figure-scripts/0_figure-functions.R.
#
# Inputs:  none (functions only)
# Outputs: none
# =============================================================================

plot_SLvim <- function(d, legend_pos="none"){
  SLvim_lab_plot <- ggplot(d, aes(x=group, y=cvAUC, group=visit_f, color=visit_f)) +
    geom_point(position = position_dodge(width = 0.5)) +
    geom_linerange(aes(ymin=ci.lb, ymax=ci.ub), position = position_dodge(width = 0.5)) +
    geom_hline(yintercept = 0.5, linetype="dashed") +
    coord_flip() +
    scale_color_manual(values = tableau10[c(8,2,10)], drop = FALSE) +
    facet_grid(studyid~arm_f) +
    labs(color = "Collection time") +
    theme(strip.background = element_blank(),
          axis.text = element_text(size = 6),
          legend.position=legend_pos)  +
    xlab("Group of\npredictors used") + ylab("CV-AUC") 
  SLvim_lab_plot
}
