# =============================================================================
# functions/plot_functions.R
#
# This script's input/output paths are assembled at runtime, so they
# cannot be listed here without executing it.
#
# Header generated from the code itself; it makes no claim about method.
# See README.md for run order and results/ARTIFACT_MANIFEST.csv for the
# exhibit each script feeds.
# =============================================================================



plot_imic_volcano <- function(res, title="", facet_type="study",
                              labels=TRUE, label_type="label", overlap_n=20,
                              n_top_vars=NULL, n_top_vars_both_directions=TRUE){

  # Defensive: ungroup so arrange() can find columns regardless of grouping,
  # and bail out early if input is empty or lacks the required p-value columns.
  res <- dplyr::ungroup(res)
  if (nrow(res) == 0 || !all(c("pval", "pval_adj") %in% colnames(res))) {
    warning("plot_imic_volcano: input has no rows or missing pval/pval_adj — returning empty plot.")
    return(ggplot2::ggplot() + labs(title = title, subtitle = "(no data)"))
  }

  tt_volcano <- res %>% arrange(pval_adj) %>%
    mutate(ATE = est, logPval = -log10(pval), 
           color = ifelse(pval < 0.05, ifelse(pval_adj < 0.05, "Significant after FDR", "Significant before FDR"),"Not Significant"),
           label_f=ifelse(color=="Significant after FDR",label_f,""),
           color = factor(color, levels = c("Significant after FDR", "Significant before FDR","Not Significant"))) 
  
  # Define the color palette based on factor levels from cbbPalette
  color_palette <- c("Not Significant"= "#000000", "Significant before FDR" = "#56B4E9", "Significant after FDR" = "#E69F00")
  
  p <- ggplot(tt_volcano, aes(x = ATE, y = logPval)) + geom_point(aes(colour = color), alpha=0.5) + 
    xlab("Average Treatment Effect") + ylab("-log10(raw p-value)") + 
    ggtitle(title) + 
    geom_vline(xintercept = 0, linetype="dashed") +
    # ggsci::scale_color_gsea() + 
    scale_color_manual(values = color_palette) +
    guides(color = guide_legend(title = NULL)) + 
    theme_bw() + theme(legend.position="none")
  
  if(!is.null(n_top_vars)){
    # Get the top N most significant points
    
    if(n_top_vars_both_directions){
      
      #NOTE! this doesnt work well with faceted plots
      df_increase <- tt_volcano %>% filter(pval_adj < 0.05, est>0) %>%
        arrange(-logPval) %>%
        head(n=n_top_vars)
      df_decrease <- tt_volcano %>% filter(pval_adj < 0.05, est<0) %>%
        arrange(-logPval) %>%
        head(n=n_top_vars)
      top_vars <- bind_rows(df_increase, df_decrease)
    }else{
      top_vars <- tt_volcano %>% filter(pval_adj < 0.05) %>%
        arrange(-logPval) %>%
        head(n=n_top_vars)
    }

    
  }else{
    top_vars = tt_volcano %>% filter(pval_adj < 0.05)
  }
  
  if(labels==TRUE){
    if(label_type=="label"){
      p <- p + geom_text_repel(data = top_vars, aes(label=label_f), max.overlaps = getOption("ggrepel.max.overlaps", default = overlap_n), alpha=0.5)
    }else{
      if(label_type=="category"){
        p <- p + geom_text_repel(data = top_vars, aes(label=category), max.overlaps = getOption("ggrepel.max.overlaps", default = overlap_n), alpha=0.5)
      }else{
        if(label_type=="label+category"){
          p <- p + geom_text_repel(data = top_vars, aes(label=paste0(label_f,"\n",category)), max.overlaps = getOption("ggrepel.max.overlaps", default = overlap_n), alpha=0.5)
        }else{
          p <- p + geom_text_repel(data = top_vars, aes(label=biomarker), max.overlaps = getOption("ggrepel.max.overlaps", default = overlap_n), alpha=0.5)
        }
      }
    }
  }
  
  
  if(facet_type=="outcome group"){
    p <- p +  facet_grid(outcome_group~studytime, scales="free_y") 
  }
  if(facet_type=="study"){
    p <- p +  facet_wrap(~studytime) 
  }
  if(facet_type=="arm"){
    p <- p +  facet_grid(visit~contrast) 
  }
  
  return(p)
}




imic_intervention_plot_function <- function(df, title=""){
  
  df <- df %>% filter(sigFDR==1)
  
  ggplot(df, aes(x=contrast, y=est)) + geom_point() +
    geom_linerange(aes(ymin=cil, ymax=ciu)) +
    geom_hline(yintercept = 0, linetype="dashed") +
    coord_flip() +
    facet_wrap(studytime~label_f, scale="free") +
    ggtitle(title) +
    theme(strip.background = element_blank(),
          axis.text = element_text(size = 6),
          strip.text = element_text(size = 8),
          legend.position="none")  + ylab("Intervention arm") + xlab("Z-score difference")
}






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

