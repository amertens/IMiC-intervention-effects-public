# =============================================================================
# src/3 visualizations/sig_heatmaps.R
#
# Reads:  results/adjusted_combined_arms_intervention_effects_results_clean.RDS
# Writes: figures/perc_sig_heatmap.png
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
res <- readRDS(file=paste0(here::here(),"/results/adjusted_combined_arms_intervention_effects_results_clean.RDS"))

primary_res <- res %>% filter(outcome_group=="primary", measure=="ATE")
primary_res <- res %>% filter(measure=="ATE")
#primary_res <- res %>% filter(outcome_group!="tertiary", measure=="ATE")

#To do: add in microbiome, metabolomics, and proteomics


#-------------------------------------------------------------------------------
# Heatmap to plot percent significant by study, time, and HM group
#-------------------------------------------------------------------------------
head(primary_res)
plotdf <- primary_res %>% group_by(study, visit, outcome_group, category) %>%
  summarize(perc_sig=mean(chi_pval<0.05)*100)

plotdf$pval_cat <- cut(plotdf$perc_sig, breaks = c(-1, 10, 20, 40, 60,80, 100), 
                       labels = c("<10%", "10-20%", "20-40%", "40-60%","60-80%",">80%"))

#order axes:
unique(plotdf$study)
plotdf$study <- factor(plotdf$study, levels=c("Misame", "Vital", "Elicit"))

unique(plotdf$visit)
plotdf$visit <- factor(plotdf$visit, levels=c("1 mo.","5 mo.", "14-21 days" ,"1-2 mo.","3-4 mo.","1.5 mo.", "2 mo."    ))

unique(plotdf$category[plotdf$outcome_group=="primary"])
unique(plotdf$category[plotdf$outcome_group=="secondary"])
dput(unique(plotdf$category[plotdf$outcome_group=="tertiary"]))

plotdf$category <- factor(plotdf$category, levels=rev(c(
   "B1","B2","B3 or related","B6","Other B vitamins",
  "Macronutrient","Micronutrient","Bioactive","Individual HMO containing Fucose",
   "Individual HMO containing Sialic Acid", "Other individual HMO", 
  "Acylcarnitines", "Amino Acid Related", "Amino Acids", "Bile Acids", 
  "Biogenic Amines", "Carboxylic Acids", "Ceramides", "Cholesteryl Esters", 
  "Diglycerides", "Dihexosylceramides", "Dihydroceramides", "Fatty Acids", 
  "Hexosylceramides", "Indoles and Derivatives", "Lysophosphatidylcholines", 
  "Other metabolomics", "Phosphatidylcholines", "Sphingomyelins", 
  "Triglycerides", "Trihexosylceramides"
)))

textcol = "grey20"
hm <- ggplot(plotdf, aes(x = visit, y = category, fill = pval_cat)) + 
  geom_tile(colour = "grey80", size = 0.25) + 
  facet_wrap(outcome_group~study, nrow=3, scales = "free") +
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_discrete(expand = c(0, 0)) + 
  theme_minimal(base_size = 10) + 
  scale_fill_manual(values = brewer.pal(n = 6, name = "Oranges"), drop = FALSE) + 
  theme(aspect.ratio = 1, legend.title = element_text(color = textcol, size = 8), 
        legend.margin = margin(grid::unit(0.1, "cm")), 
        legend.text = element_text(colour = textcol, size = 7, face = "bold"), 
        legend.key.height = grid::unit(0.2, "cm"), 
        legend.key.width = grid::unit(1, "cm"), 
        legend.position = "right", 
        #axis.text.x = element_text(size = 8, colour = textcol, angle=35, vjust=0.5), 
        axis.text.x = element_text(size = 8, colour = textcol), 
        axis.text.y = element_text(size = 8, vjust = 0.2, colour = textcol), 
        axis.ticks = element_line(size = 0.4), 
        plot.title = element_text(colour = textcol, hjust = 0, size = 12, face = "bold"), 
        strip.text.x = element_text(size = 10), 
        strip.text.y = element_text(angle = 0, size = 10), 
        plot.background = element_blank(), panel.border = element_blank(), 
        strip.background = element_blank(), 
        panel.background = element_rect(fill = "grey80", colour = "grey80"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) + 
  guides(fill = guide_legend("Percent significant", ncol = 1)) + 
  labs(x = "Study and time", y = "HM Category", 
       title = "Percent significant by group")
hm

ggsave(hm, file=paste0(here::here(),"/figures/perc_sig_heatmap.png"), width=14, height=10)

#-------------------------------------------------------------------------------
# Function for heatmap for individual ATEs
# XXX   Need to fix   XXX
#-------------------------------------------------------------------------------

plot_sig_heatmap <- function (d, pval_var = "Pval", title = "", Outcome = "Outcome", 
          Exposure = "Exposure", print.est = T, print.ci = F, null = 0){
  require(RColorBrewer)
  colnames(d)[colnames(d) == pval_var] <- "pval"
  dfull <- expand_grid(unique(d$Y), unique(d$X))
  colnames(dfull) <- c("Y", "X")
  d <- left_join(dfull, d, by = c("Y", "X"))
  d <- distinct(d)
  if (null == 0) {
    d$sign <- sign(d$point.diff)
  }
  else {
    d$sign <- ifelse(d$point.diff > 1, 1, -1)
  }
  d$pval_cat <- cut(d$pval, breaks = c(-1, 0.01, 0.05, 0.2, 
                                       0.5, 2), labels = c("<0.01", "<0.05", "0.05-0.2", "0.2-0.5", 
                                                           "0.5-1"))
  d$pval_cat <- ifelse(d$sign == 1, paste0(d$pval_cat, " increase"), 
                       paste0(d$pval_cat, " decrease"))
  d$pval_cat[d$pval_cat %in% c("0.5-1 decrease", "0.5-1 increase")] <- "0.5-1"
  table(d$pval_cat)
  d$pval_cat <- factor(d$pval_cat, levels = c("<0.01 decrease", 
                                              "<0.05 decrease", "0.05-0.2 decrease", "0.2-0.5 decrease", 
                                              "0.5-1", "0.05-0.2 increase", "0.2-0.5 increase", "<0.05 increase", 
                                              "<0.01 increase"))
  d$pval_cat <- addNA(d$pval_cat)
  levels(d$pval_cat) = c(levels(d$pval_cat), "Not estimated")
  d$pval_cat[is.na(d$pval_cat)] <- "Not estimated"
  table(d$pval_cat)
  table(is.na(d$pval_cat))
  d$est = ""
  if (print.est) {
    d$est = round(d$point.diff, 2)
    if (print.ci) {
      d$est = paste0(round(d$est, 2), " (", round(d$lb.diff, 
                                                  2), ", ", round(d$ub.diff, 2), ")")
    }
  }
  d$est = gsub("NA \\(NA, NA\\)", "", d$est)
  textcol = "grey20"
  cols = rev(brewer.pal(n = 9, name = "Spectral"))
  colours <- c(`<0.01 decrease` = cols[1], `<0.05 decrease` = cols[2], 
               `0.05-0.2 decrease` = cols[3], `0.2-0.5 decrease` = cols[4], 
               `0.5-1` = cols[5], `0.2-0.5 increase` = cols[6], `0.05-0.2 increase` = cols[7], 
               `<0.05 increase` = cols[8], `<0.01 increase` = cols[9], 
               `Not estimated` = "gray80")
  d <- droplevels(d)
  hm <- ggplot(d, aes(x = X, y = Y, fill = pval_cat)) + geom_tile(colour = "grey80", 
                                                                  size = 0.25) + scale_x_discrete(expand = c(0, 0), limits = rev(levels(d$X))) + 
    scale_y_discrete(expand = c(0, 0)) + theme_minimal(base_size = 10) + 
    scale_fill_manual(values = colours, drop = FALSE) + geom_text(aes(label = est)) + 
    theme(aspect.ratio = 1, legend.title = element_text(color = textcol, 
                                                        size = 8), legend.margin = margin(grid::unit(0.1, 
                                                                                                     "cm")), legend.text = element_text(colour = textcol, 
                                                                                                                                        size = 7, face = "bold"), legend.key.height = grid::unit(0.2, 
                                                                                                                                                                                                 "cm"), legend.key.width = grid::unit(1, "cm"), legend.position = "right", 
          axis.text.x = element_text(size = 8, colour = textcol), 
          axis.text.y = element_text(size = 8, vjust = 0.2, 
                                     colour = textcol), axis.ticks = element_line(size = 0.4), 
          plot.title = element_text(colour = textcol, hjust = 0, 
                                    size = 12, face = "bold"), strip.text.x = element_text(size = 10), 
          strip.text.y = element_text(angle = 0, size = 10), 
          plot.background = element_blank(), panel.border = element_blank(), 
          strip.background = element_blank(), panel.background = element_rect(fill = "grey80", 
                                                                              colour = "grey80"), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) + guides(fill = guide_legend("P-value strength", 
                                                                           ncol = 1)) + labs(x = Exposure, y = Outcome, title = title)
  hm
  return(hm)
}



plot_sig_heatmap(primary_res)
