# =============================================================================
# src/3 visualizations/2-volcano_plots-combined-arms-triglycerides.R
#
# Reads:  data/merged_analysis_datasets.RDS
#         results/adjusted_combined_arms_intervention_effects_results_clean.RDS
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

resfull <- readRDS(file=paste0(here::here(),"/results/adjusted_combined_arms_intervention_effects_results_clean.RDS"))


#better way to group: group by total number of double bonds (the 2nd numbers after the colons summed)
# the first number is not totally informative, because it isn't necessarily the longest
#could also sum the length of all 3 chains


#-------------------------------------------------------------------------------
# Clean and label triglycerides
#-------------------------------------------------------------------------------

res_tertiary <- resfull %>% filter(outcome_group=="tertiary") 

head(res_tertiary)
res_tri <- res_tertiary %>% filter(category=="Triglycerides")

head(res_tri)

# Extract the first number before the colon
res_tri$first_number <- as.numeric(sub(".*\\((\\d+):.*", "\\1", res_tri$label))

str_extract(res_tri$label, "\\((\\d+):")

# Extract the second number before the colon

res_tri <- res_tri %>%
  mutate(double_bonds = map_dbl(label, function(x) {
    # Extract all numbers after colons
    numbers <- str_extract_all(x, "(?<=:)\\d+") %>% 
      unlist() %>% 
      as.numeric()
    # Sum the numbers
    sum(numbers)
  }))
res_tri %>% select(label, double_bonds)
table(res_tri$double_bonds)

# Create a new column with the category label
res_tri$tri_category <- cut(res_tri$first_number, 
                   breaks = c(0, 6, 12, 20, Inf), 
                   labels = c("SCFA", "MCFA", "LCFA", "VLCFA"), 
                   right = TRUE)

#category with specific length
res_tri$detailed_tri_category <- paste0(res_tri$tri_category, " (C", res_tri$first_number, ")")

#clean up dataset
res_tri$visit <- factor(res_tri$visit, levels=c("14-21 days", "1-2 mo.","3-4 mo."))

res_tri$color <- ifelse(res_tri$chi_pval < 0.05, ifelse(res_tri$chi_pval_adj < 0.05, "Significant after FDR", "Significant before FDR"),"Not Significant")
res_tri$color <- factor(res_tri$color, levels = c("Significant after FDR", "Significant before FDR","Not Significant"))

res_tsm=res_tri %>% filter(study=="Misame", measure=="MN")
res=res_tri %>% filter(study=="Misame", measure=="ATE")

#-------------------------------------------------------------------------------
#  Explore Misame triglyceride wave
#-------------------------------------------------------------------------------

title="" 
facet_type="study"
labels=TRUE
label_type="label"
overlap_n=20
n_top_vars=NULL
n_top_vars_both_directions=TRUE
  
  tt_volcano <- res %>% arrange(chi_pval_adj) %>% 
    mutate(ATE = est, logPval = -log10(chi_pval), 
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
    theme_bw() + theme(legend.position="none") +  
    facet_wrap(detailed_tri_category~visit, ncol=3) 
  p
  
  p <- ggplot(tt_volcano, aes(x = ATE, y = logPval)) + geom_point(aes(colour = color), alpha=0.5) + 
    xlab("Average Treatment Effect") + ylab("-log10(raw p-value)") + 
    ggtitle(title) + 
    geom_vline(xintercept = 0, linetype="dashed") +
    # ggsci::scale_color_gsea() + 
    scale_color_manual(values = color_palette) +
    guides(color = guide_legend(title = NULL)) + 
    theme_bw() + theme(legend.position="none") +  
    facet_wrap(double_bonds~visit, ncol=3) 
  p
  
  res_hm=tt_volcano %>% group_by(visit, double_bonds) %>%
    summarise(perc_sig=mean(chi_pval<0.05)*100, mean(chi_pval))
  ggplot(res_hm, aes(x=visit, y=double_bonds , fill=perc_sig)) + geom_tile() + theme(legend.position = "bottom") + 
    scale_fill_viridis_c() + xlab("Visit") + ylab("Double bonds") + ggtitle("Percentage of significant triglycerides") 
  
  
#-------------------------------------------------------------------------------
#  Trajectory plots
#-------------------------------------------------------------------------------

#to do: plot out the treatment-specific means by category, 
#and link the points over time with lines. Color the points if the intervention effects are significant
#NOTE! I think I should do this with the raw values, not the scaled means
#Either use the raw data, or estimate the treatment-specific means without scaling

  head(res_tsm)
  res_tsm$biomarker_arm <- paste0(res_tsm$biomarker, " (", res_tsm$contrast, ")")
  plot_df <- res_tsm %>% filter(double_bonds==0) 
  p_tsm_traj <- ggplot(plot_df, aes(x = visit, y = est, group = biomarker_arm, alpha = color)) + 
    geom_line(alpha=0.1) + 
    geom_point(aes(color = color)) + 
    xlab("Visit") + ylab("Treatment-specific mean") + 
    ggtitle("") + 
    scale_color_manual(values = color_palette) +
    scale_alpha_manual(values = c(1,0.5,0.1)) +
    guides(color = guide_legend(title = NULL)) + 
    theme_bw() + 
    facet_wrap(~contrast, ncol=1) 
  p_tsm_traj
  
  
  p_tsm_traj <- ggplot(plot_df, aes(x = visit, y = est, group = contrast, alpha = color)) + 
    geom_line(alpha=0.1, position = position_dodge(0.5)) + 
    geom_point(aes(color = color), position = position_dodge(0.5)) + 
    xlab("Visit") + ylab("Treatment-specific mean") + 
    ggtitle("") + 
    scale_color_manual(values = color_palette) +
    scale_alpha_manual(values = c(1,0.5,0.1)) +
    guides(color = guide_legend(title = NULL)) + 
    theme_bw() 
  p_tsm_traj
  
  
  #load in raw data to plot the biomarkers with significant treatment effects
  d<-readRDS(paste0(here::here(),"/data/merged_analysis_datasets.RDS"))
  head(d)
  
  
  sig_biomarkers <- res %>% filter(chi_pval_adj < 0.05) %>% select(biomarker) %>% distinct()
  raw_df <- d %>% filter(study=="Misame") %>%
    select(bmid, visit, arm, all_of(str_to_lower(sig_biomarkers$biomarker))) %>% 
    mutate(arm = case_when(
        arm=="Az." ~ "Control",
        arm=="BEP+ExBf+AZT" ~ "BEP",
        arm=="BEP+ExBf" ~ "BEP",
        arm=="Nico+Az." ~ "Nico",
        arm=="BEP/BEP" ~ "BEP",
        arm=="IFA/BEP" ~ "BEP",
        arm=="BEP/IFA" ~ "Control",
        arm==arm ~ arm)) %>%
    rename(contrast=arm) %>%
    gather(key = "biomarker", value = "value", -bmid, -visit, -contrast) %>%
    left_join(res_tsm %>% filter(chi_pval_adj < 0.05) %>% 
                      mutate(biomarker   = str_to_lower(biomarker)) %>%
                distinct(biomarker, contrast, double_bonds, color), raw_df, by=c("contrast","biomarker"))
  
  raw_plot_df <- raw_df %>% group_by(biomarker, double_bonds, visit, contrast,color) %>% 
    summarise(value = mean(value, na.rm=TRUE))
  head(raw_plot_df)
  
  #line plot over visits
  
  p_tsm_traj_raw <- ggplot(raw_plot_df %>% filter(double_bonds ==2), 
                       aes(x = visit, y = value, group = contrast, color = contrast)) + 
    geom_line(position=position_dodge(0.25)) + 
    geom_point(position=position_dodge(0.25)) + #jitter?
    xlab("Visit") + ylab("Treatment-specific mean") + 
    ggtitle("") + 
    scale_color_manual(values = palette_pander(2)) +
    guides(color = guide_legend(title = NULL)) + 
    theme_void() + theme(legend.position="none") + facet_wrap(~biomarker) + 
    theme(legend.position="bottom")
  #+ facet_wrap(~contrast, ncol=1)
  p_tsm_traj_raw
  
  p_tsm_traj_raw_free_axis <- ggplot(raw_plot_df %>% filter(double_bonds ==2), 
                           aes(x = visit, y = value, group = contrast, color = contrast)) + 
    geom_line(position=position_dodge(0.25)) + 
    geom_point(position=position_dodge(0.25)) + #jitter?
    xlab("Visit") + ylab("Treatment-specific mean") + 
    ggtitle("") + 
    scale_color_manual(values = palette_pander(2)) +
    guides(color = guide_legend(title = NULL)) + 
    theme_void() + theme(legend.position="none") + facet_wrap(~biomarker, scales= "free_y") + 
    theme(legend.position="bottom")
  #+ facet_wrap(~contrast, ncol=1)
  p_tsm_traj_raw_free_axis
  
  
#NOTE! Need to look into intervention effect and see where/why the treatment effects 
#are significantly different and compare the scaled vs raw data

  res_tri %>% filter(study=="Misame", biomarker=="Tg.20.1_30.1.", measure    =="ATE") 
  df <-  d %>% filter(study=="Misame") %>% rename(biomarker="tg.20.1_30.1.") %>%
    mutate(arm = case_when(
      arm=="Az." ~ "Control",
      arm=="BEP+ExBf+AZT" ~ "BEP",
      arm=="BEP+ExBf" ~ "BEP",
      arm=="Nico+Az." ~ "Nico",
      arm=="BEP/BEP" ~ "BEP",
      arm=="IFA/BEP" ~ "BEP",
      arm=="BEP/IFA" ~ "Control",
      arm==arm ~ arm)) %>%
    rename(contrast=arm)
  table(df$contrast)
  df2 <- df %>% filter(visit==3)
  t.test(biomarker ~ contrast, data=df2)
  t.test(scale(biomarker) ~ contrast, data=df2)
  
  
  table(res_tsm$double_bonds)
  
  palette_tableau <- c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC949", "#AF7AA1", "#FF9DA7", "#9C755F", "#BAB0AB")
  traj_df <- res_tsm %>% group_by(biomarker) %>% mutate(min_sig = min(chi_pval_adj)) %>%
    filter(double_bonds ==0#, min_sig < 0.05
           ) %>% arrange(min_sig) %>%
    mutate(biomarker = str_to_lower(biomarker), biomarker=factor(biomarker, levels=unique(biomarker)))
  p_tsm_traj <- ggplot(traj_df, 
                       aes(x = studytime, y = est, group = contrast)) + 
    geom_line(aes(linetype=contrast), position=position_dodge(0.25), color = "grey30") + 
    geom_point(aes(color = color, shape=contrast), ,position=position_dodge(0.25), size=2) + 
    xlab("Visit") + ylab("Treatment-specific mean") + 
    ggtitle("") + 
    scale_color_manual(values = palette_tableau[c(2,5,1)]) +
    scale_shape_manual(values = c(19,0)) +
    guides(color = guide_legend(title = NULL)) + 
    theme_void() + theme(legend.position="bottom") + facet_wrap(~biomarker )
  p_tsm_traj
  
  #to do: arrange by the ATE not pval
  #make same plot on raw scale
  
