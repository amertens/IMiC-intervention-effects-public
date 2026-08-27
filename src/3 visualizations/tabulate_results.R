# =============================================================================
# src/3 visualizations/tabulate_results.R
#
# Reads:  results/adjusted_combined_arms_and_visits_intervention_effects_results_clean.RDS
#         results/adjusted_combined_arms_intervention_effects_results_clean.RDS
#         results/adjusted_intervention_effects_results_clean.RDS
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

res_pooled <- readRDS(file=paste0(here::here(),"/results/adjusted_combined_arms_and_visits_intervention_effects_results_clean.RDS"))
res <- readRDS(file=paste0(here::here(),"/results/adjusted_combined_arms_intervention_effects_results_clean.RDS"))
res_strat_arms <- readRDS(file=paste0(here::here(),"/results/adjusted_intervention_effects_results_clean.RDS"))


table(res$outcome_group, res$category)


table(res$outcome_group, res$category)

res_pooled <- extract_bioTMLE_results(res_pooled)
table(res$sig)
table(res$sigFDR)
prop.table(table(res$sig)) * 100
prop.table(table(res$sig, res$studytime)) * 100


res_primary <- res %>% filter(outcome_group=="primary", measure=="ATE") 
res_secondary <- res %>% filter(outcome_group=="secondary", measure=="ATE") 
res_tertiary <- res %>% filter(outcome_group=="tertiary", measure=="ATE") 
res_tertiary_strat <- res_strat_arms %>% filter(outcome_group=="tertiary", measure=="ATE") 


#for paper:

#In contrast, there was a consistent effect of BEP interventions on the targeted metabolites, with 

(table(res_tertiary$sig)) 
prop.table(table(res_tertiary$sig))*100
#994/3911 (25.4%) estimates impacted including 
(table(res_tertiary$sigFDR)) 
prop.table(table(res_tertiary$sigFDR)) 
#360 (9.2%) 
#significant after FDR correction (Figure 3). 
#Across all arms and trials, the bulk of significant hits 

#(n/N, X%) clustered in lipid pathways, particularly triglycerides and diglycerides.
table(res_tertiary$description[res_tertiary$sigFDR==1])
prop.table(table(res_tertiary$description[res_tertiary$sigFDR==1])) *100

#As expected from the fat sources of BEP, the largest impact was on triglycerides, the largest category of metabolites, 
#representing 43% (n/N) of all targeted metabolites, stand 75.0% of all significant hits after FDR correction. In particular, triglycerides were impacted at the 1.5 month Mumta-LW timepoint (XX significant hits after FDR correction), at the 1-month ELICIT visit (XX significant hits), and at all 3 MISAME timepoints (XX, YY and ZZ16%, 31%, and 12% significant hits, respectively). Interestingly, while the affected triglycerides were mostly increased in ELICIT and Mumta-LW, as well as at the first and third timepoint in MISAME, they were largely decreased at the second timepoint in MISAME. The same individual triglycerides were increased at 14-21 days and 3-4 months, while a different set was decreased at the 1-2 month visit (Figure S5), but when adjusted for total fat concentrations, the same set of triglycerides were increased at all timepoints (Figure S6), indicating that the apparent inversion in triglyceride responses at 1-2 months was explained by the transient decrease in total fat at this timepoint. Beyond lipids, the most affected metabolites were ceramides (n=X hits, % of total metabolites affected), free fatty acids (n, X%), and amino-acid related compounds (n, X%) (Figure 3). 
unique_bio <- res_tertiary %>% distinct(biomarker, .keep_all =TRUE)
length(unique(unique_bio$biomarker))
prop.table(table(unique_bio$category)) * 100
(table(unique_bio$category)) 

#
sig_tri <- res_tertiary %>% filter(category=="Triglycerides", sigFDR==1)
table(sig_tri$study, sig_tri$visit)


#Across trials and arms, most significant effects clustered in lipid pathways, 
#particularly triglycerides and diglycerides, accounting for 335 of 359 significant hits (93.3%).

sig_fdr <- res_tertiary %>% filter(sigFDR==1)
table(sig_fdr$category=="Triglycerides"|  sig_fdr$category=="Diglycerides" )
prop.table(table(sig_fdr$category=="Triglycerides"|  sig_fdr$category=="Diglycerides" ))*100
prop.table(table(sig_fdr$category=="Triglycerides", sig_fdr$studytime),2)*100





(table(res_tertiary$sigFDR)) 
prop.table(table(res_tertiary$category)) * 100
table(res_tertiary$sigFDR, res_tertiary$category)
prop.table(table(res_tertiary$sigFDR, res_tertiary$category),2) * 100
prop.table(table(res_tertiary$category[res_tertiary$sigFDR==1])) * 100

prop.table(table(res_primary$sigFDR, res_primary$category),2) * 100
prop.table(table(res_secondary$sigFDR, res_secondary$category),2) * 100
prop.table(table(res_tertiary$sigFDR, res_tertiary$category),2) * 100


#make heatmap, with column of percent sig by category, percent of significant that are positive, and distribution of effect sizes
#maybe volcano plots?
res_full <- bind_rows(res_primary, res_secondary, res_tertiary)
head(res_full)

table(res_full$outcome_group, res_full$category)

res_summary <- res_full %>% group_by(study, visit, studytime, category) %>% 
  summarise(n=n(), n_sig=sum(sigFDR==1), n_sig_pos=sum(sigFDR==1 & est>0), n_sig_neg=sum(sigFDR==1 & est<0), 
            n_sig_pos_prop=n_sig_pos/n_sig, n_sig_neg_prop=n_sig_neg/n_sig, 
            n_sig_prop=n_sig/n,
            n_sig_prop_total=n_sig/nrow(.), 
            mean_est=mean(est, na.rm=T), sd_est=sd(est, na.rm=T), 
            mean_est_abs=mean(abs(est), na.rm=T), sd_est_abs=sd(abs(est), na.rm=T)) %>% ungroup() %>%
  as.data.frame()
table(res_summary$category, res_summary$n)
for(i in 1:ncol(res_summary)){
  res_summary[is.nan(res_summary[,i]),i] <- NA
  #res_summary[is.na(res_summary[,i]),i] <- 0
}

p_n <- ggplot(res_summary, aes(x=studytime, y= category, fill=n)) + 
  geom_tile() + 
  scale_fill_viridis_c() + 
  theme_minimal() + coord_fixed() +
  xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + 
  ggtitle("N in category") 

#to do: organize by primary/secondary than type of category.

p_sig <- ggplot(res_summary, aes(x=studytime, y= category, fill=n_sig_prop)) + 
  geom_tile() + scale_fill_viridis_c() + 
  theme_minimal() + coord_fixed() +
  xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "none") + 
  ggtitle("Significant Effects") 

p_sig_pos <- ggplot(res_summary, aes(x=studytime, y= category, fill=n_sig_pos_prop)) + 
  geom_tile() + scale_fill_viridis_c() + 
  theme_minimal() + coord_fixed() +
  xlab("") + ylab("") +
  theme(#axis.text.y = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.position = "none") + 
  ggtitle("Significant Effects") 

p_sig_neg <- ggplot(res_summary, aes(x=studytime, y= category, fill=n_sig_neg_prop)) + 
  geom_tile() + scale_fill_viridis_c() + 
  theme_minimal() + coord_fixed() +
  xlab("") + ylab("") +
  theme(#axis.text.y = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.position = "none") + 
  ggtitle("Significant Effects") 

p_mean_est <- ggplot(res_summary, aes(x=studytime, y= category, fill=mean_est     )) + 
  geom_tile() + scale_fill_viridis_c() + 
  theme_minimal() + coord_fixed() +
  xlab("") + ylab("") +
  theme(#axis.text.y = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.position = "none") + 
  ggtitle("Mean ATE") 

p_mean_abs_est <- ggplot(res_summary, aes(x=studytime, y= category, fill=mean_est_abs      )) + 
  geom_tile() + scale_fill_viridis_c() + 
  theme_minimal() + coord_fixed() +
  xlab("") + ylab("") +
  theme(#axis.text.y = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.position = "none") + 
  ggtitle("Mean Absolute ATE") 


ggplot(res_summary %>% filter(!(category %in% c("B1","B2","B3 or related","B6","Other B vitamins"))), 
       aes(x=studytime, y= category, fill=mean_est_abs)) +
  geom_tile() + scale_fill_viridis_c() + 
  theme_minimal() + coord_fixed() +
  xlab("") + ylab("") +
  theme(#axis.text.y = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.position = "none") + 
  ggtitle("Mean Absolute ATE") 


res_full$biomarker[is.na(res_full$category)]
res_full$outcome_group[is.na(res_full$category)]
unique(res_full$biomarker[is.na(res_full$category)])
