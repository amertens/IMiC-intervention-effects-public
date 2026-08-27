
rm(list=ls())
source(paste0(here::here(),"/src/0-config.R"))

# Load the data
load(file=paste0(here::here(),"/metadata/milk_component.Rdata"))
feature_importance <- read.csv("C:\\Users\\andre\\Documents\\IMiC\\imic_intervention_effects\\results\\andrew_all_site_feature_importance.csv")

head(feature_importance)
colnames(feature_importance) <- c("biomarker", "shapley", "cor", "adj.p")

# Load the study-specific data
feature_importance_ELICIT <- read.csv("C:\\Users\\andre\\Documents\\IMiC\\imic_intervention_effects\\results\\andrew_ELICIT_feature_importance.csv")
feature_importance_MISAME <- read.csv("C:\\Users\\andre\\Documents\\IMiC\\imic_intervention_effects\\results\\andrew_MISAME_feature_importance.csv")
feature_importance_VITAL <- read.csv("C:\\Users\\andre\\Documents\\IMiC\\imic_intervention_effects\\results\\andrew_VITAL_feature_importance.csv")

int_res <- readRDS(file=paste0(here::here(),"/results/adjusted_combined_arms_and_visits_intervention_effects_results.RDS"))
int_res <- extract_bioTMLE_results(int_res) %>% filter(measure=="ATE")
head(int_res)

#heatmap:
# macro nutrients
#top corresponding nutrients across group


head(feature_importance)
head(int_res)

feature_importance <- feature_importance %>% mutate(biomarker =str_to_lower(biomarker), group= "integrated") %>% rename(est=cor)
feature_importance_ELICIT <- feature_importance_ELICIT %>% mutate(biomarker =str_to_lower(Component), study="Elicit", group= "growth") %>% rename(adj.p = Adj..P, est=Spearman.rho)
feature_importance_MISAME <- feature_importance_MISAME %>% mutate(biomarker =str_to_lower(Component), study="Misame", group= "growth") %>% rename(adj.p = Adj..P, est=Spearman.rho)
feature_importance_VITAL <- feature_importance_VITAL %>% mutate(biomarker =str_to_lower(Component), study="Vital", group= "growth") %>% rename(adj.p = Adj..P, est=Spearman.rho)
int_res <- int_res %>% mutate(biomarker =str_to_lower(biomarker), group= "intervention") %>% rename(adj.p = chi_pval_adj)

dim(feature_importance_ELICIT)
dim(int_res)

dfull <- bind_rows(feature_importance_ELICIT, feature_importance_MISAME, feature_importance_VITAL, int_res)
dfull <- bind_rows(feature_importance, int_res)
d <- dfull %>% group_by(study, group) %>% arrange(adj.p) %>% 
  slice(1:200) %>% 
  #filter(adj.p<0.05) %>%
  group_by(biomarker) %>%
  mutate(anygrowth = 1*(sum(group=="integrated")>0)) %>% 
  filter(n() > 2 & anygrowth==1)
  #filter(n() > 2 & sum(group=="growth")>0)
  
head(d)
dim(d)
  
table(d$biomarker)
table(d$label)
table(d$biomarker, d$label)
  
  
plotdf <- bind_rows(
  dfull %>% filter(biomarker %in% all_milk_components$macro),
  d)

plotdf$modality <- "Impt."
plotdf$modality[plotdf$biomarker %in% all_milk_components$macro] <- "Macro"
  
ggplot(plotdf, aes(x=biomarker, y=paste0(group,"-",study), fill=est)) +
  geom_tile() + facet_wrap(~modality, scales = "free")


#to do! Get correlation between each arm and biomarker so it is on the same scale?
#mediation?
  
  
  
  
  
  