

#https://www.bioconductor.org/packages/devel/bioc/vignettes/biotmle/inst/doc/exposureBiomarkers.html
#https://joss.theoj.org/papers/10.21105/joss.00295

rm(list=ls())
source(paste0(here::here(),"/src/0-config.R"))
library(phyloseq)
library(gtools)
library(dplyr)



#------------------------------------------------------------------------------
# Merge microbiome datasets
# For Per-site analysis only - but incorporating samples present in <8000 reads, so including flag variable
## i.e - includes all taxa present in >10% prevalence Per-Site (so many taxa will Not overlap b/t sites)
#------------------------------------------------------------------------------

all_microbiome_raw <- readRDS(paste0(here::here(),"/data/microbiome/ASV_clrNorm_perSiteVisit10pct_py_ls.rds"))
all_microbiome_raw

#This can be used to get the taxonomy of the ASVs - useful for labeling/grouping the plots
metadata<-as.data.frame(tax_table(all_microbiome_raw$`CHILD.3 Months`))

# Function to make merged OTU table with selected variables from the sample data, from a list of phyloseq objects
otu_df_bind <- function(py_ls, vars, all_vars = F) {
  
  if(all_vars){
    var_select <- lapply(py_ls, function(x) colnames(sample_data(x))) %>% unlist %>% unique
  }else{
    var_select <- vars
  }
  metadata <- py_ls %>% lapply(function(x) { subset(as(sample_data(x), "data.frame"), select = var_select) } ) %>%
    bind_rows(.)
  
  otus <- py_ls %>% lapply(function(x) { data.frame(as(otu_table(x), "matrix"))} ) %>%
    bind_rows(.)
  # could also use data.frame but merging just in case 
  otu_df <- merge(metadata, otus, by = "row.names") 
  ### could add check here stopifnot nrow(otu_df) == nrow(otus)
  
  return(otu_df)
  
}

all_otu_df <- otu_df_bind(py_ls = all_microbiome_raw, vars = c("BMID", "Study", "depth_below8000")) %>% 
  rename(bmid=BMID, studyid=Study) %>% 
  mutate(bmid=gsub("lama_","",bmid), study=case_when(studyid=="MISAME" ~ "Misame", studyid=="VITAL-LW" ~ "Vital",  studyid=="ELICIT" ~ "Elicit"))
dim(all_otu_df) # 385 ASVs, 1770 samples

#table(all_otu_df$study, all_otu_df$visit)

all_otu_df$bmid[all_otu_df$study=="MISAME"] 
#IF ASV is NA, means it's present in <10% prevalence in that study/site combination. 
# Tend not to analyze these. If doing a site-visit-integrated analysis - suggest using 'ASV_clrNorm_10pctAllSiteVisit.csv' and/or AlphaDiv_Genus_ClrNorm_10pctAllSiteVisit.csv instead


Yvars <- colnames(all_otu_df %>% subset(., select = -c(Row.names, studyid, study, bmid, depth_below8000))) 

#------------------------------------------------------------------------------
# Merge in microbiome data
#------------------------------------------------------------------------------
diversity_df <- read.csv(paste0(here::here(),"/data/microbiome/AlphaDiv_Genus_ClrNorm_10pctAllSiteVisit.csv")) %>%
  select(X, Observed_imp, Shannon_imp) %>% rename(bmid=X) %>% mutate(bmid=gsub("lama_","",bmid)) 
summary(diversity_df$Observed_imp)

dim(diversity_df)
dim(all_otu_df)
d_microbiome <- left_join(diversity_df, all_otu_df, by="bmid")
dim(d_microbiome)

d_microbiome <- d_microbiome %>% filter(studyid!="CHILD") 
dim(d_microbiome)

#------------------------------------------------------------------------------
# Load imic data and merge
#------------------------------------------------------------------------------

d<-readRDS(paste0(here::here(),"/data/merged_analysis_datasets.RDS")) %>%
  select(study, visit, subjid, subjido, bmid,  all_of(Wvars)) 

d_id <- d %>% distinct(study, bmid) 
d_microbiome_id <- d_microbiome %>% distinct(study, bmid)
dim(d)
dim(d_microbiome)

#XXXXXXXXXXXXX
#NOTE! Check the failed merges
#XXXXXXXXXXXXX

test <- anti_join(d, d_microbiome, by=c("study", "bmid"))
test2 <- anti_join(d_microbiome, d, by=c("study", "bmid"))

dim(d)
dim(d_microbiome)
d <- left_join(d, d_microbiome, by=c("study", "bmid")) 
dim(d)

table(d$study,d$visit)


table(is.na(d$arm))

save(d, file=paste0(here::here(),"/data/clean milk data/merged_analysis_datasets_microbiome_arm_strat.Rdata"))

#------------------------------------------------------------------------------
# run analysis - diversity
#------------------------------------------------------------------------------

SL.lib  = c("SL.mean","SL.glm","SL.glmnet","SL.xgboost")

res_diversity<- d %>% group_by(study, visit) %>%
  do(res=try(run_bioTMLE(d=.,  Wvars = Wvars, bppar.debug=T, g_lib = SL.lib, Q_lib = SL.lib,
                         Yvars=c("Observed_imp", "Shannon_imp"),
                         scale = TRUE,
                         bppar.type = BiocParallel::SnowParam())))
names(res_diversity$res) <- paste0(res_diversity$study, "-", res_diversity$visit)
res_diversity$res
saveRDS(res_diversity, file=paste0(here::here(),"/results/microbiome_diversity_intervention_effects_results_arm_strat.RDS"))

#------------------------------------------------------------------------------
# run analysis - CLR
#------------------------------------------------------------------------------

SL.lib  = c("SL.glm")

res<- d %>% group_by(study, visit) %>%
  do(res=try(run_bioTMLE(d=.,  Wvars = Wvars, bppar.debug=T, g_lib = SL.lib, Q_lib = SL.lib,
                     Yvars=Yvars,
                     scale = TRUE,
                     bppar.type = BiocParallel::SnowParam())))
names(res$res) <- paste0(res$study, "-", res$visit)

saveRDS(res, file=paste0(here::here(),"/results/microbiome_intervention_effects_results_arm_strat.RDS"))

