# =============================================================================
# src/2 analysis/3_adjusted_analysis_microbiome.R
#
# Merges the milk microbiome data (per-study/visit CLR-normalized ASV abundances
# and genus-level alpha diversity) with the analysis covariates, pools trial
# arms into Control, BEP and Nico, and fits adjusted biotmle intervention
# effects by study and visit on (1) observed richness and Shannon diversity
# (Fig S5) and (2) each CLR-normalized ASV (GLM-only library), which
# clean_results.R adds to the combined result tables. The merged dataset it
# writes is also read by 1 data prep/4-pca-reductions.R, so run this script
# before that one.
#
# Inputs:  data/microbiome/ASV_clrNorm_perSiteVisit10pct_py_ls.rds
#          data/microbiome/AlphaDiv_Genus_ClrNorm_10pctAllSiteVisit.csv
#          data/merged_analysis_datasets.RDS
# Outputs: data/clean milk data/merged_analysis_datasets_microbiome.RDS
#          results/microbiome_diversity_intervention_effects_results.RDS (estimate tables only)
#          results/microbiome_intervention_effects_results.RDS
# [needs restricted data]
# =============================================================================

# Method references: biotmle vignette
# https://www.bioconductor.org/packages/devel/bioc/vignettes/biotmle/inst/doc/exposureBiomarkers.html
# and https://joss.theoj.org/papers/10.21105/joss.00295

rm(list=ls())
source(paste0(here::here(),"/src/0-config.R"))
library(phyloseq)
library(gtools)
library(dplyr)

#------------------------------------------------------------------------------
# Merge microbiome datasets
# Per-site ASV tables: each study/visit keeps the taxa present in >10% of its
# samples, so many taxa do not overlap between sites. Samples with <8000 reads
# are included and flagged by depth_below8000.
#------------------------------------------------------------------------------

all_microbiome_raw <- readRDS(paste0(here::here(),"/data/microbiome/ASV_clrNorm_perSiteVisit10pct_py_ls.rds"))
all_microbiome_raw

# ASV taxonomy table, for labeling (not used below)
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
  # merge on row names (sample IDs) so sample data and OTU rows stay aligned
  otu_df <- merge(metadata, otus, by = "row.names") 
  
  return(otu_df)
  
}

all_otu_df <- otu_df_bind(py_ls = all_microbiome_raw, vars = c("BMID", "Study", "depth_below8000")) %>% 
  rename(bmid=BMID, studyid=Study) %>% 
  mutate(bmid=gsub("lama_","",bmid), study=case_when(studyid=="MISAME" ~ "Misame", studyid=="VITAL-LW" ~ "Vital",  studyid=="ELICIT" ~ "Elicit"))
dim(all_otu_df) # 385 ASVs, 1770 samples

all_otu_df$bmid[all_otu_df$study=="MISAME"] 
# An NA ASV value means the taxon is below 10% prevalence in that study/visit.

Yvars <- colnames(all_otu_df %>% subset(., select = -c(Row.names, studyid, study, bmid, depth_below8000))) 

#------------------------------------------------------------------------------
# Merge in alpha diversity (imputed observed richness and Shannon index)
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

# Samples that fail to merge in either direction (inspection only)
test <- anti_join(d, d_microbiome, by=c("study", "bmid"))
test2 <- anti_join(d_microbiome, d, by=c("study", "bmid"))

dim(d)
dim(d_microbiome)
d <- left_join(d, d_microbiome, by=c("study", "bmid")) 
dim(d)

table(d$study,d$visit)

# Pool trial arms by the nutritional supplement received during lactation:
# Misame BEP/BEP and IFA/BEP -> BEP, BEP/IFA (prenatal BEP only) -> Control;
# Vital BEP arms -> BEP; Elicit Nico+Az. -> Nico, Az. (azithromycin only) -> Control.
table(d$arm)
d <- d %>% mutate(
  arm = case_when(
    arm=="Az." ~ "Control",
    arm=="BEP+ExBf+AZT" ~ "BEP",
    arm=="BEP+ExBf" ~ "BEP",
    arm=="Nico+Az." ~ "Nico",
    arm=="BEP/BEP" ~ "BEP",
    arm=="IFA/BEP" ~ "BEP",
    arm=="BEP/IFA" ~ "Control",
    arm==arm ~ arm
  )
)
table(d$arm)
d$arm <- factor(d$arm, levels=c("Control","BEP","Nico"))
levels(d$arm)
table(d$study, d$arm)
table(is.na(d$arm))

saveRDS(d, file=paste0(here::here(),"/data/clean milk data/merged_analysis_datasets_microbiome.RDS"))

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
# Keep only the estimate table: the fitted bioTMLE object embeds participant-level data.
res_diversity$res <- lapply(res_diversity$res, function(r) if (inherits(r, "try-error")) r else list(res = r$res))
saveRDS(res_diversity, file=paste0(here::here(),"/results/microbiome_diversity_intervention_effects_results.RDS"))

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

saveRDS(res, file=paste0(here::here(),"/results/microbiome_intervention_effects_results.RDS"))
