# =============================================================================
# src/2 analysis/3_adjusted_analysis_proteomics.R
#
# Untargeted milk proteomics (DIA-NN protein-level output): converts proteins
# with heavy missingness to detection indicators, merges with the covariates,
# pools trial arms into Control, BEP and Nico, and fits adjusted biotmle
# intervention effects (GLM-only library) for each protein by study and visit
# in Misame and Vital. clean_results.R turns the output into the proteomics
# clean results used by 55-proteomics-go-uniprot.R (Misame cells of Fig 6C and
# Table S7) and the cross-compartment scripts (Table S8).
#
# Inputs:  data/milk/PBL_{MISAME,VITAL_L,ELICIT}_DIANN_Protein.csv
#          data/merged_analysis_datasets.RDS
# Outputs: data/clean milk data/merged_proteomics.RData
#          results/proteomics_intervention_effects_results_combined_arms.RDS
# [needs restricted data]
# =============================================================================

# Method references: biotmle vignette
# https://www.bioconductor.org/packages/devel/bioc/vignettes/biotmle/inst/doc/exposureBiomarkers.html
# and https://joss.theoj.org/papers/10.21105/joss.00295

rm(list=ls())
source(paste0(here::here(),"/src/0-config.R"))

#-------------------------------------------------------------------------------
# Load and clean data
#-------------------------------------------------------------------------------

misame <- read.csv(here("data/milk/PBL_MISAME_DIANN_Protein.csv")) %>% mutate(study="Misame")
# Vital (Mumta-LW) also uses the protein-level DIA-NN output (UniProt protein
# groups), so its features match Misame and Elicit and map to GO terms; the
# peptide-precursor-level Vital file does not.
vital <- read.csv(here("data/milk/PBL_VITAL_L_DIANN_Protein.csv")) %>% mutate(study="Vital")
elicit <- read.csv(here("data/milk/PBL_ELICIT_DIANN_Protein.csv")) %>% mutate(study="Elicit")

colnames(misame)[1] <- "bmid"
colnames(vital)[1] <- "bmid"
colnames(elicit)[1] <- "bmid"


elicit[,-1] <- convert_high_missing_to_indicators(elicit[,-1])
vital[,-1] <- convert_high_missing_to_indicators(vital[,-1])
misame[,-1] <- convert_high_missing_to_indicators(misame[,-1])


#clean vital bmid's
vital$bmid <- str_extract(vital$bmid, "A\\d+-LW")

proteomics <- bind_rows(misame%>% mutate(study="Misame"), 
                        vital%>% mutate(study="Vital"), 
                        elicit%>% mutate(study="Elicit"))

proteomics <- proteomics %>% select("study","bmid", everything())

Yvars <- colnames(proteomics)[-c(1:2)] 

#------------------------------------------------------------------------------
# Load imic covariates data and merge
#------------------------------------------------------------------------------

d <- readRDS(paste0(here::here(),"/data/merged_analysis_datasets.RDS")) %>%
  select(study, subjid, subjido, visit, bmid,  all_of(Wvars)) %>% distinct()

# Pool trial arms by the nutritional supplement received during lactation:
# Misame BEP/BEP and IFA/BEP -> BEP, BEP/IFA (prenatal BEP only) -> Control;
# Vital BEP arms -> BEP; Elicit Nico+Az. -> Nico, Az. (azithromycin only) -> Control.
table(d$arm)
table(d$arm, d$visit)
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

dim(d)
dim(proteomics)
d$bmid[1:100]
proteomics$bmid[1:100]

d_proteomics <- right_join(d, proteomics, by=c("study","bmid")) 
dim(d_proteomics)

table(d_proteomics$study, d_proteomics$arm)
table(d_proteomics$study, d_proteomics$visit)

summary(proteomics$A0A075B6I0)
summary(d_proteomics$A0A075B6I0)

class(proteomics$A0A075B6I0)
glimpse(d_proteomics)

save(d_proteomics, file=paste0(here::here(),"/data/clean milk data/merged_proteomics.RData"))

#------------------------------------------------------------------------------
# run analysis 
#------------------------------------------------------------------------------

# GLM-only library, as for the untargeted metabolomics
# (3_adjusted_analysis_untargeted_metabolites*.R); it is more stable at the small
# per-cell proteome sample sizes.
SL.lib  = c("SL.glm")

#Check for missingness in adjustment covariates.
missing_W <- d_proteomics %>% select(all_of(Wvars)) %>% summarise_all(funs(sum(is.na(.))))
missing_W      

# Drop proteomics samples with no matching participant record: after the right
# join above they have no arm (possibly laboratory control samples), so they
# cannot enter the arm contrast. Elicit is excluded below because its
# proteomics samples do not cover both arms; only Misame and Vital are modeled.
d_proteomics <- d_proteomics %>% filter(!is.na(arm))

res <- d_proteomics %>% filter(study!="Elicit") %>%
  group_by(study, visit) %>%
  do(res=try(run_bioTMLE(d=.,  Wvars = Wvars, bppar.debug=T, g_lib = SL.lib, Q_lib = SL.lib,
                     Yvars=Yvars,
                     scale = TRUE,
                     bppar.type = BiocParallel::SnowParam())))
names(res$res) <- paste0(res$study, "-", res$visit)
res$study

saveRDS(res, file=paste0(here::here(),"/results/proteomics_intervention_effects_results_combined_arms.RDS"))
