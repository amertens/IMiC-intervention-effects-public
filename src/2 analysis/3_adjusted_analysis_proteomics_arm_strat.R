

#https://www.bioconductor.org/packages/devel/bioc/vignettes/biotmle/inst/doc/exposureBiomarkers.html
#https://joss.theoj.org/papers/10.21105/joss.00295

rm(list=ls())
source(paste0(here::here(),"/src/0-config.R"))

#-------------------------------------------------------------------------------
# Load data
#-------------------------------------------------------------------------------

misame <- read.csv(here("data/milk/PBL_MISAME_DIANN_Protein.csv")) 
vital <- read.csv(here("data/milk/PBL_VITAL_L_DIANN_Protein.csv"))
elicit <- read.csv(here("data/milk/PBL_ELICIT_DIANN_Protein.csv")) 


#-------------------------------------------------------------------------------
# Clean data
#-------------------------------------------------------------------------------


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
  select(study, subjid, subjido,  bmid, visit, all_of(Wvars)) %>% distinct()



d_proteomics <- right_join(d, proteomics, by=c("study","bmid")) 
dim(d_proteomics)

table(d_proteomics$study, d_proteomics$arm)
#only in both arms in Misame and Vital

summary(proteomics$A0A075B6I0)
summary(d_proteomics$A0A075B6I0)

class(proteomics$A0A075B6I0)
glimpse(d_proteomics)

save(d_proteomics, file=paste0(here::here(),"/data/clean milk data/merged_proteomics.RData"))


#------------------------------------------------------------------------------
# run analysis 
#------------------------------------------------------------------------------

# GLM alone, consistent with the other exploratory untargeted omic (untargeted
# metabolomics) and more stable at the small per-cell proteome n. Reverted from
# the full ensemble.
SL.lib  = c("SL.glm")

#Check for missingness in adjustment covariates.
missing_W <- d_proteomics %>% select(all_of(Wvars)) %>% summarise_all(funs(sum(is.na(.))))
missing_W      


#TEMP- DEBUG unmatched (maybe lab) samples
d_proteomics <- d_proteomics %>% filter(!is.na(arm))


table(d_proteomics$study, d_proteomics$arm)
table(d_proteomics$study,  d_proteomics$arm, is.na(d_proteomics$A0A075B6H7))
table(d_proteomics$study,  d_proteomics$arm, is.na(d_proteomics$A0A075B6H9))

#number of dyads
length(unique(d_proteomics$subjid))

res <- d_proteomics %>% filter(study!="Elicit") %>%
  group_by(study, visit) %>%
  do(res=try(run_bioTMLE(d=.,  Wvars = Wvars, bppar.debug=T, g_lib = SL.lib, Q_lib = SL.lib,
                     Yvars=Yvars,
                     scale = TRUE,
                     bppar.type = BiocParallel::SnowParam())))
names(res$res) <- paste0(res$study, "-", res$visit)

#res <- res$res$Misame$res %>% mutate(study="Misame") %>% filter(measure=="ATE")


saveRDS(res, file=paste0(here::here(),"/results/proteomics_intervention_effects_results.RDS"))

