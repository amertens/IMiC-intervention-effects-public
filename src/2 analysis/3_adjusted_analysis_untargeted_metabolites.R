

#https://www.bioconductor.org/packages/devel/bioc/vignettes/biotmle/inst/doc/exposureBiomarkers.html
#https://joss.theoj.org/papers/10.21105/joss.00295

rm(list=ls())
source(paste0(here::here(),"/src/0-config.R"))

# Load the data.table package
library(data.table)

# Use fread to read the large CSV file
# Replace "path/to/your/large_file.csv" with the actual path to your file
misame <- fread(paste0(here::here(),"/data/clean milk data/MISAME/M_metabolite_s.csv")) %>% mutate(study="Misame")
gc()
colnames(misame)[1] <- "bmid"
misame[,1:20]

vital <- fread(paste0(here::here(),"/data/clean milk data/VITAL/V_metabolite_s.csv")) %>% mutate(study="Vital")
elicit <- fread(paste0(here::here(),"/data/clean milk data/ELICIT/E_metabolite_s.csv")) %>% mutate(study="Elicit")
colnames(vital)[1] <- "bmid"
colnames(elicit)[1] <- "bmid"

gc()

vital$bmid[1:30]
elicit$bmid[1:30]


#metabolomics <- bind_rows(misame[,1:1000], vital[,1:1000], elicit[,1:1000])
metabolomics <- bind_rows(misame, vital, elicit)
# saveRDS(metabolomics, paste0(here::here(),"/data/clean milk data/untargeted_metabolite_subset.RDS"))
# 
# metabolomics <- readRDS(paste0(here::here(),"/data/clean milk data/untargeted_metabolite_subset.RDS"))
# colnames(metabolomics)[1] <- "bmid"

metabolomics <- metabolomics %>% select("study","bmid","visit", everything())
Yvars <- colnames(metabolomics)[-c(1:3)] 

#------------------------------------------------------------------------------
# Load imic covariates data and merge
#------------------------------------------------------------------------------

d <- readRDS(paste0(here::here(),"/data/merged_analysis_datasets.RDS")) %>%
  select(study, visit, subjid, subjido, bmid,  all_of(Wvars)) 

# d$bmid[d$study=="Elicit"][1:30]
# d$bmid[d$study=="Vital"][1:30]
# 
# dim(d)
# dim(metabolomics)
# test <- anti_join(d, metabolomics, by=c("study","bmid","visit"))
# test2 <- anti_join(metabolomics, d, by=c("study","bmid","visit"))
# dim(test)
# dim(test2)

dim(d)
dim(metabolomics)
d$visit <- as.character(d$visit)
metabolomics$visit <- as.character(metabolomics$visit)
d_metabolomics <- left_join(d, metabolomics, by=c("study","bmid","visit")) 
dim(d_metabolomics)

save(d_metabolomics, Yvars, file=paste0(here::here(),"/data/clean milk data/untargeted_metabolites.RData"))

#------------------------------------------------------------------------------
# run analysis 
#------------------------------------------------------------------------------

SL.lib  = c("SL.glm")


res <- d_metabolomics %>% group_by(study, visit) %>%
  do(res=try(run_bioTMLE(d=.,  Wvars = Wvars, bppar.debug=T, g_lib = SL.lib, Q_lib = SL.lib,
                     Yvars=Yvars,
                     scale = TRUE,
                     bppar.type = BiocParallel::SnowParam())))
names(res$res) <- paste0(res$study, "-", res$visit)

saveRDS(res, file=paste0(here::here(),"/large-file-results/metabalomics_intervention_effects_results.RDS"))


res$res[4]
res

