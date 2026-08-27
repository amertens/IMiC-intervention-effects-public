

#https://www.bioconductor.org/packages/devel/bioc/vignettes/biotmle/inst/doc/exposureBiomarkers.html
#https://joss.theoj.org/papers/10.21105/joss.00295

rm(list=ls())
source(paste0(here::here(),"/src/0-config.R"))

# Load the data.table package
library(data.table)

load(paste0(here::here(),"/data/clean milk data/untargeted_metabolites.RData"))

table(d_metabolomics$arm)
d_metabolomics <- d_metabolomics %>% mutate(
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
d_metabolomics$arm <- factor(d_metabolomics$arm, levels=c("Control","BEP","Nico"))
levels(d_metabolomics$arm)
table(d_metabolomics$study, d_metabolomics$arm)
table(is.na(d_metabolomics$arm))


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

saveRDS(res, file=paste0(here::here(),"/large-file-results/metabalomics_intervention_effects_results_combined_arms.RDS"))


res$res[4]
res

res =readRDS(paste0(here::here(),"/large-file-results/metabalomics_intervention_effects_results_combined_arms.RDS"))

resdf=extract_bioTMLE_results(res, single_group = TRUE)
head(resdf)

summary(resdf$est)
summary(resdf$est[resdf$sigFDR==1])

tab=resdf %>% group_by(study) %>%
  summarise(mean(sigFDR), mean(est), sd(est), median(est), IQR(est),
            mean(est[sigFDR==1]), sd(est[sigFDR==1]), median(est[sigFDR==1]), IQR(est[sigFDR==1]),
            n_sig=sum(sigFDR==1)
  )

knitr::kable(tab, digits=3)


res =readRDS(paste0(here::here(),"/large-file-results/metabalomics_intervention_effects_results_combined_arms.RDS"))

resdf=extract_bioTMLE_results(res, single_group = TRUE)
head(resdf)

summary(resdf$est)
summary(resdf$est[resdf$sigFDR==1])

tab=resdf %>% group_by(study) %>%
  summarise(mean(sigFDR), mean(est), sd(est), median(est), IQR(est),
            mean(est[sigFDR==1]), sd(est[sigFDR==1]), median(est[sigFDR==1]), IQR(est[sigFDR==1]),
            n_sig=sum(sigFDR==1)
  )

knitr::kable(tab, digits=3)



resdf <- readRDS(file=paste0(here::here(),"/results/adjusted_combined_arms_intervention_effects_results_clean.RDS"))


summary(resdf$est)
summary(resdf$est[resdf$sigFDR==1])

tab=resdf %>% group_by(study, outcome_group) %>%
  summarise(mean(sigFDR)*100, mean(est),  median(est),
            mean(est[sigFDR==1]), median(est[sigFDR==1]),
            n_sig=sum(sigFDR==1)
  )

knitr::kable(tab, digits=3)

