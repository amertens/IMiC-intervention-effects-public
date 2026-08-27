# =============================================================================
# src/0-config.R
#
# Reads:  metadata/milk_component.Rdata
#
# Paths above were recovered from this script's syntax tree and are
# repo-relative; they resolve from the repo root via here::here().
#
# Header generated from the code itself; it makes no claim about method.
# See README.md for run order and results/ARTIFACT_MANIFEST.csv for the
# exhibit each script feeds.
# =============================================================================



library(washb)
library(knitr)
library(biotmle)
#library(biotmleData)
library(BiocParallel)
library(SuperLearner)
library(caret)
library(SummarizedExperiment, quietly=TRUE)
library(metafor)
library(drtmle)
library(tidyverse)
library(future)
library(tlverse)
library(origami)
library(cvAUC)
library(sl3)
library(data.table)
library(ggrepel)
library(ggthemes)
library(here)


#-------------------------------------------------------------------------------
# Functions
#-------------------------------------------------------------------------------

source(paste0(here(),"/functions/SL_functions.R"))
source(paste0(here(),"/functions/plot_functions.R"))
source(paste0(here(),"/functions/bioTMLE_functions.R"))
source(paste0(here(),"/functions/data_cleaning_functions.R"))

source(paste0(here::here(),"/figure-scripts/0_figure-functions.R"))


#components by primary/secondary/tertiary
load(file=paste0(here(),"/metadata/milk_component.Rdata"))

#-------------------------------------------------------------------------------
# confounder list
#-------------------------------------------------------------------------------

Wvars = c("arm","sex","mage", "meducyrs", 
          "mhtcm", 
          "parity", "nperson","nrooms", "imp_water_src",
          "impfloor", "cookplac",
          "dvseason", "dlvloc","hhwealth", "hhfoodsecure",
          "nperson_miss","parity_miss","meducyrs_miss","hhfoodsecure_miss")

#dont adjust for gestational age at birth or maternal anthro because prenatal interventions may affect these
#"mbmi","mmuaccm" ,gagebrth
