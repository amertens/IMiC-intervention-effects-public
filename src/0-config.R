# =============================================================================
# src/0-config.R
#
# Shared setup sourced by every analysis and figure script: loads the R
# packages, sources the helper functions in functions/ and
# figure-scripts/0_figure-functions.R, loads the milk-component lists by panel
# (all_milk_components: primary macro/micro/bvit, secondary hmo/protein,
# tertiary metabolomics), and defines Wvars, the adjustment set used in every
# adjusted model ("arm" first, then the baseline covariates).
#
# Inputs:  metadata/milk_component.Rdata (shipped)
# Outputs: none (objects in the calling session)
# =============================================================================

library(washb)
library(knitr)
library(biotmle)
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

# Gestational age at birth (gagebrth), maternal BMI (mbmi) and maternal MUAC
# (mmuaccm) are not adjusted for because the prenatal interventions may affect them.
