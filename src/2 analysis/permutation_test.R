# =============================================================================
# src/2 analysis/permutation_test.R
#
# Reads:  data/merged_analysis_datasets.RDS
#         metadata/milk_component.Rdata
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
load(file=paste0(here::here(),"/metadata/milk_component.Rdata"))
d<-readRDS(paste0(here::here(),"/data/merged_analysis_datasets.RDS"))

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
  ),
  tr=ifelse(arm=="Control", 0, 1)
)
table(d$arm)

d1 <- d %>% filter(study=="Misame", visit==1)
d2 <- d %>% filter(study=="Misame", visit==2)
d3 <- d %>% filter(study=="Misame", visit==3)
d1_2 <- d %>% filter(study=="Misame", visit %in% c(1,2)) %>% 
  group_by(study, subjid, subjido) %>% arrange(visit) %>%  
  mutate_at(vars(all_milk_components$macro, all_milk_components$micro, all_milk_components$bvit, all_milk_components$hmo, all_milk_components$protein, all_milk_components$metabolomics), funs(. - lag(.))) %>%
  filter(visit == 2) 
d2_3 <- d %>% filter(study=="Misame", visit %in% c(2,3)) %>% 
  group_by(study, subjid, subjido) %>% arrange(visit) %>%  
  mutate_at(vars(all_milk_components$macro, all_milk_components$micro, all_milk_components$bvit, all_milk_components$hmo, all_milk_components$protein, all_milk_components$metabolomics), funs(. - lag(.))) %>%
  filter(visit == 3) 



permutation_test <- function(data, outcome_vars, treatment_var, n_permutations = 1000) {
  # Calculate observed test statistic
  observed_stat <- calculate_test_statistic(data, outcome_vars, treatment_var)
  
  # Perform permutations
  perm_stats <- numeric(n_permutations)
  for (i in 1:n_permutations) {
    perm_data <- data
    perm_data[[treatment_var]] <- sample(perm_data[[treatment_var]])
    perm_stats[i] <- calculate_test_statistic(perm_data, outcome_vars, treatment_var)
  }
  
  # Calculate p-value
  p_value <- mean(perm_stats >= observed_stat)
  
  return(list(observed_stat = observed_stat, p_value = p_value))
}


data=d1
outcome_vars = c( "protein","fat","cho","kcal.l"  )
treatment_var = "tr"
n_permutations = 1000

# Helper function to calculate test statistic
calculate_test_statistic <- function(data, outcome_vars, treatment_var) {
  # Calculate standardized mean differences for each outcome
  std_diffs <- sapply(outcome_vars, function(var) {
    treat_mean <- mean(data[[var]][data[[treatment_var]] == 1])
    control_mean <- mean(data[[var]][data[[treatment_var]] == 0])
    pooled_sd <- sqrt(var(data[[var]][data[[treatment_var]] == 1]) + 
                        var(data[[var]][data[[treatment_var]] == 0])) / 2
    (treat_mean - control_mean) / pooled_sd
  })
  
  # Use sum of absolute standardized differences as test statistic
  sum(abs(std_diffs), na.rm = TRUE)
}


#NOTE! Permutation test code not currently setup to handle missing data


colnames(d1)
result <- permutation_test(d1, 
                           outcome_vars = colnames(d)[grepl("tg.", colnames(d))], 
                           treatment_var = "tr",
                           n_permutations = 1000)

print(paste("Observed test statistic:", round(result$observed_stat, 4)))
print(paste("P-value:", round(result$p_value, 4)))


result <- permutation_test(d1_2, 
                           outcome_vars = colnames(d)[grepl("tg.", colnames(d))], 
                           treatment_var = "tr",
                           n_permutations = 1000)

print(paste("Observed test statistic:", round(result$observed_stat, 4)))
print(paste("P-value:", round(result$p_value, 4)))
