# install.packages(c("tmle","dplyr","purrr"))
library(dplyr)
library(purrr)
library(tmle)

# --- Helper: fit TMLE for one biomarker & return raw p-value ---
get_pval <- function(df_subset, W_names,
                     Q.SL.library = c("SL.glm"), 
                     g.SL.library = c("SL.glm")) {
  # df_subset should have columns: value, arm, plus covariates W_names
  tm <- tmle(
    Y = df_subset$value,
    A = df_subset$arm,
    W = df_subset %>% select(all_of(W_names)),
    Q.SL.library = Q.SL.library,
    g.SL.library = g.SL.library
  )
  # two-sided Wald p-value for ATE (psi) 
  pval <- 2 * pnorm(-abs(tm$estimates$ATE$psi / tm$estimates$ATE$var.psi^0.5))
  return(pval)
}

# --- Compute observed Fisher statistic ---
compute_fisher <- function(pvals) {
  -2 * sum(log(pvals), na.rm = TRUE)
}

# --- Main permutation test function ---
global_permutation_test <- function(df_long, W_names,
                                    n_perm = 5000,
                                    seed = 2025) {
  set.seed(seed)
  
  # 1. compute observed p-values by biomarker
  obs_p <- df_long %>% 
    group_by(study, visit, biomarker) %>% 
    group_modify(~ tibble(p = get_pval(.x, W_names))) %>% 
    pull(p)
  
  T_obs <- compute_fisher(obs_p)
  
  # 2. permutation loop
  T_perm <- replicate(n_perm, {
    # permute arm within each (study,visit) stratum
    df_perm <- df_long %>%
      group_by(study, visit) %>%
      mutate(arm = sample(arm)) %>%
      ungroup()
    
    # recompute p-values
    p_perm <- df_perm %>%
      group_by(study, visit, biomarker) %>%
      group_modify(~ tibble(p = get_pval(.x, W_names))) %>%
      pull(p)
    
    compute_fisher(p_perm)
  })
  
  # 3. compute global p-value
  p_global <- (1 + sum(T_perm >= T_obs)) / (1 + n_perm)
  
  list(
    T_obs = T_obs,
    T_perm = T_perm,
    p_global = p_global
  )
}

# --- Example usage ---

# Suppose your data frame is `hmo_long` and your covariates are:
covariates <- c("sex", "mage", "meducyrs", "mhtcm", "parity", 
                "nperson","nrooms", "imp_water_src", "impfloor", 
                "cookplac", "dvseason", "dlvloc", "hhwealth", 
                "hhfoodsecure")

# Run the permutation test
res <- global_permutation_test(
  df_long   = hmo_long,
  W_names   = covariates,
  n_perm    = 5000,
  seed      = 2025
)

cat("Observed Fisher statistic: ", res$T_obs, "\n")
cat("Global p-value: ", res$p_global, "\n")
