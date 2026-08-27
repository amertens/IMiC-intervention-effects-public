# =============================================================================
# 13-cross-compartment-proteome-overlap.R
#
# Cross-compartment concordance of BEP intervention effects between HUMAN MILK and
# MATERNAL BLOOD proteomes, matched on shared UniProt accession. This is the
# reviewer-facing "are the same features moving the same way across compartments?"
# analysis -- the proteomics arm (metabolomics needs m/z-RT alignment, deferred).
#
# Inputs (combined-arms, MISAME-3):
#   milk  : results/adjusted_combined_arms_intervention_effects_proteomics_results_clean_ATE.RDS
#   blood : results/blood_compartment_combined_arms_intervention_effects_results_clean.RDS
#           (datasets ProteomicsDepleted [primary] + ProteomicsNaive [sensitivity])
#
# Outputs:
#   results/cross_compartment_proteome_overlap.csv      per-protein milk vs blood ATE
#   results/cross_compartment_proteome_summary.RDS      concordance summary stats
#   figures/cross_compartment/proteome_milk_vs_blood_<depleted|naive>.png
# =============================================================================

suppressMessages({library(dplyr); library(ggplot2); library(ggrepel)})
root <- paste0(here::here(), "/")
outfig <- paste0(root, "figures/cross_compartment"); dir.create(outfig, showWarnings = FALSE, recursive = TRUE)

up <- function(x) toupper(trimws(x))   # normalize UniProt accession case

# --- milk proteome (MISAME-3), one representative ATE per protein (most significant) ---
milk <- readRDS(paste0(root, "results/adjusted_combined_arms_intervention_effects_proteomics_results_clean_ATE.RDS")) %>%
  filter(measure == "ATE", study == "Misame") %>%
  mutate(uniprot = up(biomarker)) %>%
  group_by(uniprot) %>% slice_min(pval, n = 1, with_ties = FALSE) %>% ungroup() %>%
  transmute(uniprot, milk_est = est, milk_pval = pval, milk_sigFDR = sigFDR, milk_visit = visit)

# --- maternal blood proteome, postnatal visits, representative ATE per protein ---
if (!exists("BLOOD_ADJUST")) BLOOD_ADJUST <- FALSE
.btag <- if (BLOOD_ADJUST) "adjusted_" else ""   # read adjusted_ blood results when set
.osuf <- if (BLOOD_ADJUST) "_adjusted" else ""   # output suffix (keep adjusted/unadjusted separate)
blood_all <- readRDS(paste0(root, "results/blood_compartment_", .btag,
  "combined_arms_intervention_effects_results_clean.RDS")) %>%
  filter(measure == "ATE", visit %in% c("pn12", "pn34", "pn56", "acco"))   # postnatal (vs postnatal milk)

# Selenoproteins to highlight (BEP raised selenium in Burkina Faso -> selenoproteins)
SELENO <- c("P22352" = "GPX3", "P49908" = "SELENOP")

analyze <- function(blood_label) {
  blood <- blood_all %>% filter(dataset == blood_label) %>%
    mutate(uniprot = up(biomarker)) %>%
    group_by(uniprot) %>% slice_min(pval, n = 1, with_ties = FALSE) %>% ungroup() %>%
    transmute(uniprot, blood_est = est, blood_pval = pval, blood_sigFDR = sigFDR, blood_visit = visit)

  ov <- inner_join(milk, blood, by = "uniprot") %>%
    mutate(concordant   = sign(milk_est) == sign(blood_est),
           both_sig     = milk_sigFDR == 1 & blood_sigFDR == 1,
           either_sig   = milk_sigFDR == 1 | blood_sigFDR == 1,
           gene         = SELENO[uniprot],
           sig_class    = case_when(both_sig ~ "FDR-sig in both",
                                    either_sig ~ "FDR-sig in one",
                                    TRUE ~ "ns"))

  n_shared <- nrow(ov)
  rho <- suppressWarnings(cor(ov$milk_est, ov$blood_est, method = "spearman"))
  conc_all  <- mean(ov$concordant)
  conc_eith <- with(filter(ov, either_sig), mean(concordant))
  conc_both <- with(filter(ov, both_sig),  mean(concordant))
  cat(sprintf("\n=== MILK vs %s ===\n", blood_label))
  cat(sprintf("shared proteins: %d | Spearman rho(ATE): %.2f\n", n_shared, rho))
  cat(sprintf("concordant direction: all %.0f%% | among either-sig %.0f%% (n=%d) | among both-sig %.0f%% (n=%d)\n",
              100*conc_all, 100*conc_eith, sum(ov$either_sig), 100*conc_both, sum(ov$both_sig)))
  cat("selenoproteins present:\n")
  print(ov %>% filter(!is.na(gene)) %>% select(uniprot, gene, milk_est, milk_sigFDR, blood_est, blood_sigFDR, concordant))

  p <- ggplot(ov, aes(milk_est, blood_est)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey70") +
    geom_point(aes(colour = sig_class), alpha = 0.6) +
    geom_text_repel(data = filter(ov, either_sig | !is.na(gene)),
                    aes(label = ifelse(!is.na(gene), gene, uniprot)), size = 3, max.overlaps = 20) +
    scale_colour_manual(values = c("FDR-sig in both" = "#E69F00", "FDR-sig in one" = "#56B4E9", "ns" = "grey75")) +
    labs(x = "Milk proteome ATE (BEP vs control)", y = paste0("Maternal blood ATE (", blood_label, ")"),
         colour = NULL,
         title = paste0("Cross-compartment proteome concordance: milk vs maternal blood (", blood_label, ")"),
         subtitle = sprintf("%d shared proteins; Spearman rho=%.2f; %.0f%% same-direction among either-FDR-sig",
                            n_shared, rho, 100*conc_eith)) +
    theme_bw() + theme(legend.position = "top")
  ggsave(p, file = paste0(outfig, "/proteome_milk_vs_blood_", tolower(sub("Proteomics", "", blood_label)), .osuf, ".png"),
         width = 8, height = 7)

  ov %>% mutate(blood_dataset = blood_label)
}

res <- bind_rows(analyze("ProteomicsDepleted"), analyze("ProteomicsNaive"))
write.csv(res, paste0(root, "results/cross_compartment_proteome_overlap", .osuf, ".csv"), row.names = FALSE)
saveRDS(res, paste0(root, "results/cross_compartment_proteome_summary", .osuf, ".RDS"))
cat("\nWrote results/cross_compartment_proteome_overlap.csv and figures/cross_compartment/\n")
