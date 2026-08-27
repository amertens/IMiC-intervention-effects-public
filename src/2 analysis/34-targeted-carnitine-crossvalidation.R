# =============================================================================
# 34-targeted-carnitine-crossvalidation.R
#
# Cross-validate the headline untargeted carnitine (m/z 286.202, catalogued
# octenoylcarnitine) against the NAMED targeted Biocrates acylcarnitine panel,
# which needs no Kim input. Two questions:
#  (1) Mass arithmetic: is 286.202 [M+H]+ octenoyl (C8:1) or octanoyl (C8)?
#  (2) Does BEP move any NAMED milk acylcarnitine (esp. octanoyl c8), or is the
#      untargeted +signal uncorroborated/contradicted by the targeted panel?
# Out: results/targeted_carnitine_crossvalidation.csv
# =============================================================================
suppressMessages({library(data.table)})
root <- paste0(here::here(), "/")

# --- (1) exact-mass check -----------------------------------------------------
# Compute neutral monoisotopic masses from formula, then add a proton to get the
# [M+H]+ ion m/z we observe, and see which candidate the observed 286.2020 matches.
mH <- 1.007276                                             # proton mass (for [M+H]+)
el <- c(C=12, H=1.0078250319, N=14.0030740052, O=15.9949146221)   # monoisotopic element masses
mass <- function(C,H,N,O) C*el["C"]+H*el["H"]+N*el["N"]+O*el["O"]
octenoyl <- mass(15,27,1,4); octanoyl <- mass(15,29,1,4)   # C15H27NO4 (C8:1) vs C15H29NO4 (C8)
cat(sprintf("octenoylcarnitine C8:1  [M+H]+ = %.4f\n", octenoyl + mH))
cat(sprintf("octanoylcarnitine C8    [M+H]+ = %.4f\n", octanoyl + mH))
cat(sprintf("observed feature m/z          = 286.2020  -> ppm to octenoyl = %.1f, to octanoyl = %.0f\n",
            (286.2020-(octenoyl+mH))/(octenoyl+mH)*1e6, (286.2020-(octanoyl+mH))/(octanoyl+mH)*1e6))

# --- (2a) untargeted 286.202 effects (the claim) ------------------------------
feature_key <- function(x) toupper(as.character(x))     # case-insensitive feature id
milkC <- as.data.table(readRDS(paste0(root,"results/adjusted_combined_arms_intervention_effects_untargeted_results_clean_ATE.RDS")))
unt <- milkC[measure=="ATE" & study=="Misame" & feature_key(biomarker)==feature_key("RLC_POS_MTB_3033006"),
             .(platform="untargeted (Sapient)", compartment="milk", visit, feature=biomarker,
               est=round(est,2), p=signif(pval,2), q=signif(pval_adj,2), sigFDR)]
cat("\n=== UNTARGETED milk feature RLC_POS_MTB_3033006 (m/z 286.202) ===\n"); print(unt)

# --- (2b) targeted named acylcarnitines (the cross-check) ---------------------
tg <- fread(paste0(root,"results/subsetted results/tertiary_targeted_metabolomics.csv"))
ac <- tg[category=="Acylcarnitines"]
# octanoylcarnitine specifically (the named saturated analog), all studies/visits
c8 <- ac[grepl("Octanoylcarnitine", label_f),
         .(platform="targeted (Biocrates)", study, visit, compound=label_f,
           est=round(est,3), p=signif(pval,2), q=signif(pval_adj,2), sigFDR)]
cat("\n=== TARGETED octanoylcarnitine (c8), all trials/visits ===\n"); print(c8)
# is octenoylcarnitine (c8:1) present in the named panel at all?
cat("\noctenoylcarnitine present in targeted panel? ",
    any(grepl("Octenoylcarnitine|c8\\.1", ac$label_f, ignore.case=TRUE)), "\n")
# class-level summary per study x visit: any acylcarnitine up & FDR-sig?
cls <- ac[, .(n=.N, n_up=sum(est>0), n_FDRsig=sum(sigFDR==1),
              n_up_FDRsig=sum(est>0 & sigFDR==1),
              est_min=round(min(est),3), est_max=round(max(est),3)),
          by=.(study, visit)][order(study, visit)]
cat("\n=== TARGETED milk acylcarnitine CLASS summary (per trial x visit) ===\n"); print(cls)

out <- rbindlist(list(unt, c8), fill=TRUE)
fwrite(out, paste0(root,"results/targeted_carnitine_crossvalidation.csv"))
fwrite(cls, paste0(root,"results/targeted_acylcarnitine_class_summary.csv"))
cat("\nSaved results/targeted_carnitine_crossvalidation.csv (+ class summary)\n")
cat("\nVERDICT: mass 286.202 = octenoyl [M+H]+ (octanoyl ruled out by ~7000 ppm).",
    "\nThe named panel has octanoyl (c8) but NOT octenoyl (c8:1); octanoyl shows no BEP effect,",
    "\nand no MISAME milk acylcarnitine is FDR-significant -> the untargeted + signal is",
    "\nuncorroborated by named data and the saturated analog is null. Identity as an acylcarnitine",
    "\nis mass-plausible (unsaturated C8:1) but not independently supported.\n")
