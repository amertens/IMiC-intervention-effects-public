# =============================================================================
# 43-acylcarnitine-composition.R
#
# April's suggestion (2026-06-30 IMiC meeting): look at the COMPOSITION of the
# targeted (Biocrates) acylcarnitines and free carnitine -- the ratios of
# acylcarnitines to free carnitine and of acylcarnitines to each other -- to tell
# apart two readings of the untargeted carnitine signal:
#   (a) "more carnitine substrate available"  -> acylcarnitines and free carnitine
#       move together, the ACYLCARNITINE-TO-FREE-CARNITINE RATIO is roughly flat; vs
#   (b) "a bottleneck in lipid processing"     -> acylcarnitines accumulate relative
#       to free carnitine, the RATIO rises.
#
# Uses the per-sample milk Biocrates codes in merged_analysis_datasets.RDS (c0 = free
# carnitine; c2..c18 = acylcarnitines, by acyl-chain carbon number). MISAME milk only
# (the BEP trial). Exposure = postnatal BEP (received BEP after birth: arms IFA/BEP +
# BEP/BEP), matching the postnatal, milk-borne carnitine signal. Per milk visit
# (1 = 14-21 d, 2 = 1-2 mo, 3 = 3-4 mo). Exploratory unadjusted contrast on the
# log-ratio; FDR across the ratio set per visit.
#
# Out: results/acylcarnitine_composition.csv
# =============================================================================
suppressMessages({library(data.table)})
root <- paste0(here::here(), "/")

milk <- as.data.table(readRDS(paste0(root, "data/merged_analysis_datasets.RDS")))[study == "Misame"]
codes <- grep("^c[0-9]", names(milk), value = TRUE, ignore.case = TRUE)
carbon <- as.integer(sub("^c([0-9]+).*", "\\1", codes))             # acyl-chain carbon number
free   <- "c0"
ac_all <- codes[codes != free]
scac <- codes[carbon %in% 2:5]                                       # short-chain (C2-C5)
mcac <- codes[carbon %in% 6:13]                                      # medium-chain (C6-C12, incl octanoyl c8)
lcac <- codes[carbon >= 14]                                          # long-chain (C14-C18)

vlab <- c(`1` = "14-21 d", `2` = "1-2 mo", `3` = "3-4 mo")
milk[, bep_post := as.integer(arm %in% c("IFA/BEP", "BEP/BEP"))]     # received BEP after birth

# per-sample pooled concentrations by chain-length class (sum the member Biocrates codes)
row_total <- function(cols) rowSums(as.matrix(milk[, ..cols]), na.rm = TRUE)
milk[, `:=`(C0   = get(free),         # free carnitine
            AC   = row_total(ac_all), # all acylcarnitines
            SCAC = row_total(scac),   # short-chain sum
            MCAC = row_total(mcac),   # medium-chain sum
            LCAC = row_total(lcac),   # long-chain sum
            C8   = if ("c8" %in% codes) get("c8") else NA_real_)]  # octanoylcarnitine (named C8 analog)

# ratios April named: acylcarnitine-to-free-carnitine (overall + by chain length),
# the specific medium-chain octanoyl, and acylcarnitine-to-acylcarnitine balances.
ratios <- list(
  `total acylcarnitine / free carnitine` = quote(AC / C0),
  `short-chain / free carnitine`         = quote(SCAC / C0),
  `medium-chain / free carnitine`        = quote(MCAC / C0),
  `long-chain / free carnitine`          = quote(LCAC / C0),
  `octanoylcarnitine (C8) / free carnitine` = quote(C8 / C0),
  `medium-chain / short-chain`           = quote(MCAC / SCAC),
  `long-chain / short-chain`             = quote(LCAC / SCAC),
  `acylcarnitine fraction AC/(AC+C0)`    = quote(AC / (AC + C0)))

# BEP (postnatal) effect on each log-ratio, per visit.
# The linear model is on log(ratio), so the bep_post coefficient is a log-fold-change
# (exp - 1 gives the % difference vs control). Require >=20 samples and both arms present.
out <- rbindlist(lapply(names(ratios), function(rn) {
  milk[, r := eval(ratios[[rn]])]                          # evaluate this ratio for every sample
  rbindlist(lapply(sort(unique(milk$visit)), function(vv) {
    d <- milk[visit == vv & is.finite(r) & r > 0]
    if (nrow(d) < 20 || uniqueN(d$bep_post) < 2) return(NULL)
    d[, lr := log(r)]
    fit <- summary(lm(lr ~ bep_post, d))$coefficients
    log_fc <- fit["bep_post", "Estimate"]
    data.table(ratio = rn, visit = vlab[as.character(vv)], n = nrow(d),
               mean_ctrl = round(mean(d[bep_post == 0]$r), 3),
               mean_bep  = round(mean(d[bep_post == 1]$r), 3),
               log_fold_change = round(log_fc, 3),
               pct_change = round(100 * (exp(log_fc) - 1), 1),
               p = signif(fit["bep_post", "Pr(>|t|)"], 2))
  }))
}))
out[, p_adj := p.adjust(p, "BH"), by = visit]                        # FDR across the ratio set, per visit
out[, sigFDR := as.integer(p_adj < 0.05)]
fwrite(out, paste0(root, "results/acylcarnitine_composition.csv"))

cat("=== Targeted milk acylcarnitine composition: postnatal-BEP effect on ratios ===\n")
cat(sprintf("MISAME milk; free carnitine=c0; %d acylcarnitines (SC %d / MC %d / LC %d)\n",
            length(ac_all), length(scac), length(mcac), length(lcac)))
print(out[order(ratio, visit)])
cat("\nRead: a POSITIVE log-fold-change in 'total acylcarnitine / free carnitine' = acylcarnitines\n",
    "accumulate relative to free carnitine (processing-bottleneck reading); ~zero = the pool simply\n",
    "scales (substrate-availability reading). The C8/C0 row is the named analog of the untargeted\n",
    "octenoylcarnitine (m/z 286.202). FDR is across the ratio set within each visit.\n")
cat("\nSaved: results/acylcarnitine_composition.csv\n")
