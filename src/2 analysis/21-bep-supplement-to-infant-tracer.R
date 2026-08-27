# =============================================================================
# 21-bep-supplement-to-infant-tracer.R   [cross-platform version]
#
# CORRECTION (2026-06-25): the BEP-supplement metabolomics are on the **V1 (LCV1)**
# rLC catalogue (recovered from prenatalVamsImputedCappedScaled.csv -> 27 supplement
# samples; profile in results/bep_supplement_profile_V1.csv). They are the SAME
# method as milk / maternal plasma / prenatal VAMS, but feature IDs (`RLC_*_MTB_*`)
# are assigned PER DATASET, so even within V1 there is no shared-ID join -> matching
# is by accurate mass (25 ppm + ionization mode, ppm-only first pass, like script 14).
#   - supplement -> milk   : V1<->V1  (linear RT; the cleanest, highest-confidence leg)
#   - supplement -> plasma : V1<->V1
#   - supplement -> infant : V1<->V3  (postnatal VAMS; non-linear RT, ISOBARIC-limited
#                                      -> reported with that caveat, NOT a direct-ID tracer)
#
# QUESTION: are the metabolites ABUNDANT IN THE BEP PRODUCT (high supplement
# enrichment) preferentially ELEVATED (BEP-up) in milk / maternal plasma / infant
# blood? i.e. does the supplement's own metabolite signature show up as a positive
# intervention effect downstream -- a compositional "tracer" of supplement intake.
#
# Measurement constraint (Kim): every downstream effect is a WITHIN-dataset z-scored
# BEP-vs-control ATE; we compare DIRECTION only, never absolute intensity across
# datasets. supp_mean = enrichment of a feature in the supplement vs maternal VAMS.
#
# Out: results/bep_supplement_tracer.csv  (per target: matched n, %up, enrichment p)
# =============================================================================

suppressMessages({library(data.table)})
root <- paste0(here::here(), "/")
source(paste0(root, "src/2 analysis/_blood_helpers.R"))

SUPP_HI <- 1.0   # "abundant in supplement" = supp_mean >= 1 SD above the maternal-VAMS mean

# -- supplement profile + its m/z-RT (from the prenatal-VAMS V1 catalogue) --------
supp    <- fread(paste0(root, "results/bep_supplement_profile_V1.csv"))
supp_hi <- toupper(supp[supp_mean >= SUPP_HI]$feature)
pv_mz   <- mzrt_rlc("ProcessedDataMISAME3_VAMS.csv")          # supplement & prenatal VAMS share this V1 catalogue
supp_mz <- pv_mz[feature %in% supp_hi]
cat(sprintf("Supplement-abundant features (supp_mean >= %.1f): %d  (with m/z-RT: %d)\n",
            SUPP_HI, length(supp_hi), nrow(supp_mz)))

# -- downstream BEP ATEs (combined-arms primary, adjusted) -----------------------
rep_feat <- function(d) { d <- copy(d); d[, feature := toupper(biomarker)]; d[order(pval)][!duplicated(feature)] }
milk <- rep_feat(as.data.table(readRDS(paste0(root,"results/adjusted_combined_arms_intervention_effects_untargeted_results_clean_ATE.RDS")))[
          measure=="ATE" & study=="Misame" & contrast=="BEP" & visit %in% c("1-2 mo.","3-4 mo.")])
blood<- as.data.table(readRDS(paste0(root,"results/blood_compartment_adjusted_combined_arms_intervention_effects_results_clean.RDS")))[measure=="ATE" & contrast=="BEP"]
plas <- rep_feat(blood[dataset=="MaternalPlasma" & visit=="pn12"])
inf  <- rep_feat(blood[dataset=="VamsPostnatalInfant" & visit %in% c("pn12","pn34","pn56")])

mk_mz   <- mzrt_milk(misame_only=FALSE)
plas_mz <- mzrt_rlc("ProcessedDataMISAME3_plasma.csv")
vm_mz   <- mzrt_vams()

# enrichment of supplement-matched features among BEP-UP features of a target
tracer <- function(tgt, tgt_mz, label, v1v3) {
  rt <- if (v1v3) CC_RT_V1V3 else CC_RT_V1V1
  pr <- match_mzrt(supp_mz, tgt_mz[feature %in% tgt$feature], rt)   # 25 ppm + mode, ppm-only
  matched <- unique(pr$feature_B)
  tt <- tgt[feature %in% matched]
  bg_up <- mean(tgt$est > 0); m_up <- if (nrow(tt)) mean(tt$est > 0) else NA_real_
  tgt2 <- copy(tgt)[, `:=`(mtch = feature %in% matched, up = est > 0)]
  ft <- tryCatch(fisher.test(table(factor(tgt2$mtch,c(FALSE,TRUE)), factor(tgt2$up,c(FALSE,TRUE))), alternative="greater"),
                 error=function(e) list(p.value=NA_real_, estimate=NA_real_))
  hits <- tt[sigFDR==1][order(-est)]
  cat(sprintf("\n== supplement -> %-15s (%s) ==\n", label, if(v1v3)"V1<->V3 isobaric-limited" else "V1<->V1 clean"))
  cat(sprintf("  target features: %d | matched to a supplement-abundant feature (25 ppm): %d\n", nrow(tgt), nrow(tt)))
  cat(sprintf("  %% BEP-up among matched: %s%%  vs background %% up: %d%%  (Fisher enriched-among-up p=%s, OR=%s)\n",
              ifelse(is.na(m_up),"NA",round(100*m_up)), round(100*bg_up),
              ifelse(is.na(ft$p.value),"NA",signif(ft$p.value,2)),
              ifelse(is.null(ft$estimate)||is.na(ft$estimate),"NA",round(as.numeric(ft$estimate),2))))
  cat(sprintf("  matched & FDR-sig downstream: %d\n", nrow(hits)))
  data.table(target=label, platform=if(v1v3)"V1<->V3 (isobaric)" else "V1<->V1 (clean)",
             n_target=nrow(tgt), n_matched=nrow(tt),
             pct_up_matched=ifelse(is.na(m_up),NA,round(100*m_up)), pct_up_background=round(100*bg_up),
             fisher_up_p=ifelse(is.na(ft$p.value),NA,signif(ft$p.value,3)),
             n_matched_FDRsig=nrow(hits))
}

res <- rbindlist(list(
  tracer(milk, mk_mz,   "milk",            v1v3=FALSE),
  tracer(plas, plas_mz, "maternal plasma", v1v3=FALSE),
  tracer(inf,  vm_mz,   "infant blood",    v1v3=TRUE)
))
fwrite(res, paste0(root, "results/bep_supplement_tracer.csv"))
cat("\nSaved: results/bep_supplement_tracer.csv\n")
cat("NOTE: supplement->milk and ->plasma are V1<->V1 (clean mass-match, ~linear RT).\n")
cat("supplement->infant is V1<->V3 (isobaric-limited) -> directional hint only, needs Kim annotation.\n")
