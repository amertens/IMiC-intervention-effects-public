# =============================================================================
# 19b-cross-compartment-fdr-first-lists.R
#
# Trenton's FDR-FIRST reorientation (Andrew x Trenton, last session): instead of
# the bulk p<0.05 %-concordance (script 19) we want the small, reviewer-proof
# lists --- the features that are FDR-significant on BOTH sides of an arrow,
# matched within 25 ppm (cross-platform) or by shared id (postnatal VAMS), each
# annotated individually so they can be interpreted / sent to Kim.
#
# Two layers:
#   (A) COMBINED-ARMS PRIMARY (any-postnatal-BEP vs control) --- read the
#       both_sig pairs already written by script 14 (ppm-only) and annotate.
#   (B) GROUP-STRATIFIED (BEP/BEP, IFA/BEP, BEP/IFA vs control) --- compute the
#       FDR-sig-in-both matched pairs per contrast for the identity-anchored
#       maternal<->infant VAMS arrow + the cross-platform milk/plasma arrows.
#
# Annotation is metadata-only: the postnatal-VAMS catalogue carries names for
# ~444/38,761 features; milk + plasma (V1 rLC) features are UNannotated in-repo
# (-> "-"), which is exactly what §9 asks Kim to fill. No numbers are invented.
#
# Out: results/cross_compartment_fdr_first_lists.csv   (every both-FDR-sig pair, annotated)
#      console: the per-pair / per-contrast small lists for the memo
# =============================================================================

suppressMessages({library(data.table)})
root <- paste0(here::here(), "/")
source(paste0(root, "src/2 analysis/_blood_helpers.R"))

# -- annotation lookup (postnatal VAMS only; everything else stays unannotated) --
.vam_ann <- {
  d <- fread(paste0(root, "data/additional datasets/metabolite_description_vam_with_global_id.csv"))
  data.table(feature = toupper(d$feature_label),
             mz = d$mz, rt = d$rt_minute,
             annotation = ifelse(is.na(d$ID) | d$ID == "", NA_character_, d$ID),
             cclass = d[["Compound Class"]])
}
annotate <- function(features) {
  out <- .vam_ann[data.table(feature = toupper(features)), on = "feature",
                  .(feature, mz, rt, annotation, cclass)]
  out[is.na(annotation), annotation := "-"]
  out
}

# =============================================================================
# (A) COMBINED-ARMS PRIMARY --- both-FDR-sig matched features (from script 14)
# =============================================================================
pair_files <- c(
  "Milk<->MaternalPlasma" = "cross_compartment_metab_MilkVsMaternalPlasma_adjusted.csv",
  "Milk<->InfantVAMS"     = "cross_compartment_metab_MilkVsInfantVAMS_adjusted.csv",
  "MaternalVAMS<->InfantVAMS" = "cross_compartment_metab_MaternalVAMSvsInfantVAMS_adjusted.csv",
  "MaternalPlasma<->InfantVAMS"= "cross_compartment_metab_MaternalPlasmaVsInfantVAMS_adjusted.csv")

primary <- list()
cat("================ (A) COMBINED-ARMS PRIMARY (any-postnatal-BEP) ================\n")
for (nm in names(pair_files)) {
  f <- paste0(root, "results/", pair_files[nm])
  if (!file.exists(f)) { cat("  [missing]", pair_files[nm], "\n"); next }
  d <- fread(f)
  if (!"feature_A" %in% names(d) && "feature" %in% names(d))   # id-matched pair: one shared `feature` column
    d[, `:=`(feature_A = feature, feature_B = feature)]
  bs <- d[both_sig == TRUE]
  # de-duplicate to unique A-B feature pairs (a feature can match >1 partner)
  bs <- unique(bs, by = c("feature_A", "feature_B"))
  conc <- if (nrow(bs)) round(100*mean(bs$concordant)) else NA_integer_
  cat(sprintf("\n%-28s  both-FDR-sig pairs: %d   concordant: %s%%   (up-up %d / dn-dn %d / discordant %d)\n",
              nm, nrow(bs), ifelse(is.na(conc),"NA",conc),
              sum(bs$est_A>0 & bs$est_B>0), sum(bs$est_A<0 & bs$est_B<0), sum(!bs$concordant)))
  if (nrow(bs)) {
    ann <- annotate(bs$feature_B)                          # B is the postnatal-VAMS side where applicable
    bs[, `:=`(mz_B = ann$mz, annotation_B = ann$annotation)]
    print(bs[order(-abs(est_B))][, .(feature_A, est_A = round(est_A,2),
            feature_B, est_B = round(est_B,2), concordant, mz_B, annotation_B)], nrow = 60)
    primary[[nm]] <- cbind(layer = "combined-arms", arrow = nm, contrast = "any-postnatal-BEP", bs)
  }
}

# =============================================================================
# (B) GROUP-STRATIFIED --- FDR-sig-in-both per contrast
# =============================================================================
milk  <- as.data.table(readRDS(paste0(root,"results/adjusted_intervention_effects_res_untargeted_metabolomics_clean.RDS")))[measure=="ATE" & study=="Misame"]
blood <- as.data.table(readRDS(paste0(root,"results/blood_compartment_adjusted_intervention_effects_results_clean.RDS")))[measure=="ATE"]
# Feature-level tables, one representative (smallest-pval) row per feature.
get_milk  <- function(cc,vv){
  d <- milk[contrast==cc & visit==vv]
  d[, feature := toupper(biomarker)]
  d <- d[order(pval)][!duplicated(feature)]
  d[, .(feature, est, sigFDR)]
}
get_blood <- function(dd,vv,cc){
  d <- blood[dataset==dd & visit %in% vv & contrast==cc]
  d[, feature := toupper(biomarker)]
  d <- d[order(pval)][!duplicated(feature)]
  d[, .(feature, est, sigFDR)]
}
mk_mz<-mzrt_milk(misame_only=FALSE); pl_mz<-mzrt_rlc("ProcessedDataMISAME3_plasma.csv"); vm_mz<-mzrt_vams()

# matched FDR-sig-in-both pairs for one arrow + contrast
fdr_pairs <- function(A, B, mtype, mzA, mzB, cc) {
  selA <- (if (A[[1]]=="milk") get_milk(cc,A[[3]]) else get_blood(A[[2]],A[[3]],cc))[sigFDR==1]
  selB <- (if (B[[1]]=="milk") get_milk(cc,B[[3]]) else get_blood(B[[2]],B[[3]],cc))[sigFDR==1]
  if (!nrow(selA) || !nrow(selB)) return(data.table())
  if (mtype=="id") {
    m <- merge(selA, selB, by="feature", suffixes=c("_A","_B"))
    if (!nrow(m)) return(data.table())
    return(m[, .(feature_A=feature, est_A, feature_B=feature, est_B, concordant=sign(est_A)==sign(est_B))])
  }
  pr <- match_mzrt(mzA[feature %in% selA$feature], mzB[feature %in% selB$feature])  # 25 ppm + mode, ppm-only
  if (!nrow(pr)) return(data.table())
  pr <- merge(pr, selA[,.(feature_A=feature, est_A=est)], by="feature_A")
  pr <- merge(pr, selB[,.(feature_B=feature, est_B=est)], by="feature_B")
  unique(pr[, .(feature_A, est_A, feature_B, est_B, concordant=sign(est_A)==sign(est_B))])
}

strat_arrows <- list(
  list(name="MatVAMS pn56 <-> Infant pn56", A=list("blood","VamsPostnatalMaternal","pn56",vm_mz), B=list("blood","VamsPostnatalInfant","pn56",vm_mz), match="id"),
  list(name="Milk 1-2mo -> Infant pn12",    A=list("milk",NA,"1-2 mo.",mk_mz),                   B=list("blood","VamsPostnatalInfant","pn12",vm_mz), match="v1v3"),
  list(name="Milk 3-4mo -> Infant pn34",    A=list("milk",NA,"3-4 mo.",mk_mz),                   B=list("blood","VamsPostnatalInfant","pn34",vm_mz), match="v1v3"),
  list(name="Plasma pn12 -> Milk 1-2mo",    A=list("blood","MaternalPlasma","pn12",pl_mz),        B=list("milk",NA,"1-2 mo.",mk_mz),                  match="v1v1"))
CONTRASTS <- c("BEP/BEP","IFA/BEP","BEP/IFA")

strat <- list()
cat("\n================ (B) GROUP-STRATIFIED (FDR-sig in both) ================\n")
for (ar in strat_arrows) for (cc in CONTRASTS) {
  pr <- fdr_pairs(ar$A, ar$B, ar$match, ar$A[[4]], ar$B[[4]], cc)
  if (!nrow(pr)) next
  ann <- annotate(pr$feature_B)
  pr[, `:=`(mz_B = ann$mz, annotation_B = ann$annotation)]
  cat(sprintf("\n%-28s  %-8s  pairs: %d  concordant: %d%%\n", ar$name, cc, nrow(pr), round(100*mean(pr$concordant))))
  print(pr[order(-abs(est_B))][, .(feature_A, est_A=round(est_A,2), feature_B, est_B=round(est_B,2), concordant, mz_B, annotation_B)])
  strat[[length(strat)+1]] <- cbind(layer="stratified", arrow=ar$name, contrast=cc, pr)
}

# -- combined tidy output ------------------------------------------------------
allp <- rbindlist(c(
  lapply(primary, function(x) x[, .(layer, arrow, contrast, feature_A, est_A, feature_B, est_B, concordant)]),
  lapply(strat,   function(x) x[, .(layer, arrow, contrast, feature_A, est_A, feature_B, est_B, concordant)])),
  use.names = TRUE, fill = TRUE)
# carry the annotations through to the saved table (postnatal-VAMS names; "-" else)
vam_name <- function(f) {                                   # look up the VAMS-catalogue name, "-" if none
  a <- .vam_ann$annotation[match(toupper(f), .vam_ann$feature)]
  a[is.na(a)] <- "-"
  a
}
vam_mz <- function(f) .vam_ann$mz[match(toupper(f), .vam_ann$feature)]
allp[, `:=`(name_A = vam_name(feature_A), mz_A = vam_mz(feature_A),
            name_B = vam_name(feature_B), mz_B = vam_mz(feature_B))]
setcolorder(allp, c("layer","arrow","contrast","feature_A","name_A","est_A","mz_A",
                    "feature_B","name_B","est_B","mz_B","concordant"))
fwrite(allp, paste0(root, "results/cross_compartment_fdr_first_lists.csv"))
cat("\nSaved: results/cross_compartment_fdr_first_lists.csv  (", nrow(allp), "both-FDR-sig matched pairs;",
    sum(allp$name_A != "-" | allp$name_B != "-"), "with a named feature )\n")
