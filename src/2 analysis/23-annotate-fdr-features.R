# =============================================================================
# 23-annotate-fdr-features.R
#
# Putative annotation of the FDR-significant intervention-effect features, using
# the best available source per feature:
#   1. Sapient catalogue name  (V3 VAMS only; the ~444/38,761 named features in
#      metabolite_description_vam_with_global_id.csv) -> higher confidence
#   2. else Mummichog EmpiricalCompound candidates (m/z -> KEGG compounds via the
#      human_mfn network, MUM_PPM, adduct-aware) -> PUTATIVE, often AMBIGUOUS
#      (n_candidates column = how many distinct compounds share that mass).
#
# Mummichog is a PATHWAY tool, not an identifier: these per-feature candidates are
# MSI level ~3 (mass-only), frequently multiple per m/z, and include some
# malformed/partial names straight from the mummichog output. Confident IDs
# (MSI 1-2) still require Kim/Sapient MS/MS. The carnitine sibling vam_1005524 is
# an example that no source can place.
#
# Inputs : results/blood_compartment_adjusted_combined_arms_..._clean.RDS (FDR-sig ATEs)
#          results/mummichog_output_adjusted/<comp>_<visit>_<mode>/tables/userInput_to_EmpiricalCompounds.tsv
#          data/additional datasets/metabolite_description_vam_with_global_id.csv (Sapient names)
# Output : results/fdr_sig_putative_annotation.csv
# Scope  : combined-arms primary, MISAME-III, the compartments that have a mummichog
#          run (MaternalPlasma pn12; VamsPostnatalInfant pn12/pn34/pn56; VamsPostnatalMaternal pn56).
# =============================================================================

suppressMessages({library(data.table)})
root <- paste0(here::here(), "/"); source(paste0(root, "src/2 analysis/_blood_helpers.R"))
mcgroot <- paste0(root, "results/mummichog_output_adjusted/")

sap <- fread(paste0(root, "data/additional datasets/metabolite_description_vam_with_global_id.csv"))
sap <- sap[, .(feature = toupper(feature_label), sapient_name = ifelse(is.na(ID) | ID == "", NA_character_, ID))]

ca <- as.data.table(readRDS(paste0(root, "results/blood_compartment_adjusted_combined_arms_intervention_effects_results_clean.RDS")))[measure == "ATE" & contrast == "BEP" & sigFDR == 1]
vm <- mzrt_vams(); pl <- mzrt_rlc("ProcessedDataMISAME3_plasma.csv")

feat_tbl <- function(ds, vis, mz) {
  d <- ca[dataset == ds & visit == vis]; d[, feature := toupper(biomarker)]
  d <- d[order(pval)][!duplicated(feature)]
  merge(d[, .(dataset, visit, feature, est = round(est, 2), q = signif(pval_adj, 2))], mz[, .(feature, mz, mode)], by = "feature")
}
# read a mummichog run's per-feature candidate compounds (matched by directory suffix)
read_mcg <- function(ds, vis, mode) {
  alld <- list.files(mcgroot, full.names = TRUE)
  dir <- alld[endsWith(basename(alld), paste0(ds, "_", vis, "_", mode))]
  if (length(dir) == 0) return(NULL)
  f <- paste0(dir[1], "/tables/userInput_to_EmpiricalCompounds.tsv"); if (!file.exists(f)) return(NULL)
  e <- fread(f, sep = "\t", quote = ""); setnames(e, "m/z", "mz", skip_absent = TRUE)
  e[, .(mz = as.numeric(mz), compound_names)][!is.na(mz)]
}
# distinct candidate compounds (mummichog uses `$` between compounds, `;` between synonyms)
parse_names <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (!length(x)) return(character())
  compounds <- unlist(strsplit(paste(x, collapse = "$"), "$", fixed = TRUE))
  compounds <- trimws(sub(";.*", "", compounds))   # keep only the first synonym of each compound
  compounds <- compounds[compounds != ""]
  unique(compounds)
}
annotate <- function(fe) {
  out <- list()
  for (md in c("positive", "negative")) {
    sub <- fe[mode == md]; if (!nrow(sub)) next
    ec <- read_mcg(sub$dataset[1], sub$visit[1], md); if (is.null(ec)) next
    for (i in seq_len(nrow(sub))) {
      hit <- ec[abs(mz - sub$mz[i]) <= 0.003]; nm <- parse_names(hit$compound_names)
      out[[length(out) + 1]] <- data.table(feature = sub$feature[i], n_candidates = length(nm),
        putative_candidates = paste(head(nm, 5), collapse = " | "))
    }
  }
  rbindlist(out)
}

sets <- list(c("MaternalPlasma","pn12"), c("VamsPostnatalInfant","pn12"), c("VamsPostnatalInfant","pn34"),
             c("VamsPostnatalInfant","pn56"), c("VamsPostnatalMaternal","pn56"))
res <- rbindlist(lapply(sets, function(s) {
  mz <- if (grepl("Vams", s[1])) vm else pl
  fe <- feat_tbl(s[1], s[2], mz); a <- annotate(fe)
  m <- merge(fe, sap, by = "feature", all.x = TRUE); m <- merge(m, a, by = "feature", all.x = TRUE)
  m[is.na(n_candidates), n_candidates := 0L]; m[is.na(putative_candidates), putative_candidates := ""]
  m[order(q)]
}), fill = TRUE)
first_cand <- function(x) vapply(strsplit(x, " | ", fixed = TRUE), function(z) z[1], character(1))
res[, best_annotation := fifelse(!is.na(sapient_name), sapient_name,
       fifelse(n_candidates > 0, first_cand(putative_candidates), "- none -"))]
setcolorder(res, c("dataset","visit","feature","est","q","mode","sapient_name","best_annotation","n_candidates","putative_candidates"))
fwrite(res, paste0(root, "results/fdr_sig_putative_annotation.csv"))
cat("Saved results/fdr_sig_putative_annotation.csv (", nrow(res), "FDR-sig features )\n")
