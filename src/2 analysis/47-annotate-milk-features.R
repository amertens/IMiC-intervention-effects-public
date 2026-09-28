# =============================================================================
# 47-annotate-milk-features.R
#
# Putative names for the untargeted milk metabolite features with a BEP effect
# (covariate-adjusted combined arms, all studies), in two foregrounds: nominal
# P < 0.05, the enrichment query behind Fig. 6A and Table S5
# (src/metaboanalyst/run-untargeted-msea.R), and FDR-significant, used by script 54
# for supplement matching and cross-compartment linkage. Mirrors the blood version
# (23-annotate-fdr-features.R). Best available source per feature:
#   1. the reconciled `name_final` from
#      data/untargeted_annotation/untargeted_to_annotate_manually_annotated.xlsx;
#   2. otherwise the curated compound-ID name (`ID` column of the per-study
#      compound_ID tables, e.g. "(Iso)Butyrylcarnitine_[M+H]+");
#   3. otherwise Mummichog EmpiricalCompound candidates from the directional milk
#      runs of script 18, matched within 0.003 m/z; putative and often ambiguous
#      (n_candidates = number of distinct compounds sharing that mass).
# Mummichog candidates are MSI level ~3 (mass only); confident identities need
# tandem MS. Unresolved features are labelled "- none -". Reads existing Mummichog
# output only; it does not run Mummichog.
#
# Inputs : results/adjusted_combined_arms_intervention_effects_untargeted_results_clean_ATE.RDS
#          data/untargeted_annotation/untargeted_to_annotate_manually_annotated.xlsx
#          data/additional datasets/MISAME3_metabolite_description_compound_ID_combined_20250731.csv
#          data/additional datasets/CHILD_ELICIT_VITAL_metabolite_description_compound_ID_FINAL_20250730.csv
#          results/mummichog_output_directional/<run>.Milk_{positive,negative}_{up,down}/tables/
#            userInput_to_EmpiricalCompounds.tsv (script 18)
# Outputs: results/milk_nominal_putative_annotation.csv  (nominal P < 0.05 foreground)
#          results/milk_fdr_sig_putative_annotation.csv  (FDR-significant foreground)
# [needs restricted data]
# =============================================================================

suppressMessages({library(data.table); library(readxl)})
root <- paste0(here::here(), "/")
mcgroot <- paste0(root, "results/mummichog_output_directional/")

# --- primary feature annotation: the reconciled `name_final` ------------------------
# untargeted_to_annotate_manually_annotated.xlsx combines the identifier pipelines
# (MS2-verified, "global" library IDs, MetID, and the annotation provider, Sapient) into
# one reconciled `name_final` per feature. It is the annotation key the Fig. 6A MSEA
# foreground is named with (Methods: metabolites "mapped to matched metabolite names
# using the corresponding annotation keys"), and it comes from instrument/software
# pipelines, not the MetaboAnalyst web tool. Only about 506 of 36,853 features carry a
# confident name_final (untargeted annotation is sparse), so the per-study compound-ID
# `ID` is the next tier and local Mummichog candidates the last resort.
name_final_lookup <- {
  nf <- as.data.table(readxl::read_excel(
    paste0(root, "data/untargeted_annotation/untargeted_to_annotate_manually_annotated.xlsx"),
    sheet = 1, col_types = "text"))
  nf[, .(feature = toupper(Metabolite_Feature_Label),
         name_final = ifelse(is.na(name_final) | trimws(name_final) == "", NA_character_, trimws(name_final)))
     ][!is.na(name_final)][!duplicated(feature)]
}

# --- per-study curated compound-ID tables (feature -> curated name + m/z) -------
# ID is the curated compound name where present (blank/NA when unknown).
load_cid <- function(file) {
  d <- fread(paste0(root, "data/additional datasets/", file),
             select = c("Metabolite_Feature_Label", "MZ", "ID"))
  d[, .(feature = toupper(Metabolite_Feature_Label), mz = as.numeric(MZ),
        curated_name = ifelse(is.na(ID) | trimws(ID) == "", NA_character_, trimws(ID)))]
}
cid_mis <- load_cid("MISAME3_metabolite_description_compound_ID_combined_20250731.csv")
cid_cev <- load_cid("CHILD_ELICIT_VITAL_metabolite_description_compound_ID_FINAL_20250730.csv")
# de-duplicate to one curated row per feature (prefer a non-NA curated name)
dedup_cid <- function(x) x[order(is.na(curated_name))][!duplicated(feature)]
cid_mis <- dedup_cid(cid_mis); cid_cev <- dedup_cid(cid_cev)

# --- shared naming helpers (local only) -----------------------------------------
# Milk untargeted combined-arms ATE. Feature = biomarker (rLC_pos_mtb_* / rLC_neg_mtb_*).
ca_all <- as.data.table(readRDS(paste0(root,
  "results/adjusted_combined_arms_intervention_effects_untargeted_results_clean_ATE.RDS")))

# mummichog fallback: read one directional milk EmpiricalCompound table (local output)
read_mcg <- function(mode, dir) {
  alld <- list.files(mcgroot, full.names = TRUE)
  d <- alld[endsWith(basename(alld), paste0("Milk_", mode, "_", dir))]
  if (length(d) == 0) return(NULL)
  f <- paste0(d[1], "/tables/userInput_to_EmpiricalCompounds.tsv")
  if (!file.exists(f)) return(NULL)
  e <- fread(f, sep = "\t", quote = ""); setnames(e, "m/z", "mz", skip_absent = TRUE)
  e[, .(mz = as.numeric(mz), compound_names)][!is.na(mz)]
}
# mummichog: `$` between compounds, `;` between synonyms -> first synonym of each compound
parse_names <- function(x) {
  x <- x[!is.na(x) & x != ""]; if (!length(x)) return(character())
  compounds <- unlist(strsplit(paste(x, collapse = "$"), "$", fixed = TRUE))
  compounds <- trimws(sub(";.*", "", compounds)); compounds <- compounds[compounds != ""]
  unique(compounds)
}
first_cand <- function(x) vapply(strsplit(x, " | ", fixed = TRUE), function(z) z[1], character(1))
mcg <- list()
for (md in c("positive", "negative")) for (dr in c("up", "down"))
  mcg[[paste(md, dr)]] <- read_mcg(md, dr)

# --- build one annotation table: features with `pcol` < 0.05, named locally ------
# pcol = "pval_adj" -> FDR-significant foreground; "pval" -> nominal P<0.05 foreground
# (the Fig. 6A over-representation foreground: "unadjusted P < 0.05").
build_annotation <- function(pcol, outfile, label) {
  ca <- ca_all[measure == "ATE" & !is.na(get(pcol)) & get(pcol) < 0.05 & !is.na(est)]
  ca[, feature := toupper(biomarker)]
  ca[, mode := fifelse(grepl("_POS_", feature), "positive",
               fifelse(grepl("_NEG_", feature), "negative", NA_character_))]
  ca <- ca[!is.na(mode)]                                     # drop non-metabolite rows (e.g. "Arm")
  ca <- ca[order(pval)][!duplicated(paste(study, feature))]  # best occurrence per study x feature
  ca[, direction := fifelse(est > 0, "up", "down")]

  fe <- ca[, .(study, studytime, visit, contrast, feature,
               est = round(est, 2), q = signif(pval_adj, 2), mode, direction)]
  fe[, `:=`(curated_name = NA_character_, mz = NA_real_)]
  fe[study == "Misame", c("curated_name", "mz") :=
       cid_mis[.SD, .(curated_name, mz), on = "feature"], .SDcols = "feature"]
  fe[study != "Misame", c("curated_name", "mz") :=
       cid_cev[.SD, .(curated_name, mz), on = "feature"], .SDcols = "feature"]
  # name_final (the reconciled MS2/global/MetID/Sapient identity) overrides the per-study
  # lab ID where both exist: it is the paper's annotation key, and on the 33 features
  # where the two disagree, name_final is the one the paper used. Features name_final
  # does not cover keep the lab ID (then the Mummichog fallback).
  fe[name_final_lookup, name_final := i.name_final, on = "feature"]
  fe[!is.na(name_final), curated_name := name_final]
  fe[, name_final := NULL]

  fe[, `:=`(n_candidates = 0L, putative_candidates = "")]
  for (i in seq_len(nrow(fe))) {
    ec <- mcg[[paste(fe$mode[i], fe$direction[i])]]
    if (is.null(ec) || is.na(fe$mz[i])) next
    hit <- ec[abs(mz - fe$mz[i]) <= 0.003]
    cand_names <- parse_names(hit$compound_names)
    set(fe, i, "n_candidates", length(cand_names))
    set(fe, i, "putative_candidates", paste(head(cand_names, 5), collapse = " | "))
  }
  fe[, best_annotation := fifelse(!is.na(curated_name), curated_name,
        fifelse(n_candidates > 0, first_cand(putative_candidates), "- none -"))]
  setcolorder(fe, c("study", "studytime", "visit", "contrast", "feature", "est", "q",
                    "mode", "direction", "curated_name", "n_candidates",
                    "putative_candidates", "best_annotation"))
  fe <- fe[order(q)]
  fwrite(fe, paste0(root, outfile))
  cat(sprintf("Saved %s (%d %s features; curated=%d, mummichog-only=%d, none=%d)\n",
      outfile, nrow(fe), label, sum(!is.na(fe$curated_name)),
      sum(is.na(fe$curated_name) & fe$n_candidates > 0),
      sum(is.na(fe$curated_name) & fe$n_candidates == 0)))
  invisible(fe)
}

# FDR-significant foreground -> script 54 (supplement matching, Fig. 6D linkage).
build_annotation("pval_adj", "results/milk_fdr_sig_putative_annotation.csv", "FDR-sig")
# Nominal P < 0.05 foreground -> the Fig. 6A / Table S5 enrichment query.
build_annotation("pval",     "results/milk_nominal_putative_annotation.csv", "nominal P<0.05")
