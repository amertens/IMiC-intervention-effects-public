# =============================================================================
# 32-signed-pathway-direction.R
#
# Single-null signed confirmation of the §5 directional Mummichog arrows. The §5
# table assigns up/down by running Mummichog twice (up- and down-subsets, with the
# other direction's p set to 1), which gives anti-conservative directional p-values
# because direction correlates with feature class. Here we instead take each
# pathway's member features (the features Mummichog matched into that pathway in the
# standard, non-directional run), look up the SIGN of their BEP effect from our ATE
# results, and test net direction with one binomial sign test (a single null). The
# direction call should reproduce §5; the p-value is now honest.
# Blood compartments (the new data): maternal plasma pn12, postnatal maternal VAMS
# pn56, infant VAMS pn12/34/56. Combined-arms (pooled postnatal-BEP) effects.
# Out: results/signed_pathway_direction.csv
# =============================================================================
suppressMessages({library(data.table)})
root <- paste0(here::here(), "/")
source(paste0(root, "src/2 analysis/_blood_helpers.R"))
up <- function(x) toupper(as.character(x))

bloodC <- as.data.table(readRDS(paste0(root,"results/blood_compartment_adjusted_combined_arms_intervention_effects_results_clean.RDS")))
plmz <- mzrt_rlc("ProcessedDataMISAME3_plasma.csv"); vmz <- mzrt_vams()

# Our per-feature signed statistic for a compartment, keyed by m/z so it can be
# joined to Mummichog's features. Sign = BEP direction, magnitude = -log10(p).
# pmax(pval, 1e-300) floors p away from 0 so log10 stays finite.
sig_tab <- function(ds, vv, mzt){
  d <- bloodC[measure=="ATE" & dataset==ds & visit==vv]; d[, feature := up(biomarker)]
  d <- merge(d, mzt[, .(feature, mz)], by="feature")
  d[, signed := sign(est) * -log10(pmax(pval, 1e-300))]
  d[is.finite(signed) & is.finite(mz), .(mz=round(mz,4), signed)][!duplicated(mz)]
}

# Parse one Mummichog run folder into a pathway -> member-feature-m/z table.
# Mummichog reports pathways in terms of "empirical compound" (EC) ids, so we
# bridge EC id -> m/z using its userInput table, then attach it to each pathway.
parse_run <- function(dir){
  tables_dir <- paste0(dir, "/tables")
  ui_files <- list.files(tables_dir, pattern="userInput_to_EmpiricalCompounds", full.names=TRUE)
  pa_files <- list.files(tables_dir, pattern="mcg_pathwayanalysis", full.names=TRUE)
  if(!length(ui_files) || !length(pa_files)) return(NULL)

  # EC id -> feature m/z map (the "m/z" column name needs get() because of the slash)
  user_input <- fread(ui_files[1], sep="\t")
  eid2mz <- user_input[, .(EID, mz=round(as.numeric(get("m/z")),4))][!is.na(mz)]

  # Pathway table: explode each pathway's comma-separated EC id list into one row per id
  pathways <- fread(pa_files[1], sep="\t")
  setnames(pathways, c("pathway","p-value","overlap_EmpiricalCompounds (id)"), c("pathway","mcg_p","ecids"), skip_absent=TRUE)
  pathway_members <- pathways[, .(EID = trimws(unlist(strsplit(ecids, ",")))), by=pathway]

  merge(pathway_members, eid2mz, by="EID", allow.cartesian=TRUE)
}

# compartment -> list of run dirs (pos+neg) under a base
runs_for <- function(base, token){
  d <- list.dirs(paste0(root, base), recursive=FALSE)
  d[grepl(token, d)]
}

COMPARTMENTS <- list(
  list(lab="Maternal plasma",        ds="MaternalPlasma",       vv="pn12", mzt=plmz, base="results/mummichog_output_adjusted", tok="MaternalPlasma_pn12"),
  list(lab="Postnatal maternal VAMS",ds="VamsPostnatalMaternal",vv="pn56", mzt=vmz,  base="results/mummichog_output_adjusted", tok="VamsPostnatalMaternal_pn56"),
  list(lab="Infant VAMS pn12",       ds="VamsPostnatalInfant",  vv="pn12", mzt=vmz,  base="results/mummichog_output_adjusted", tok="VamsPostnatalInfant_pn12"),
  list(lab="Infant VAMS pn34",       ds="VamsPostnatalInfant",  vv="pn34", mzt=vmz,  base="results/mummichog_output_adjusted", tok="VamsPostnatalInfant_pn34"),
  list(lab="Infant VAMS pn56",       ds="VamsPostnatalInfant",  vv="pn56", mzt=vmz,  base="results/mummichog_output_adjusted", tok="VamsPostnatalInfant_pn56"))

out <- list()
for (cc in COMPARTMENTS){
  st <- sig_tab(cc$ds, cc$vv, cc$mzt)
  mem <- rbindlist(lapply(runs_for(cc$base, cc$tok), parse_run), fill=TRUE)
  if(!nrow(mem)) next
  mem <- merge(unique(mem[, .(pathway, mz)]), st, by="mz")        # member features with our signed stat
  # For each pathway, test whether its members lean up or down with one binomial
  # sign test on the count of up-features (single null -> honest p-value).
  agg <- mem[, {
      mean_signed_stat <- mean(signed)
      nonzero <- signed[signed != 0]           # drop signed==0 (est==0 or p==1): uninformative for a sign test
      n_nonzero <- length(nonzero); n_up <- sum(nonzero > 0)
      p <- if(n_nonzero>=4) signif(binom.test(n_up, n_nonzero, 0.5)$p.value,2) else NA_real_
      .(n_members=n_nonzero, frac_up=round(n_up/n_nonzero,2), mean_signed=round(mean_signed_stat,2),
        direction=fifelse(mean_signed_stat>0,"up",fifelse(mean_signed_stat<0,"down","flat")), sign_test_p=p)
    }, by=pathway]
  agg[, compartment := cc$lab]
  out[[length(out)+1]] <- agg
}
res <- rbindlist(out, fill=TRUE)[order(compartment, sign_test_p)]
fwrite(res, paste0(root, "results/signed_pathway_direction.csv"))

KEY <- "Fatty Acid Biosynthesis|Arachidonic|Leukotriene|Ascorbate|Carnitine|Prostaglandin"
cat("=== SIGNED DIRECTION for the key §5/§6 pathways (single-null sign test) ===\n")
print(res[grepl(KEY, pathway, ignore.case=TRUE) & n_members>=4,
          .(compartment, pathway, n_members, frac_up, mean_signed, direction, sign_test_p)], nrow=60)
cat("\nSaved results/signed_pathway_direction.csv (", nrow(res), "compartment x pathway rows )\n")
