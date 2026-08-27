# =============================================================================
# 33-signed-pathway-direction-milk.R
#
# Milk side of the single-null signed confirmation (companion to script 32, which
# did the blood compartments). Confirms the §5/§6 milk directions (lipid pathways
# DOWN in milk) with a binomial sign test on each pathway's member features,
# instead of the two-run directional Mummichog. Milk Mummichog runs are per visit
# (results/mummichog_output/Milk_{1421days,12mo,34mo}_{pos,neg}); milk effects are
# the combined-arms (pooled postnatal-BEP) ATEs, matched to each run's visit.
# Out: results/signed_pathway_direction_milk.csv
# =============================================================================
suppressMessages({library(data.table)})
root <- paste0(here::here(), "/")
source(paste0(root, "src/2 analysis/_blood_helpers.R"))
up <- function(x) toupper(as.character(x))

milkC <- as.data.table(readRDS(paste0(root,"results/adjusted_combined_arms_intervention_effects_untargeted_results_clean_ATE.RDS")))
milkC <- milkC[study=="Misame" & measure=="ATE"]
cat("milk combined-arms visits:\n"); print(unique(as.character(milkC$visit)))
mkmz <- mzrt_milk(misame_only=FALSE)

# Our milk signed statistic, keyed by m/z: sign = BEP direction, magnitude = -log10(p).
# pmax(pval, 1e-300) floors p away from 0 so log10 stays finite.
sig_tab <- function(vv){
  d <- milkC[visit==vv]; d[, feature := up(biomarker)]
  d <- merge(d, mkmz[, .(feature, mz)], by="feature")
  d[, signed := sign(est) * -log10(pmax(pval, 1e-300))]
  d[is.finite(signed) & is.finite(mz), .(mz=round(mz,4), signed)][!duplicated(mz)]
}

# Parse one Mummichog run folder into pathway -> member-feature-m/z (see script 32).
# Bridge is EC id -> m/z from the userInput table, then attach to each pathway.
parse_run <- function(dir){
  tables_dir <- paste0(dir, "/tables")
  ui_files <- list.files(tables_dir, pattern="userInput_to_EmpiricalCompounds", full.names=TRUE)
  pa_files <- list.files(tables_dir, pattern="mcg_pathwayanalysis", full.names=TRUE)
  if(!length(ui_files) || !length(pa_files)) return(NULL)

  user_input <- fread(ui_files[1], sep="\t")
  eid2mz <- user_input[, .(EID, mz=round(as.numeric(get("m/z")),4))][!is.na(mz)]

  pathways <- fread(pa_files[1], sep="\t")
  setnames(pathways, c("overlap_EmpiricalCompounds (id)"), c("ecids"), skip_absent=TRUE)
  pathway_members <- pathways[, .(EID = trimws(unlist(strsplit(ecids, ",")))), by=pathway]

  merge(pathway_members, eid2mz, by="EID", allow.cartesian=TRUE)
}
runs_for <- function(token){ d <- list.dirs(paste0(root,"results/mummichog_output"), recursive=FALSE); d[grepl(token, d)] }

# map run token -> milk visit string (keyword match against milkC visits)
vlevels <- unique(as.character(milkC$visit))
pick_visit <- function(rx) { m <- vlevels[grepl(rx, vlevels, ignore.case=TRUE)]; if(length(m)) m[1] else NA }
RUNS <- list(
  list(tok="Milk_12mo",     visit=pick_visit("1.?2|1-2")),
  list(tok="Milk_34mo",     visit=pick_visit("3.?4|3-4")),
  list(tok="Milk_1421days", visit=pick_visit("14|21|day|acco|colos")))

out <- list()
for (r in RUNS){
  if(is.na(r$visit)){ cat("no milk visit matched for", r$tok, "\n"); next }
  st  <- sig_tab(r$visit)
  mem <- rbindlist(lapply(runs_for(r$tok), parse_run), fill=TRUE)
  if(!nrow(mem)) next
  mem <- merge(unique(mem[, .(pathway, mz)]), st, by="mz")
  # Per pathway: net direction from the mean signed stat, tested with one binomial
  # sign test on the count of up-features (all members counted here).
  agg <- mem[, { n_members <- .N; n_up <- sum(signed>0); mean_signed_stat <- mean(signed)
      .(visit=r$visit, n_members=n_members, frac_up=round(n_up/n_members,2), mean_signed=round(mean_signed_stat,2),
        direction=fifelse(mean_signed_stat>0,"up",fifelse(mean_signed_stat<0,"down","flat")),
        sign_test_p=if(n_members>=4) signif(binom.test(n_up,n_members,0.5)$p.value,2) else NA_real_) }, by=pathway]
  out[[length(out)+1]] <- agg
}
res <- rbindlist(out, fill=TRUE)
fwrite(res, paste0(root,"results/signed_pathway_direction_milk.csv"))

KEY <- "Fatty Acid Biosynthesis|Arachidonic|Leukotriene|Ascorbate|Carnitine|Prostaglandin"
cat("\n=== MILK signed direction, key §5/§6 pathways (single-null sign test) ===\n")
print(res[grepl(KEY, pathway, ignore.case=TRUE) & n_members>=4,
          .(visit, pathway, n_members, frac_up, mean_signed, direction, sign_test_p)][order(pathway, visit)], nrow=60)
cat("\nSaved results/signed_pathway_direction_milk.csv\n")
