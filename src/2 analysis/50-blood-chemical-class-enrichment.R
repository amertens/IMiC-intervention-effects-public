# =============================================================================
# 50-blood-chemical-class-enrichment.R
#
# "MSEA-style" enrichment for the BLOOD compartments, using CHEMICAL CLASS
# instead of pathway names (blood untargeted features are mass-based; only ~1%
# carry a name, so a name-based MSEA is not feasible - see 51 pilot). The
# Biocrates/Sapient "Compound Class" annotation, however, covers the V3 VAMS
# catalogue, so we can ask a class-level enrichment question directly:
#
#   Among measured+classed features in a compartment, are the BEP-responsive
#   ones (combined-arms ATE, p<0.05) over-represented in any chemical class?
#
# Compartments (all three requested):
#   VamsPostnatalInfant   - vam_ ids -> Compound Class DIRECT (cleanest)
#   VamsPostnatalMaternal - vam_ ids -> Compound Class DIRECT
#   MaternalPlasma        - rLC_ (V1) ids -> class BORROWED by m/z match to the
#                           V3 class table (25 ppm + ion mode); APPROXIMATE, flagged.
#
# Fisher exact per class (foreground = responsive & classed, background = classed),
# BH-FDR across classes within compartment. Also reports direction (% up).
# EXPLORATORY: p<0.05 foreground is permissive (matches the milk untargeted ORA
# convention); class annotation is putative; plasma borrowing is isobaric.
# =============================================================================
suppressMessages({library(data.table); library(biotmle)})
root <- paste0(here::here(), "/")
add  <- paste0(root, "data/additional datasets/")
PPM <- 25; MIN_CLASS <- 5; SIG_P <- 0.05

# ---- BEP combined-arms ATE per compartment ----------------------------------
r <- as.data.table(readRDS(paste0(root,
  "results/blood_compartment_adjusted_combined_arms_intervention_effects_results_clean.RDS")))
ate <- r[measure == "ATE" & contrast == "BEP" & !is.na(pval)]
# collapse to one row per (dataset, biomarker): responsive if p<0.05 at any visit
feat <- ate[order(pval), .(pval = pval[1], est = est[1],
                           responsive = as.integer(min(pval) < SIG_P)),
            by = .(dataset, biomarker)]

# ---- class annotation -------------------------------------------------------
vam <- fread(paste0(add, "metabolite_description_vam_with_global_id.csv"))
setnames(vam, c("feature_label", "Compound Class", "ionization_mode"),
              c("biomarker", "class", "mode"), skip_absent = TRUE)
vam_cls <- vam[!is.na(class) & class != "" & class != "-",
               .(biomarker, mz, mode = tolower(mode), class)]

# plasma (V1) m/z + ion mode, to borrow class by mass
pl <- fread(paste0(add, "ProcessedDataMISAME3_plasma.csv"), select = c("MZ","RT","Metabolite_Feature_Label"))
setnames(pl, c("mz","rt","biomarker"))
pl[, mode := fifelse(grepl("_pos_", biomarker), "positive",
              fifelse(grepl("_neg_", biomarker), "negative", NA_character_))]

borrow_class_by_mz <- function(q_mz, q_mode) {   # nearest same-mode classed V3 feature within PPM
  cand <- vam_cls[mode == q_mode]
  if (!nrow(cand)) return(NA_character_)
  tol <- q_mz * PPM * 1e-6
  hit <- cand[abs(mz - q_mz) <= tol]
  if (!nrow(hit)) return(NA_character_)
  hit[which.min(abs(mz - q_mz))]$class
}

# ---- assemble per-compartment feature table with class ----------------------
attach_class <- function(ds) {
  f <- feat[dataset == ds]
  if (grepl("^Vams", ds)) {                     # V3: direct id join
    f <- merge(f, vam_cls[, .(biomarker, class)], by = "biomarker", all.x = TRUE)
    f[, class_src := "direct (V3 id)"]
  } else if (ds == "MaternalPlasma") {          # V1: borrow by m/z
    f <- merge(f, pl[, .(biomarker, mz, mode)], by = "biomarker", all.x = TRUE)
    f[, class := mapply(function(m, io) if (is.na(m) || is.na(io)) NA_character_
                        else borrow_class_by_mz(m, io), mz, mode)]
    f[, class_src := "borrowed (V1->V3 m/z, isobaric)"]
  }
  f[!is.na(class)]
}

# ---- Fisher enrichment per class --------------------------------------------
enrich <- function(f, ds) {
  n_resp <- sum(f$responsive); n_tot <- nrow(f)
  if (n_resp < 3) return(NULL)
  cls <- names(which(table(f$class) >= MIN_CLASS))
  rbindlist(lapply(cls, function(cc) {
    inC <- f$class == cc
    a <- sum(inC & f$responsive == 1); b <- sum(!inC & f$responsive == 1)
    c <- sum(inC & f$responsive == 0); d <- sum(!inC & f$responsive == 0)
    ft <- fisher.test(matrix(c(a, b, c, d), 2, byrow = TRUE))
    up <- if (a > 0) mean(f[inC & responsive == 1]$est > 0) else NA_real_
    data.table(compartment = ds, class = cc, n_class = a + c, n_responsive = a,
               pct_class_responsive = round(100 * a / (a + c), 1),
               OR = unname(ft$estimate), p = ft$p.value, pct_up = round(100 * up, 0))
  }))
}

comps <- c("VamsPostnatalInfant", "VamsPostnatalMaternal", "MaternalPlasma")
res <- rbindlist(lapply(comps, function(ds) {
  f <- attach_class(ds)
  cat(sprintf("%-22s classed features: %d (responsive p<.05: %d) [%s]\n",
              ds, nrow(f), sum(f$responsive), f$class_src[1]))
  enrich(f, ds)
}), fill = TRUE)
res[, fdr := p.adjust(p, "BH"), by = compartment]
setorder(res, compartment, p)
fwrite(res, paste0(root, "results/blood_chemical_class_enrichment.csv"))

cat("\n==== BLOOD chemical-class enrichment (BEP-responsive features, combined-arms) ====\n")
for (ds in comps) {
  cat("\n---", ds, "---\n")
  x <- res[compartment == ds, .(class, n_class, n_resp = n_responsive,
           `%resp` = pct_class_responsive, OR = round(OR, 2),
           p = signif(p, 2), fdr = signif(fdr, 2), `%up` = pct_up)]
  if (nrow(x)) print(head(x, 10), row.names = FALSE) else cat("  (too few responsive/classed features)\n")
}
cat("\nSaved results/blood_chemical_class_enrichment.csv\n")
