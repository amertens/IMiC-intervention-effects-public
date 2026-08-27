# =============================================================================
# 52-blood-class-enrichment-direction-sensitivity.R
#
# Two refinements of 50-blood-chemical-class-enrichment.R:
#  (a) DIRECTION-SPLIT foregrounds: enrich up-responsive (ATE>0, p<0.05) and
#      down-responsive (ATE<0, p<0.05) features separately, per compartment.
#  (b) PLASMA MATCHING-STRINGENCY LADDER: plasma class is borrowed from the V3
#      VAMS table by m/z (V1->V3, isobaric). Re-run the plasma enrichment under
#      progressively stricter matching to see whether the acylcarnitine hit
#      survives or is an artifact of promiscuous matching:
#        25 ppm  ->  10 ppm  ->  10 ppm + unique class  ->  + RT window.
#      NOTE: V1 (plasma) and V3 (VAMS) run different LC gradients, so their RT is
#      NOT linearly comparable (Kim). The RT filter is therefore an OVER-STRICT
#      stress test, not a correct matcher: if the signal survives even it, the
#      result is not purely a loose-match artifact; if it collapses, the plasma
#      class result cannot be verified from mass alone.
# =============================================================================
suppressMessages({library(data.table); library(biotmle)})
root <- paste0(here::here(), "/")
add  <- paste0(root, "data/additional datasets/")
MIN_CLASS <- 5; SIG_P <- 0.05

# ---- BEP combined-arms ATE, one row per (dataset, biomarker) -----------------
r <- as.data.table(readRDS(paste0(root,
  "results/blood_compartment_adjusted_combined_arms_intervention_effects_results_clean.RDS")))
ate <- r[measure == "ATE" & contrast == "BEP" & !is.na(pval)]
feat <- ate[order(pval), .(pval = pval[1], est = est[1]), by = .(dataset, biomarker)]
feat[, `:=`(up = as.integer(pval < SIG_P & est > 0),
            down = as.integer(pval < SIG_P & est < 0))]

# ---- class annotation (V3 vam_ table) ---------------------------------------
vam <- fread(paste0(add, "metabolite_description_vam_with_global_id.csv"))
setnames(vam, c("feature_label","Compound Class","ionization_mode","rt_minute"),
              c("biomarker","class","mode","rt"), skip_absent = TRUE)
vam_cls <- vam[!is.na(class) & class != "" & class != "-",
               .(biomarker, mz, rt, mode = tolower(mode), class)]

# ---- Fisher enrichment for one directional foreground -----------------------
enrich_dir <- function(f, ds, dir) {
  fg <- f[[dir]]; if (sum(fg) < 3) return(NULL)
  cls <- names(which(table(f$class) >= MIN_CLASS))
  rbindlist(lapply(cls, function(cc) {
    inC <- f$class == cc
    a <- sum(inC & fg == 1); b <- sum(!inC & fg == 1)
    c <- sum(inC & fg == 0); d <- sum(!inC & fg == 0)
    ft <- fisher.test(matrix(c(a, b, c, d), 2, byrow = TRUE))
    data.table(compartment = ds, direction = dir, class = cc,
               n_class = a + c, n_fg = a, OR = unname(ft$estimate), p = ft$p.value)
  }))
}

# =============================================================================
# (a) DIRECTION-SPLIT enrichment (VAMS = direct V3 id; plasma = 25 ppm borrow)
# =============================================================================
classify_plasma <- function(q_mz, q_rt, q_mode, ppm, uniq, rt_win) {
  cand <- vam_cls[mode == q_mode]
  cand <- cand[abs(mz - q_mz) <= q_mz * ppm * 1e-6]
  if (!is.na(rt_win)) cand <- cand[abs(rt - q_rt) <= rt_win]
  if (!nrow(cand)) return(NA_character_)
  cl <- unique(cand$class)
  if (uniq && length(cl) > 1) return(NA_character_)      # ambiguous -> drop
  if (length(cl) == 1) return(cl)
  cand[which.min(abs(mz - q_mz))]$class                  # nearest-mass class
}

# V1 (rLC catalogue) m/z sources for the borrow path. BOTH plasma and PRENATAL
# VAMS are V1, so their class must be borrowed from the V3 table by m/z (isobaric).
# (Only the two POSTNATAL VAMS compartments are V3 and join class by direct id.)
load_v1 <- function(file) {
  x <- fread(paste0(add, file), select = c("MZ","RT","Metabolite_Feature_Label"))
  setnames(x, c("mz","rt","biomarker"))
  x[, mode := fifelse(grepl("_pos_", biomarker), "positive",
              fifelse(grepl("_neg_", biomarker), "negative", NA_character_))][]
}
v1tab <- list(MaternalPlasma = load_v1("ProcessedDataMISAME3_plasma.csv"),
              VamsPrenatal   = load_v1("ProcessedDataMISAME3_VAMS.csv"))

attach_class <- function(ds, ppm = 25, uniq = FALSE, rt_win = NA) {
  f <- feat[dataset == ds]
  if (ds %in% c("VamsPostnatalMaternal","VamsPostnatalInfant")) {  # V3: direct id join
    f <- merge(f, vam_cls[, .(biomarker, class)], by = "biomarker", all.x = TRUE)
  } else if (ds %in% names(v1tab)) {                              # V1: borrow by mass
    f <- merge(f, v1tab[[ds]][, .(biomarker, mz, rt, mode)], by = "biomarker", all.x = TRUE)
    f[, class := mapply(function(m, t, io) if (is.na(m) || is.na(io)) NA_character_
                        else classify_plasma(m, t, io, ppm, uniq, rt_win), mz, rt, mode)]
  }
  f[!is.na(class)]
}

# VamsPrenatal = NEGATIVE CONTROL (V1; reflects the prenatal-BEP contrast, where
# the postnatal lipid remodeling is not expected and was empirically ~null).
comps <- c("VamsPostnatalInfant","VamsPostnatalMaternal","MaternalPlasma","VamsPrenatal")
dir_res <- rbindlist(lapply(comps, function(ds) {
  f <- attach_class(ds)
  rbind(enrich_dir(f, ds, "up"), enrich_dir(f, ds, "down"))
}), fill = TRUE)
dir_res[, fdr := p.adjust(p, "BH"), by = .(compartment, direction)]
setorder(dir_res, compartment, direction, p)
fwrite(dir_res, paste0(root, "results/blood_chemical_class_enrichment_directional.csv"))

cat("==== (a) DIRECTION-SPLIT class enrichment (top per compartment x direction) ====\n")
for (ds in comps) for (dr in c("up","down")) {
  x <- dir_res[compartment == ds & direction == dr]
  if (!nrow(x)) { cat(sprintf("\n--- %s / %s: <3 features, skipped ---\n", ds, dr)); next }
  cat(sprintf("\n--- %s / %s (n foreground=%d) ---\n", ds, dr, max(0, sum(feat[dataset==ds][[dr]]))))
  print(head(x[, .(class, n_class, n_fg, OR = round(OR,2), p = signif(p,2), fdr = signif(fdr,2))], 5), row.names = FALSE)
}

# =============================================================================
# (b) PLASMA matching-stringency ladder (focus: Acylcarnitines)
# =============================================================================
cat("\n\n==== (b) PLASMA class enrichment vs matching stringency (per direction) ====\n")
configs <- list(
  list(lab = "25 ppm (base)",            ppm = 25, uniq = FALSE, rt = NA),
  list(lab = "10 ppm",                   ppm = 10, uniq = FALSE, rt = NA),
  list(lab = "10 ppm + unique class",    ppm = 10, uniq = TRUE,  rt = NA),
  list(lab = "10 ppm + unique + RT<0.3", ppm = 10, uniq = TRUE,  rt = 0.3))
# the FDR-significant plasma (borrowed-class) hits from part (a), by direction
targets <- list(list(cls = "Acylcarnitines",             dir = "down"),
                list(cls = "Oxylipins",                  dir = "up"),
                list(cls = "Fatty acids and conjugates", dir = "up"))
lad <- rbindlist(lapply(c("MaternalPlasma", "VamsPrenatal"), function(cp) {
  rbindlist(lapply(configs, function(cf) {
    f <- attach_class(cp, cf$ppm, cf$uniq, cf$rt)
    rbindlist(lapply(targets, function(tg) {
      fg <- f[[tg$dir]]; inC <- f$class == tg$cls
      a <- sum(inC & fg == 1); b <- sum(!inC & fg == 1)
      c <- sum(inC & fg == 0); d <- sum(!inC & fg == 0)
      tgt <- sprintf("%s (%s)", tg$cls, tg$dir)
      if (a + c < 3) return(data.table(compartment = cp, target = tgt, config = cf$lab,
                       n_classed = nrow(f), n_class = a + c, n_fg = a, OR = NA_real_, p = NA_real_))
      ft <- fisher.test(matrix(c(a, b, c, d), 2, byrow = TRUE))
      data.table(compartment = cp, target = tgt, config = cf$lab, n_classed = nrow(f),
                 n_class = a + c, n_fg = a, OR = round(unname(ft$estimate), 2), p = signif(ft$p.value, 2))
    }))
  }))
}))
setorder(lad, target, compartment)
cat("(MaternalPlasma = the hits to firm; VamsPrenatal = negative control, expect null)\n")
for (tg in unique(lad$target)) { cat("\n--", tg, "--\n")
  print(lad[target == tg, .(compartment, config, n_class, n_fg, OR, p)], row.names = FALSE) }
fwrite(lad, paste0(root, "results/blood_plasma_class_match_sensitivity.csv"))
cat("\nSaved: blood_chemical_class_enrichment_directional.csv + blood_plasma_class_match_sensitivity.csv\n")
