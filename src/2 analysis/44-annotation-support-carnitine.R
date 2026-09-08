# =============================================================================
# 44-annotation-support-carnitine.R
#
# LOCAL, OFFLINE annotation support for the transfer feature at m/z 286.2013
# (catalogued octenoylcarnitine; Kim read "octanoyl"). Uses ONLY the m/z-RT
# feature catalogues already on disk -- no raw MS2, no internet, no Sapient file.
# It cannot reach an MSI Level 1/2 identification (that needs MS2 / a standard),
# but it bolsters the C15H27NO4 (octenoylcarnitine) molecular-formula assignment
# and the acylcarnitine-class reading with three orthogonal, reproducible checks:
#
#   (1) Isotopologue / adduct co-features: for the [M+H]+ target, is there a
#       co-eluting 13C isotopologue (+1.00335), Na adduct (+21.9819 vs [M+H]+),
#       and in-source water-loss ion ([M+H-H2O]+)? Their presence at the right
#       mass + same RT corroborates the formula.
#   (2) Acylcarnitine diagnostic fragments: co-eluting in-source m/z 85.0284
#       (C4H5O2+) and 60.0808 (C3H10N+) -> class marker (weak; MS1 in-source only).
#   (3) Homolog RT order: place the target on the acylcarnitine C0..C18 [M+H]+
#       ladder per platform and check it elutes where an octenoyl (C8:1) should.
#
# Out: results/carnitine_annotation_support.csv   (checks 1-2, per platform)
#      results/carnitine_homolog_rt_order.csv      (check 3, per platform)
# =============================================================================

suppressMessages({library(data.table); library(here)})
root <- paste0(here::here(), "/")
source(paste0(root, "src/2 analysis/_blood_helpers.R"))

# -- exact masses -------------------------------------------------------------
# mE = monoisotopic element masses; formulas below are named integer-count vectors
# c(C=, H=, N=, O=). mono() sums a formula to its neutral mass; mh() adds a proton
# to get the observed [M+H]+ ion m/z.
mE  <- c(C = 12.0000000, H = 1.0078250319, N = 14.0030740052, O = 15.9949146221)
PROTON <- 1.0072764669; NA_ION <- 22.9892207; H2O <- 2*mE["H"] + mE["O"]
D13C   <- 1.0033548378                       # 13C - 12C mass gap (isotopologue spacing)
mono   <- function(f) sum(mapply(function(el, n) mE[[el]] * n, names(f), f))  # neutral monoisotopic
mh     <- function(f) mono(f) + PROTON                                        # [M+H]+

# octenoylcarnitine C15H27NO4 (the catalogue call), the compound under test
oct1 <- c(C = 15, H = 27, N = 1, O = 4)
TGT  <- mh(oct1)                                                              # 286.2013

# ion species we expect to co-elute with the [M+H]+ target if it really is this compound
species <- data.table(
  species  = c("[M+H]+ (target)", "13C isotopologue [M+H+1]", "[M+2]",
               "[M+Na]+", "[M+H-H2O]+", "acylcarnitine frag C4H5O2+", "acylcarnitine frag C3H10N+"),
  exp_mz   = c(TGT, TGT + D13C, TGT + 2*D13C,
               mono(oct1) + NA_ION, TGT - H2O,
               (4*mE["C"] + 5*mE["H"] + 2*mE["O"]) - 0.0005486,      # 85.0284
               (3*mE["C"] + 10*mE["H"] + mE["N"]) - 0.0005486))      # 60.0808

# -- acylcarnitine homolog ladder (neutral formulas, verified) ----------------
acns <- data.table(
  label = c("C0 (free carnitine)","C2","C3","C4","C5","C6",
            "C8:1 (octenoyl, TARGET)","C8 (octanoyl)","C10","C12","C14","C16","C18"),
  C = c(7,9,10,11,12,13,15,15,17,19,21,23,25),
  H = c(15,17,19,21,23,25,27,29,33,37,41,45,49),
  N = 1, O = c(3,4,4,4,4,4,4,4,4,4,4,4,4))
acns[, exp_mh := mapply(function(c,h,n,o) mh(c(C=c,H=h,N=n,O=o)), C, H, N, O)]

# -- platform catalogues (all positive-mode features) -------------------------
loaders <- list(
  Milk_V1          = function() mzrt_milk(misame_only = TRUE),
  MaternalPlasma_V1= function() mzrt_rlc("ProcessedDataMISAME3_plasma.csv"),
  PrenatalVAMS_V1  = function() mzrt_rlc("ProcessedDataMISAME3_VAMS.csv"),
  PostnatalVAMS_V3 = function() mzrt_vams())

PPM   <- 15                 # within-platform accurate-mass tolerance for reporting
within_ppm <- function(mz, target, ppm = PPM) abs(mz - target)/target*1e6 <= ppm

sup_rows <- list(); hom_rows <- list()

for (pl in names(loaders)) {
  cat("\n================  ", pl, "  ================\n")
  tab <- tryCatch(loaders[[pl]](), error = function(e) { cat("  LOAD FAILED:", conditionMessage(e), "\n"); NULL })
  if (is.null(tab) || !nrow(tab)) next
  pos <- tab[mode == "positive" & !is.na(mz)]
  cat(sprintf("  %d positive-mode features (%d total)\n", nrow(pos), nrow(tab)))

  # locate the target [M+H]+ (may be >1 feature, e.g. the vam_1005523/1005524 pair)
  tgt <- pos[within_ppm(mz, TGT, 25)][order(abs(mz - TGT))]
  if (!nrow(tgt)) { cat("  target 286.2013 NOT found within 25 ppm on this platform\n"); }
  else {
    cat(sprintf("  target feature(s) at ~286.2013: %s\n",
                paste(sprintf("%s (mz=%.4f, rt=%.3f, %.1f ppm)", tgt$feature, tgt$mz, tgt$rt,
                              (tgt$mz - TGT)/TGT*1e6), collapse = "; ")))
  }
  rt0 <- if (nrow(tgt)) tgt$rt[1] else NA_real_        # anchor RT = nearest-ppm target feature

  # (1)-(2) isotope / adduct / fragment co-features: for each expected ion species,
  # find catalogue features within tolerance of its m/z, then pick the one that best
  # co-elutes with the target (smallest RT gap) -- or, absent an RT anchor, closest in mass.
  for (i in seq_len(nrow(species))) {
    expected_mz <- species$exp_mz[i]
    hit <- pos[within_ppm(mz, expected_mz)]
    if (nrow(hit)) hit[, drt := rt - rt0]                 # RT gap from the target anchor
    if (nrow(hit)) {
      hit <- if (!is.na(rt0)) hit[order(abs(drt))] else hit[order(abs(mz - expected_mz))]
      best <- hit[1]
      sup_rows[[length(sup_rows)+1]] <- data.table(
        platform = pl, species = species$species[i], exp_mz = round(expected_mz,4),
        matched_feature = best$feature, matched_mz = round(best$mz,4),
        ppm = round((best$mz - expected_mz)/expected_mz*1e6,1),
        matched_rt = round(best$rt,3), drt = if(!is.na(rt0)) round(best$rt - rt0,3) else NA_real_,
        n_within_ppm = nrow(hit))
    } else {
      # no feature within tolerance -> record an explicit "not found" row
      sup_rows[[length(sup_rows)+1]] <- data.table(
        platform = pl, species = species$species[i], exp_mz = round(expected_mz,4),
        matched_feature = NA_character_, matched_mz = NA_real_, ppm = NA_real_,
        matched_rt = NA_real_, drt = NA_real_, n_within_ppm = 0L)
    }
  }

  # (3) homolog RT ladder, nearest-ppm feature per acylcarnitine [M+H]+
  for (i in seq_len(nrow(acns))) {
    hit <- pos[within_ppm(mz, acns$exp_mh[i], 10)][order(abs(mz - acns$exp_mh[i]))]
    hom_rows[[length(hom_rows)+1]] <- data.table(
      platform = pl, acn = acns$label[i], exp_mh = round(acns$exp_mh[i],4),
      matched_feature = if(nrow(hit)) hit$feature[1] else NA_character_,
      matched_mz = if(nrow(hit)) round(hit$mz[1],4) else NA_real_,
      ppm = if(nrow(hit)) round((hit$mz[1]-acns$exp_mh[i])/acns$exp_mh[i]*1e6,1) else NA_real_,
      rt = if(nrow(hit)) round(hit$rt[1],3) else NA_real_, n_hits = nrow(hit))
  }
}

sup <- rbindlist(sup_rows); hom <- rbindlist(hom_rows)
fwrite(sup, paste0(root, "results/carnitine_annotation_support.csv"))
fwrite(hom, paste0(root, "results/carnitine_homolog_rt_order.csv"))

cat("\n\n========== ISOTOPE / ADDUCT / FRAGMENT SUPPORT ==========\n")
print(sup[, .(platform, species, exp_mz, matched_mz, ppm, drt, n_within_ppm)])
cat("\n========== ACYLCARNITINE HOMOLOG RT LADDER ==========\n")
print(hom[, .(platform, acn, exp_mh, matched_mz, ppm, rt, n_hits)])
cat("\nSaved: results/carnitine_annotation_support.csv, results/carnitine_homolog_rt_order.csv\n")
