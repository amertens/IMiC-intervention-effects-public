# =============================================================================
# 31-detectability-proxy.R
#
# PROXY detectability sensitivity analysis for the headline infant features
# (vam_1005523, vam_1005524, vam_1003290). We do NOT have a missingness mask or
# pre-imputation data (everything is ImputedCappedScaled), so this is a proxy,
# pending Kim's per-feature detection rates (memo §9-A1). Question: is the BEP
# "increase" in infant blood driven by more BEP infants escaping the imputed
# detection floor (detectability), or by a shift among already-detected samples
# (abundance)? We compare by-arm distributions, the fraction near the floor, and
# re-estimate the contrast among detected-only (above-floor) infants.
#
# Sample->arm mapping replicates src/1 data prep/7-blood-compartment-prep.R.
# Out: results/detectability_proxy_infant_carnitine.csv
# =============================================================================
suppressMessages({library(data.table); library(haven)})
root <- paste0(here::here(), "/")
add  <- paste0(root, "data/additional datasets/")

# --- arm crosswalk (idBiospe -> 2x2 arm), exactly as the prep script ----------
ak <- read_dta(paste0(add, "metadata_for_sharing.dta"))
ak <- data.table(idBiospe = as.numeric(ak$idbs),
                 code_chr = as.character(haven::as_factor(ak$code_bep_n)))
ak <- ak[!is.na(idBiospe)]
ak[, bep_postpartum := as.integer(grepl("post:BEP", code_chr, fixed = TRUE))]
ak[, bep_prenatal   := as.integer(grepl("pre:BEP",  code_chr, fixed = TRUE))]
ak[, arm := fifelse(bep_prenatal==0 & bep_postpartum==0, "Control",
             fifelse(bep_prenatal==0 & bep_postpartum==1, "IFA/BEP",
              fifelse(bep_prenatal==1 & bep_postpartum==0, "BEP/IFA", "BEP/BEP")))]
ak <- unique(ak[, .(idBiospe, arm)])

# --- carnitine extract (features x samples) -> long infant table --------------
ex <- fread(paste0(root, "results/_carnitine_infant_extract.csv"))
setnames(ex, 1, "feature"); ex[, feature := gsub('"','',feature)]
samples <- names(ex)[-1]
# Sample column names encode metadata as "<timepoint>_<..>_<idBiospe>_<dyadcode>";
# the 4th token is "e" for the infant ("enfant") and otherwise the mother.
parse_sample_ids <- function(x){ p <- tstrsplit(x, "_", fixed=TRUE)
  data.table(sample=x, timePoint=p[[1]], idBiospe=as.numeric(p[[3]]),
             dyad=fifelse(p[[4]]=="e","infant","mother")) }
meta <- parse_sample_ids(samples)
long <- melt(ex, id.vars="feature", variable.name="sample", value.name="val")
long[, sample := as.character(sample)]
long <- merge(long, meta, by="sample")
long <- merge(long, ak, by="idBiospe", all.x=TRUE)
inf  <- long[dyad=="infant" & !is.na(arm)]

# --- per-feature detection floor (across ALL postnatal samples for that feat) --
floor_tab <- long[, .(fmin=min(val), frng=max(val)-min(val)), by=feature]
inf <- merge(inf, floor_tab, by="feature")
inf[, near_floor := val <= fmin + 0.10*frng]          # within 10% of range above min
inf[, detected   := !near_floor]                       # proxy "detected"

contrast <- function(d, exposed){                      # mean diff exposed - Control, Wilcoxon p
  a <- d[arm==exposed, val]; c0 <- d[arm=="Control", val]
  if(length(a)<3 || length(c0)<3) return(list(diff=NA, p=NA, na=length(a), nc=length(c0)))
  list(diff=round(mean(a)-mean(c0),3), p=signif(suppressWarnings(wilcox.test(a,c0)$p.value),2),
       na=length(a), nc=length(c0))
}

out <- list()
for (ft in unique(inf$feature)) for (vv in c("pn12","pn34","pn56")) {
  d <- inf[feature==ft & timePoint==vv]
  if(!nrow(d)) next
  # by-arm summary (printed)
  cat(sprintf("\n== %s  %s ==\n", ft, vv))
  print(d[, .(n=.N, mean=round(mean(val),3), med=round(median(val),3),
              min=round(min(val),3), max=round(max(val),3),
              pct_near_floor=round(100*mean(near_floor))), by=arm][order(arm)])
  for (ex_arm in c("IFA/BEP","BEP/BEP")) {
    all_c <- contrast(d, ex_arm)
    det_c <- contrast(d[detected==TRUE], ex_arm)              # detected-only
    fl    <- d[, .(pf=mean(near_floor)), by=arm]
    out[[length(out)+1]] <- data.table(
      feature=ft, visit=vv, exposed=ex_arm,
      n_exp=all_c$na, n_ctrl=all_c$nc,
      diff_all=all_c$diff, p_all=all_c$p,
      pct_floor_ctrl=round(100*fl[arm=="Control", pf]),
      pct_floor_exp =round(100*fl[arm==ex_arm,  pf]),
      n_exp_detected=det_c$na, n_ctrl_detected=det_c$nc,
      diff_detected=det_c$diff, p_detected=det_c$p)
  }
}
res <- rbindlist(out)
fwrite(res, paste0(root, "results/detectability_proxy_infant_carnitine.csv"))
cat("\n\n=== DETECTABILITY PROXY SUMMARY (effect on all infants vs detected-only) ===\n")
print(res)
cat("\nIf diff_detected collapses toward 0 vs diff_all, and pct_floor_exp << pct_floor_ctrl,\n")
cat("the increase is largely a floor-escape (detectability) effect. Proxy only; Kim's mask is definitive.\n")
