# =============================================================================
# 46-build-table-s8-cross-compartment.R
#
# Builds Table S8, cross-compartment BEP-responsive features in MISAME-III
# (covariate-adjusted combined-arm ATEs, per dataset x visit BH q-values):
#   Panel A - named metabolite features increased in both maternal and infant blood
#             (postnatal VAMS, 5-6 mo; matched by shared feature id)
#   Panel B - octenoylcarnitine (m/z 286.202) in milk, maternal plasma, maternal
#             blood and infant blood at each postnatal visit
#   Panel C - selenoproteins (SELENOP, GPX3) raised in both milk and maternal blood
#
# Inputs : results/adjusted_combined_arms_intervention_effects_untargeted_results_clean_ATE.RDS
#          results/adjusted_combined_arms_intervention_effects_proteomics_results_clean_ATE.RDS
#          results/blood_compartment_adjusted_combined_arms_intervention_effects_results_clean.RDS
#          results/cross_compartment_proteome_overlap_adjusted.csv (script 13, BLOOD_ADJUST = TRUE)
# Outputs: results/table_s8_cross_compartment.csv (Panel A)
#          results/table_s8_fragment.md (Panels A-C as markdown tables)
# [needs restricted data]
# =============================================================================
suppressMessages({ library(data.table) })
root <- paste0(here::here(), "/")
low  <- function(x) tolower(as.character(x))

milkC  <- as.data.table(readRDS(paste0(root,"results/adjusted_combined_arms_intervention_effects_untargeted_results_clean_ATE.RDS")))
bloodC <- as.data.table(readRDS(paste0(root,"results/blood_compartment_adjusted_combined_arms_intervention_effects_results_clean.RDS")))
# Panel C uses the covariate-adjusted blood proteome (script 13 run with
# BLOOD_ADJUST = TRUE), matching Panels A/B and the adjusted milk proteome.
prot   <- fread(paste0(root,"results/cross_compartment_proteome_overlap_adjusted.csv"))

# Panel A rows: putatively annotated postnatal-VAMS features with their m/z and ion
# mode; "*" marks isobaric-ambiguous annotations (more than one candidate compound).
named <- data.table(
  fid = c("vam_1005523","vam_2001392","vam_2004310","vam_2002077","vam_2002332",
          "vam_1003290","vam_2003705","vam_2000289","vam_2000920","vam_1001196"),
  name = c("Octenoylcarnitine*","4-Hydroxy-2-oxoglutarate","Citrate / Isocitrate*",
           "Perillic acid*","Aconitic acid","N-Acetylserotonin","Dehydroascorbate*",
           "Methylmalonate","Threonic acid","N-Methylhistidine"),
  mz = c("286.202","142.998","228.994","165.092","173.009","220.118","208.986","99.009","135.030","153.066"),
  mode = c("pos","neg","neg","neg","neg","pos","neg","neg","neg","pos"))

# table-cell formatters: fq() renders a q-value ("<0.001" below threshold, else 2 sig figs);
# cell() renders "+est (q)" for a filled cell or an em dash when the effect is missing.
fq <- function(q) ifelse(is.na(q), "n/a", ifelse(q < 1e-3, "<0.001", formatC(q, format="g", digits=2)))
cell <- function(est, q) ifelse(is.na(est), "—", sprintf("%+.2f (%s)", est, fq(q)))
# getb(): look up one blood feature's ATE (est + adjusted p) for a given dataset/feature/visit
getb <- function(ds, k, v) {
  r <- bloodC[measure=="ATE" & dataset==ds & low(biomarker)==k & visit==v]
  if (nrow(r)==0) list(est=NA_real_, q=NA_real_) else list(est=r$est[1], q=r$pval_adj[1])
}

# ---- Panel A: maternal & infant blood (VAMS pn56) ----
A <- named[, {
  inf <- getb("VamsPostnatalInfant", fid, "pn56")
  mat <- getb("VamsPostnatalMaternal", fid, "pn56")
  .(name, mz, mode,
    maternal = cell(mat$est, mat$q), infant = cell(inf$est, inf$q),
    inf_est = inf$est)
}, by=seq_len(nrow(named))][order(-inf_est)]
A[, seq_len := NULL][, inf_est := NULL]
setnames(A, c("name","mz","mode","maternal","infant"),
         c("Putative annotation","m/z","Mode","Maternal blood, ΔSD (q)","Infant blood, ΔSD (q)"))

# ---- Panel B: octenoylcarnitine across compartments ----
oc_milk <- milkC[measure=="ATE" & study=="Misame" & low(biomarker)=="rlc_pos_mtb_3033006" & visit=="1-2 mo."]
oc_plas <- getb("MaternalPlasma","rlc_pos_mtb_2679130","pn12")
# Infant rows are shown at every postnatal visit, to 3 decimals: the 1-2 and 5-6 month
# estimates coincide at 2 decimals (+1.044 vs +1.043) and would otherwise look like one
# estimate under two visit labels. Fig. 6D plots the largest (3-4 mo.).
cell3 <- function(est, q) ifelse(is.na(est), "—", sprintf("%+.3f (%s)", est, fq(q)))
inf_oc <- function(v) { r <- getb("VamsPostnatalInfant","vam_1005523",v); cell3(r$est, r$q) }
B <- data.table(
  Compartment = c("Breast milk (1–2 mo.)","Maternal plasma (1–2 mo.)",
                  "Maternal blood, VAMS (5–6 mo.)","Infant blood, VAMS (1–2 mo.)",
                  "Infant blood, VAMS (3–4 mo.)","Infant blood, VAMS (5–6 mo.)"),
  `ΔSD (q)` = c(
    cell(oc_milk$est[1], oc_milk$pval_adj[1]),
    cell(oc_plas$est, oc_plas$q),
    cell(getb("VamsPostnatalMaternal","vam_1005523","pn56")$est, getb("VamsPostnatalMaternal","vam_1005523","pn56")$q),
    inf_oc("pn12"), inf_oc("pn34"), inf_oc("pn56")))

# ---- Panel C: selenoproteins (milk & maternal blood proteome) ----
psel <- prot[uniprot %in% c("P49908","P22352")][order(uniprot, -blood_est)]
psel <- psel[, .SD[1], by=uniprot]   # one row/protein (best blood hit)
# The overlap file carries raw p-values only (milk_pval/blood_pval), so the BH q-values
# are looked up in the same result files script 13 read (milk: adjusted combined-arm
# proteome; blood: bloodC, the adjusted blood results above).
milkP  <- as.data.table(readRDS(paste0(root,"results/adjusted_combined_arms_intervention_effects_proteomics_results_clean_ATE.RDS")))
milk_q  <- mapply(function(u, v) milkP[measure=="ATE" & study=="Misame" & toupper(biomarker)==u & visit==v, pval_adj][1],
                  psel$uniprot, psel$milk_visit)
blood_q <- mapply(function(u, d, v) bloodC[measure=="ATE" & dataset==d & toupper(biomarker)==u & visit==v, pval_adj][1],
                  psel$uniprot, psel$blood_dataset, psel$blood_visit)
C <- data.table(
  `Protein (gene, UniProt)` = ifelse(psel$uniprot=="P49908","Selenoprotein P (SELENOP, P49908)",
                                     "Glutathione peroxidase 3 (GPX3, P22352)"),
  `Milk, ΔSD (q)`          = mapply(cell, psel$milk_est, milk_q),
  `Maternal blood, ΔSD (q)`= mapply(cell, psel$blood_est, blood_q))

# ---- write ----
fwrite(rbind(
  data.table(panel="A", A[, lapply(.SD, as.character)]),
  fill=TRUE), paste0(root,"results/table_s8_cross_compartment.csv"))

# render a data.table as a GitHub-flavoured markdown table (header / separator / body rows)
md <- function(dt) {
  header    <- paste0("| ", paste(names(dt), collapse=" | "), " |")
  separator <- paste0("|", paste(rep("---", ncol(dt)), collapse="|"), "|")
  body      <- apply(dt, 1, function(r) paste0("| ", paste(r, collapse=" | "), " |"))
  paste(c(header, separator, body), collapse="\n")
}
out <- paste0(root,"results/table_s8_fragment.md")
writeLines(c("### A. Metabolites increased in maternal AND infant blood (VAMS, 5-6 mo.)","",
             md(A),"","### B. Octenoylcarnitine across compartments and infant visits","",
             md(B),"","### C. Selenoproteins concordant in milk & maternal blood","",
             md(C)), out)
cat("wrote", out, "\n\n"); cat(readLines(out), sep="\n"); cat("\n")
