# =============================================================================
# exploratory-proteomics-uniprot.R
# Faithful R port of Trenton's "Exploratory Outcomes (Proteomics - UniProt).Rmd".
#
# Mirrors his notebook: GO Biological-Process over-representation run DIRECTLY on
# UniProt identifiers (each measured protein counted ONCE), via a custom
# UniProt -> GO(BP) TERM2GENE map from org.Hs.eg.db + clusterProfiler::enricher.
# This avoids enrichGO's UniProt->Entrez expansion (one protein -> several genes),
# which inflates the hypergeometric counts.
#
#   Foreground : per-cell significant proteins (sig == 1) split by direction (est sign).
#   Background : the measured proteins in that study x contrast cell (GO-mapped).
#   Cells      : MISAME combined (BEP); Mumta-LW stratified (BEP+ExBf, BEP+ExBf+AZT).
#                No ELICIT proteomics. (Matches his Rmd exactly.)
#
# AUTOMATION: none needed -- his Rmd is already pure R (no metaboanalyst.ca step).
# This port is the in-repo, shared-RDS version of his analysis; it is identical in
# method to src/2 analysis/55-proteomics-go-uniprot.R (which feeds Fig 6C).
#
# FIDELITY vs his DOWNLOADED results (exploratory_proteomics_uniprot_go_enrichment_
# results.xlsx, 2026-08-10): the method reproduces his output with the foreground
# protein sets and counts IDENTICAL; p-values differ only in the ~3rd significant
# figure because the GO-BP background universe differs by ~1% (e.g. 14/1196 vs
# 14/1207 for GO:0000723) -- GO-annotation VERSION drift between his Bioconductor
# build and this one (org.Hs.eg.db/GO.db 3.20.0), NOT a code difference. FDR and
# significance calls are unchanged. Pin his annotation versions for a bit-exact run.
#
# Run from repo root:  Rscript src/trenton-ports/exploratory-proteomics-uniprot.R
# =============================================================================
suppressMessages({
  library(clusterProfiler); library(org.Hs.eg.db); library(GO.db)
  library(AnnotationDbi); library(data.table)
})
COMBINED_RDS   <- "results/combined_intervention_effects_results_combined_arms.RDS"
STRATIFIED_RDS <- "results/combined_intervention_effects_results_stratified_arms.RDS"
OUT <- "results/trenton-ports/proteomics_uniprot"

combined   <- as.data.table(readRDS(COMBINED_RDS))
stratified <- as.data.table(readRDS(STRATIFIED_RDS))

# -- all measured proteomics UniProt ids (foreground + background universe) -----
prot_ids <- unique(rbind(combined, stratified, fill = TRUE)[
  description == "Proteomics" & measure == "ATE" & !is.na(label_f) & label_f != "", label_f])
cat(sprintf("measured proteomics UniProt ids: %d\n", length(prot_ids)))

# -- UniProt -> GO(BP) map (GOALL = includes ancestor terms; each protein once) --
go_ann <- as.data.table(AnnotationDbi::select(org.Hs.eg.db, keys = prot_ids, keytype = "UNIPROT",
                                             columns = c("GOALL", "ONTOLOGYALL")))
go_ann <- unique(go_ann[ONTOLOGYALL == "BP" & !is.na(GOALL) & GOALL != "" & !is.na(UNIPROT) & UNIPROT != "",
                        .(term = GOALL, gene = UNIPROT)])
term2gene <- unique(go_ann[, .(term, gene)])
t2n <- as.data.table(AnnotationDbi::select(GO.db, keys = unique(term2gene$term), keytype = "GOID", columns = "TERM"))
term2name <- unique(t2n[!is.na(TERM) & TERM != "", .(term = GOID, name = TERM)])

# -- one enrichment cell (his construction) ------------------------------------
run_cell <- function(data, studytime_value, contrast_value, direction_value) {
  base <- data[studytime == studytime_value & contrast == contrast_value &
               description == "Proteomics" & measure == "ATE" & !is.na(label_f) & label_f != ""]
  background <- intersect(unique(base$label_f), term2gene$gene)
  fg <- base[sig == 1 & (if (direction_value == "Upregulated") est > 0 else est < 0)]
  foreground <- intersect(unique(fg$label_f), background)
  if (!length(foreground)) return(NULL)
  e <- tryCatch(enricher(gene = foreground, universe = background, pAdjustMethod = "BH",
                         pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 10, maxGSSize = 500,
                         TERM2GENE = term2gene, TERM2NAME = term2name), error = function(x) NULL)
  if (is.null(e) || !nrow(e@result)) return(NULL)
  as.data.table(e@result)[, `:=`(studytime = studytime_value, contrast = contrast_value,
                                 direction = direction_value)]
}

cells <- rbindlist(list(
  run_cell(combined,   "Misame (1-2 mo.)", "BEP",          "Downregulated"),
  run_cell(combined,   "Misame (1-2 mo.)", "BEP",          "Upregulated"),
  run_cell(stratified, "Vital (2 mo.)",    "BEP+ExBf",     "Downregulated"),
  run_cell(stratified, "Vital (2 mo.)",    "BEP+ExBf",     "Upregulated"),
  run_cell(stratified, "Vital (2 mo.)",    "BEP+ExBf+AZT", "Downregulated"),
  run_cell(stratified, "Vital (2 mo.)",    "BEP+ExBf+AZT", "Upregulated")
), fill = TRUE)

cells[, study := fifelse(grepl("Misame", studytime), "MISAME-III", "Mumta-LW")]
cells[, timepoint := fifelse(grepl("Misame", studytime), "1-2 mo.", "2 mo.")]
cells[, signed_fold_enrichment := fifelse(direction == "Upregulated", abs(FoldEnrichment), -abs(FoldEnrichment))]
setorder(cells, study, timepoint, contrast, direction, pvalue)

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
fwrite(cells, file.path(OUT, "proteomics_uniprot_go.csv"))
fwrite(cells[p.adjust < 0.05], file.path(OUT, "proteomics_uniprot_go_fdrsig.csv"))
cat(sprintf("wrote %s  (%d GO-BP terms, %d FDR<0.05; up FDR-sig: %d)\n",
            OUT, nrow(cells), nrow(cells[p.adjust < 0.05]),
            nrow(cells[p.adjust < 0.05 & direction == "Upregulated"])))
