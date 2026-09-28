# =============================================================================
# 55-proteomics-go-uniprot.R
#
# Gene Ontology Biological-Process over-representation of the BEP-responsive milk
# proteins, run on UniProt identifiers; the enrichment behind Fig. 6C and Table S7
# (src/metaboanalyst/run-proteomics-go.R reads this output). Reproduces the
# UniProt-based proteomics GO analysis of the original R Markdown workflow from the
# shared milk results.
#
# Method: enricher() on a custom UniProt -> GO(BP) TERM2GENE map built from
# org.Hs.eg.db, so each measured protein is counted once in foreground and
# background (clusterProfiler::enrichGO expands one UniProt id into several Entrez
# genes and inflates the hypergeometric counts). Foreground = nominally significant
# proteins (sig == 1) split by effect direction; background = all measured proteins in
# that study x contrast cell; BH-adjusted p-values. MISAME-III uses the combined BEP
# contrast; Mumta-LW the two stratified contrasts (BEP+ExBf, BEP+ExBf+AZT). ELICIT has
# no proteomics. Fold enrichment is signed negative for down-regulated results.
#
# Inputs : results/combined_intervention_effects_results_combined_arms.RDS
#          results/combined_intervention_effects_results_stratified_arms.RDS
# Outputs: results/proteomics_go_uniprot.csv         (all GO-BP results, signed fold enrichment)
#          results/proteomics_go_uniprot_fdrsig.csv  (FDR < 0.05 subset)
# [needs restricted data] (feature-level result files, not shipped)
# =============================================================================
suppressMessages({
  library(clusterProfiler); library(org.Hs.eg.db); library(GO.db)
  library(AnnotationDbi); library(data.table)
})
root <- paste0(here::here(), "/")

combined   <- as.data.table(readRDS(paste0(root, "results/combined_intervention_effects_results_combined_arms.RDS")))
stratified <- as.data.table(readRDS(paste0(root, "results/combined_intervention_effects_results_stratified_arms.RDS")))

# -- all measured proteomics UniProt ids (fore+back universe) ------------------
prot_ids <- unique(rbind(combined, stratified, fill = TRUE)[
  description == "Proteomics" & measure == "ATE" & !is.na(label_f) & label_f != "", label_f])
cat(sprintf("measured proteomics UniProt ids: %d\n", length(prot_ids)))

# -- UniProt -> GO(BP) maps (each protein counted once) -----------------------
go_ann <- as.data.table(AnnotationDbi::select(org.Hs.eg.db, keys = prot_ids, keytype = "UNIPROT",
                                              columns = c("GOALL", "ONTOLOGYALL")))
go_ann <- unique(go_ann[ONTOLOGYALL == "BP" & !is.na(GOALL) & GOALL != "" & !is.na(UNIPROT) & UNIPROT != "",
                        .(term = GOALL, gene = UNIPROT)])
term2gene <- unique(go_ann[, .(term, gene)])
t2n <- as.data.table(AnnotationDbi::select(GO.db, keys = unique(term2gene$term), keytype = "GOID", columns = "TERM"))
term2name <- unique(t2n[!is.na(TERM) & TERM != "", .(term = GOID, name = TERM)])
cat(sprintf("UniProt->GO(BP): %d proteins mapped / %d measured\n",
            length(intersect(prot_ids, term2gene$gene)), length(prot_ids)))

# -- one enrichment cell ------------------------------------------------------
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

# -- cells: MISAME combined (BEP); Mumta stratified (two contrasts) -----------
cells <- rbindlist(list(
  run_cell(combined,   "Misame (1-2 mo.)", "BEP",          "Downregulated"),
  run_cell(combined,   "Misame (1-2 mo.)", "BEP",          "Upregulated"),
  run_cell(stratified, "Vital (2 mo.)",    "BEP+ExBf",     "Downregulated"),
  run_cell(stratified, "Vital (2 mo.)",    "BEP+ExBf",     "Upregulated"),
  run_cell(stratified, "Vital (2 mo.)",    "BEP+ExBf+AZT", "Downregulated"),
  run_cell(stratified, "Vital (2 mo.)",    "BEP+ExBf+AZT", "Upregulated")
), fill = TRUE)

# -- tidy: signed fold enrichment, readable labels ----------------------------
cells[, study := fifelse(grepl("Misame", studytime), "MISAME-III", "Mumta-LW")]
cells[, timepoint := fifelse(grepl("Misame", studytime), "1-2 mo.", "2 mo.")]
cells[, signed_fold_enrichment := fifelse(direction == "Upregulated", abs(FoldEnrichment), -abs(FoldEnrichment))]
cells[, description_sentence := paste0(toupper(substr(Description, 1, 1)), substr(Description, 2, nchar(Description)))]
setorder(cells, study, timepoint, contrast, direction, pvalue)

fwrite(cells, paste0(root, "results/proteomics_go_uniprot.csv"))
fwrite(cells[p.adjust < 0.05], paste0(root, "results/proteomics_go_uniprot_fdrsig.csv"))
cat(sprintf("\nwrote results/proteomics_go_uniprot.csv  (%d GO-BP terms, %d FDR<0.05)\n",
            nrow(cells), nrow(cells[p.adjust < 0.05])))
cat("\nFDR-significant per cell:\n")
print(cells[p.adjust < 0.05, .N, by = .(study, contrast, direction)])
