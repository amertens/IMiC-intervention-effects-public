# =============================================================================
# 53-proteomics-go-ora-rerun.R  INTERIM
# Reproduce the manuscript's milk-proteome GO Biological Process over-representation
# (clusterProfiler::enrichGO) on THIS MORNING's re-run proteomics file, so Fig 6C
# can be self-consistent (fresh count 7 + fresh GO terms) rather than old narrative
# + new count. Mirrors Methods: universe = all measured UniProt; foreground =
# per (study x visit x direction) features with unadjusted p<0.05, contrast BEP;
# enrichGO(BP), BH, pvalueCutoff=0.05, qvalueCutoff=0.20.
# INTERIM: proteomics is being re-run again by Trenton; treat as provisional.
# =============================================================================
suppressMessages({library(data.table); library(biotmle); library(clusterProfiler); library(org.Hs.eg.db)})
setwd(here::here())

d <- as.data.table(readRDS("results/adjusted_combined_arms_intervention_effects_proteomics_results_clean.RDS"))
d <- d[measure == "ATE" & contrast == "BEP"]

uni2ent <- suppressWarnings(bitr(unique(d$biomarker), "UNIPROT", "ENTREZID", org.Hs.eg.db))
universe <- unique(uni2ent$ENTREZID)
ent <- setNames(uni2ent$ENTREZID, uni2ent$UNIPROT)

cells <- unique(d[, .(study, visit)])
res <- list()
for (i in seq_len(nrow(cells))) {
  for (dir in c("up","down")) {
    sub <- d[study == cells$study[i] & visit == cells$visit[i] & pval < 0.05 &
             (if (dir=="up") est>0 else est<0)]
    fg <- unique(na.omit(ent[sub$biomarker]))
    cat(sprintf("%s %s %-4s : %d p<.05 features -> %d Entrez\n",
                cells$study[i], cells$visit[i], dir, nrow(sub), length(fg)))
    if (length(fg) < 3) next
    eg <- tryCatch(enrichGO(gene=fg, universe=universe, OrgDb=org.Hs.eg.db,
                    keyType="ENTREZID", ont="BP", pAdjustMethod="BH",
                    pvalueCutoff=0.05, qvalueCutoff=0.20, readable=TRUE),
                   error=function(e) NULL)
    if (is.null(eg) || !nrow(as.data.frame(eg))) next
    t <- as.data.table(as.data.frame(eg))
    t[, `:=`(study=cells$study[i], visit=cells$visit[i], direction=dir)]
    res[[paste(i,dir)]] <- t
  }
}
if (length(res)) {
  out <- rbindlist(res, fill=TRUE)
  fwrite(out, "results/proteomics_go_ora_rerun.csv")
  cat("\n==== INTERIM GO BP terms (q<0.20) on fresh proteomics ====\n")
  for (k in names(res)) {
    t <- res[[k]]
    cat(sprintf("\n-- %s %s / %s (%d terms) --\n", t$study[1], t$visit[1], t$direction[1], nrow(t)))
    print(head(t[order(pvalue), .(Description, GeneRatio, p=signif(pvalue,2), q=signif(qvalue,2))], 6), row.names=FALSE)
  }
  cat("\nSaved results/proteomics_go_ora_rerun.csv\n")
} else cat("\nNo GO terms passed q<0.20 in any cell/direction on the fresh data.\n")
