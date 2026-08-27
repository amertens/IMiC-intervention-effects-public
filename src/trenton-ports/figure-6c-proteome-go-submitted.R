# =============================================================================
# figure-6c-proteome-go-submitted.R
# Faithful R port of the SUBMITTED Fig 6C figure code, found in Trenton's older
# "imicPaperProteomics.Rmd" (2025-07-23) -- the figure generator that the
# 2026 "2. Scripts/Exploratory Outcomes (Proteomics).Rmd" (analysis-only) lacked.
#
# This reproduces the SUBMITTED proteome-GO panel exactly: gene-level enrichGO on
# UniProt->Entrez ids against a GLOBAL universe (ALL measured proteomics proteins,
# pooled), 4 cells, then the signed-fold volcano coloured by significance tier.
#   Universe   : pbl_background_vector_entrez = every proteomics UniProt in the
#                STRATIFIED-arms RDS, mapped to Entrez (his L66-80). GLOBAL, not per-cell.
#   Cells      : MISAME 1-2mo down/up; Mumpta(=Vital) 2mo down/up. sig==1, est sign,
#                contrast != control.
#   enrichGO   : ont="BP", BH, pvalueCutoff=0.05, qvalueCutoff=0.2 (his exact params).
#   Figure     : x = NegativeFoldEnrichment (down = -abs FE), y = -log10(pvalue),
#                size = Count, colour tier grey/orange/green by (p, q); label the
#                top-10 (by p) FDR-significant terms per sign. (his "★Untargeted
#                Proteomics Figure", saved for_andrew_figure_4b.png.)
#
# CODING ERROR FIXED (not preserved): in his Rmd the "Mumpta upregulated" cell
# (L315-323) filters `studytime == "Misame (1-2 mo.)"` -- a copy-paste bug, so the
# submitted Mumpta-up panel actually plotted MISAME data. Here it correctly uses
# "Vital (2 mo.)". This is the one clear error corrected.
#
# RELATION TO THE OTHER PROTEOMICS PORTS: this is the SUBMITTED (global-universe,
# gene-level) result -- it is why the submitted up-regulated anti-proteolysis /
# immunity terms survive FDR. The corrected methods (per-cell background, and the
# UniProt-native each-protein-once test) live in exploratory-proteomics-genelevel.R
# and exploratory-proteomics-uniprot.R and OVERTURN that up-regulated result. Keep
# this port only to reproduce/annotate the submitted figure, not as the analysis.
#
# Run from repo root:  Rscript src/trenton-ports/figure-6c-proteome-go-submitted.R
# =============================================================================
suppressMessages({
  library(dplyr); library(stringr); library(ggplot2); library(ggrepel)
  library(clusterProfiler); library(org.Hs.eg.db)
})
STRATIFIED_RDS <- "results/combined_intervention_effects_results_stratified_arms.RDS"
OUT <- "results/trenton-ports/figure_6c"
PNG <- "figures/trenton_ports/figure_6c_proteome_go_submitted.png"

strat <- readRDS(STRATIFIED_RDS)

# GLOBAL universe: all measured proteomics UniProt -> Entrez (his pbl_background).
pbl_entrez <- strat %>% filter(description == "Proteomics") %>% pull(label_f) %>% unique() %>%
  { suppressWarnings(bitr(., fromType = "UNIPROT", toType = "ENTREZID", OrgDb = org.Hs.eg.db)) } %>%
  pull(ENTREZID) %>% unique()

run_cell <- function(studytime_value, direction, study, timepoint) {
  fg <- strat %>% filter(studytime == studytime_value, contrast != "control",
                         description == "Proteomics", sig == 1, measure == "ATE",
                         if (direction == "upregulated") est > 0 else est < 0) %>%
    pull(label_f) %>% unique()
  ge <- suppressWarnings(bitr(fg, fromType = "UNIPROT", toType = "ENTREZID", OrgDb = org.Hs.eg.db)) %>%
    pull(ENTREZID) %>% unique()
  if (!length(ge)) return(NULL)
  e <- enrichGO(gene = ge, universe = pbl_entrez, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
  if (is.null(e) || !nrow(as.data.frame(e))) return(NULL)
  as_tibble(e@result) %>% mutate(study = study, timepoint = timepoint, direction = direction)
}

all_enrich <- bind_rows(
  run_cell("Misame (1-2 mo.)", "downregulated", "misame", "pn12"),
  run_cell("Misame (1-2 mo.)", "upregulated",   "misame", "pn12"),
  run_cell("Vital (2 mo.)",    "downregulated", "mumpta", "pn2"),
  run_cell("Vital (2 mo.)",    "upregulated",   "mumpta", "pn2")   # FIX: was "Misame (1-2 mo.)"
) %>%
  mutate(NegativeFoldEnrichment = if_else(direction == "downregulated", -abs(FoldEnrichment), FoldEnrichment),
         sig_color = case_when(pvalue > 0.05 ~ "grey",
                               pvalue < 0.05 & qvalue > 0.05 ~ "orange",
                               qvalue < 0.05 ~ "green"))

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
readr::write_csv(all_enrich, file.path(OUT, "figure_6c_proteome_go_all_cells.csv"))

labels_to_plot <- bind_rows(
  all_enrich %>% filter(qvalue < 0.05, NegativeFoldEnrichment < 0) %>% arrange(pvalue) %>% slice_head(n = 10),
  all_enrich %>% filter(qvalue < 0.05, NegativeFoldEnrichment > 0) %>% arrange(pvalue) %>% slice_head(n = 10))

p <- ggplot(all_enrich, aes(NegativeFoldEnrichment, -log10(pvalue), size = Count, alpha = Count, color = sig_color)) +
  geom_point() +
  geom_label_repel(data = labels_to_plot, aes(label = Description), color = "black", alpha = 1,
                   size = 2, max.overlaps = 10, show.legend = FALSE) +
  scale_color_manual(values = c(grey = "grey", orange = "orange", green = "limegreen"),
                     name = "Statistical Significance",
                     labels = c(green = "Q < 0.05", grey = "Not significant", orange = "P < 0.05, Q >= 0.05")) +
  scale_alpha_continuous(range = c(0.4, 1)) +
  labs(x = "Fold Enrichment (signed by direction)", y = expression(-Log[10]*"("*italic(P)*"-value)")) +
  theme_classic(base_size = 7) +
  theme(axis.title = element_text(size = 7), axis.text = element_text(size = 7),
        legend.title = element_text(size = 7), legend.text = element_text(size = 7))

dir.create(dirname(PNG), recursive = TRUE, showWarnings = FALSE)
ggsave(PNG, p, width = 210/2, height = 297/3, units = "mm", dpi = 600, device = ragg::agg_png)
cat(sprintf("wrote %s and %s | %d GO terms; up FDR-sig: %d, down FDR-sig: %d\n",
            PNG, OUT, nrow(all_enrich),
            sum(all_enrich$qvalue < 0.05 & all_enrich$NegativeFoldEnrichment > 0, na.rm = TRUE),
            sum(all_enrich$qvalue < 0.05 & all_enrich$NegativeFoldEnrichment < 0, na.rm = TRUE)))
