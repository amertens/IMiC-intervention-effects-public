# =============================================================================
# build_all_manuscript_figures.R
#
# Regenerates the printed figures (Figs 1-6 and Figs S1-S5) in figure order. Each
# figure key lists its steps in run order: upstream data scripts first, then the
# panel scripts, then the script that writes figures/figureN.* (or figureSN_*.*).
# Keys whose steps read restricted or on-request inputs (see each script's header
# and the repository README) fail at that step when those inputs are absent.
#
# Run from the repo root:
#   Rscript figure-scripts/manuscript_figures/build_all_manuscript_figures.R          # every figure
#   Rscript figure-scripts/manuscript_figures/build_all_manuscript_figures.R 1 3 S2   # only these keys
#
# Each step runs in a fresh Rscript process (several scripts call rm(list = ls())
# or setwd(), so they must not share a session). A failed step is reported and the
# run continues; a summary is printed at the end.
# =============================================================================
root   <- here::here()
Rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
sel     <- commandArgs(trailingOnly = TRUE)          # optional: figure keys to build
# the generators write into these two folders, which are not part of the repository
for (d in c("figures", "figure-data")) dir.create(file.path(root, d), showWarnings = FALSE)

# Figure key -> ordered steps (paths relative to the repo root).
FIGS <- list(
  "1"  = list(desc = "ML classifier CV-AUC (A) + PCA effects (B) -> figures/figure1.{png,pdf,eps}",
              steps = c("src/3 visualizations/5-SL_VIM_plots.R",              # -> figure-data/SL_vim_plot_data.RDS
                        "figure-scripts/manuscript_figures/fig1-ml-vim-classifier.R")),
  "2"  = list(desc = "Primary-outcome forest -> figures/figure2.{png,pdf,eps}",
              steps = c("figure-scripts/manuscript_figures/fig2-primary-forest.R")),
  # Panel B is a local MetaboAnalystR KEGG pathway run (no web tool); the composite
  # embeds the Panel B PNG, so fig3B-pathway.R runs first.
  "3"  = list(desc = "Primary volcanoes (A) + KEGG pathway impact (B) -> figures/figure3.{png,pdf,eps}",
              steps = c("src/metaboanalyst/run-primary-pathway-local.R",           # -> primary_pathway_all_cells.csv (Panel B, Table S2)
                        "figure-scripts/manuscript_figures/fig3B-pathway.R",       # -> figures/figure3_panelB_msea.png
                        "figure-scripts/manuscript_figures/fig3-primary-volcano-composite.R")),
  "4"  = list(desc = "Milk nutrient distributions vs MILQ references -> figures/figure4.jpeg",
              steps = c("figure-scripts/manuscript_figures/fig4-milq-boxplots.R")),
  # Panel C needs the Biocrates structure file (see fig5C-triglyceride.R); the
  # composite embeds the Panel C PNG, so it runs after fig5C.
  "5"  = list(desc = "Tertiary volcanoes (A) + MSEA (B) + TG fatty acids (C) -> figures/figure5.{png,pdf,eps}",
              steps = c("src/metaboanalyst/run-tertiary-msea-dual.R",               # -> tertiary_msea_dual.csv (Panel B, Table S3)
                        "figure-scripts/manuscript_figures/fig5C-triglyceride.R",    # Panel C (+ Table S4 CSV)
                        "figure-scripts/manuscript_figures/fig5-tertiary-composite.R")),
  # Panel B reads the shipped Table S6 file (src/metaboanalyst/run-milk-mummichog.R,
  # not rerun here). Panel D reads the cross-compartment rows and pathway hits written
  # by run-compartment-pathway.R, which needs the query list built just before it.
  "6"  = list(desc = "Untargeted MSEA (A) + mummichog (B) + proteome GO (C) + cross-compartment (D) -> figures/figure6.{png,pdf}",
              steps = c("src/2 analysis/47-annotate-milk-features.R",             # -> milk_nominal_putative_annotation.csv (Panel A foreground)
                        "src/metaboanalyst/run-untargeted-msea.R",                # -> untargeted_msea_combined.csv (Panel A, Table S5)
                        "src/2 analysis/54-supplement-detection-status.R",        # -> supplement_status_fdr_features.csv (Panel D supplement labels)
                        "src/metaboanalyst/build-compartment-query-list.R",       # -> pathway_compound_list.csv (Panel D / Table S11 query list)
                        "src/metaboanalyst/run-compartment-pathway.R",            # -> linked_crosscompartment.csv + metabolite_pathways.csv (Panel D)
                        "src/2 analysis/55-proteomics-go-uniprot.R",              # -> results/proteomics_go_uniprot.csv (Panel C)
                        "figure-scripts/manuscript_figures/fig6A-untargeted-msea.R",
                        "figure-scripts/manuscript_figures/fig6B-mummichog.R",
                        "figure-scripts/manuscript_figures/fig6C-proteomics.R",
                        "figure-scripts/manuscript_figures/fig6D-crosscompartment.R",
                        "figure-scripts/manuscript_figures/fig6-composite.R")),  # A-D -> figure6.png
  "S1" = list(desc = "Child growth outcomes -> figures/figureS1_growth_outcomes.{png,pdf,eps}",
              steps = c("figure-scripts/manuscript_figures/figS1-growth-outcomes.R")),
  "S2" = list(desc = "Milk nutrient deficiency relative risks -> figures/figureS2_milq_deficiency.png",
              steps = c("figure-scripts/manuscript_figures/figS2-milq-deficiency.R")),
  "S3" = list(desc = "Secondary outcomes (HMO + bioactive) forest -> figures/figureS3_hmo_bioactives.{png,pdf,eps}",
              steps = c("figure-scripts/manuscript_figures/figS3-hmo-bioactives.R")),
  "S4" = list(desc = "Triglyceride treatment-specific means -> figures/figureS4_triglyceride_means.{png,pdf,eps}",
              steps = c("figure-scripts/manuscript_figures/figS4-triglyceride-means.R")),
  "S5" = list(desc = "Milk microbiome alpha diversity -> figures/figureS5_microbiome_diversity.{png,pdf,eps}",
              steps = c("figure-scripts/manuscript_figures/figS5-microbiome-diversity.R"))
)

keys <- if (length(sel)) intersect(sel, names(FIGS)) else names(FIGS)
results <- list()
for (k in keys) {
  cat(sprintf("\n########## FIGURE %s  %s\n", k, FIGS[[k]]$desc))
  for (s in FIGS[[k]]$steps) {
    cat(sprintf("  -> %s ... ", s))
    st <- system2(Rscript, shQuote(file.path(root, s)), stdout = FALSE, stderr = FALSE)
    cat(if (st == 0) "OK\n" else sprintf("FAILED (exit %s)\n", st))
    results[[paste(k, s)]] <- st
  }
}
cat("\n=========================== SUMMARY ===========================\n")
for (nm in names(results))
  cat(sprintf("  [%s] %s\n", if (results[[nm]] == 0) "ok  " else "FAIL", nm))
failed <- names(results)[unlist(results) != 0]
if (length(failed)) cat(sprintf("\n%d step(s) failed; see above.\n", length(failed))) else
  cat("\nAll figures rebuilt.\n")
