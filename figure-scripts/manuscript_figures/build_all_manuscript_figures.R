# =============================================================================
# build_all_manuscript_figures.R
#
# SINGLE ENTRY POINT to regenerate every IMiC manuscript figure, in figure-number
# order, with the full data -> panel -> composite chain made explicit. This driver
# is the SOURCE OF TRUTH for "which script builds which figure". As of 2026-08-13 the
# generator filenames self-document the figure they build: fig1-*, fig2-*,
# fig3-*/fig3B-*, fig4-*, fig5-*/fig5C-*, fig6A-*..fig6D-*/fig6-composite, and the
# supplements figS3-*/figS4-*/figS8-*/figS9-*. (The output PNGs keep their legacy
# names, e.g. fig3B-pathway.R still writes figures/figure4_panelB_msea.png.)
#
# Run from the repo root:
#   Rscript figure-scripts/manuscript_figures/build_all_manuscript_figures.R
#   Rscript figure-scripts/manuscript_figures/build_all_manuscript_figures.R 1 3 6   # only Figs 1, 3, 6
#
# Each step runs in a FRESH Rscript subprocess (scripts call rm(list=ls())/setwd(),
# so they must not share a session). Failures are reported but do not stop the run.
# =============================================================================
root   <- here::here()
Rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
sel     <- commandArgs(trailingOnly = TRUE)          # optional: figure keys to build

# --- the figure map: label -> ordered list of scripts (data-gen -> panels -> fig)
# Paths are relative to the repo root. A figure that needs upstream data or panel
# PNGs lists those steps FIRST, then the script that emits figures/figureN.*.
FIGS <- list(
  "1"  = list(desc = "ML variable-importance / CV-AUC classifier -> figures/figure1.png",
              steps = c("src/3 visualizations/5-SL_VIM_plots.R",              # -> figure-data/SL_vim_plot_data.RDS (incl. ELICIT Az negative control)
                        "figure-scripts/manuscript_figures/fig1-ml-vim-classifier.R")),
  "2"  = list(desc = "Primary-outcome forest -> figures/figure2.png",
              steps = c("figure-scripts/manuscript_figures/fig2-primary-forest.R")),
  "3"  = list(desc = "Primary volcano (A) + KEGG pathway-impact (B) -> figures/figure3.png (combined MAIN) + figures/figure3_stratified_supplement.png (arm-stratified SUPPLEMENT)",
              # Panel B = KEGG pathway-impact plot. WEB-INDEPENDENT: built from a LOCAL
              # MetaboAnalystR KEGG run (run-primary-pathway-local.R) that reproduces the
              # submitted impact axis exactly; fig3B-pathway.R writes
              # figures/figure4_panelB_msea.png (legacy output name), which
              # fig3-primary-volcano-composite.R embeds. That composite emits BOTH the
              # combined-arm main figure and the arm-stratified supplement; the stratified
              # Panel B needs the scripted per-arm pathway run first.
              steps = c("src/metaboanalyst/run-primary-pathway-local.R",                # -> primary_pathway_local (LOCAL KEGG combined Panel B; web-independent)
                        "src/metaboanalyst/run-primary-pathway-stratified.R",           # -> primary_stratified_pathway (stratified Panel B data)
                        "figure-scripts/manuscript_figures/fig3B-pathway.R",    # -> figures/figure4_panelB_msea.png (combined B) + figure3_panelB_stratified.png (stratified B)
                        "figure-scripts/manuscript_figures/fig3-primary-volcano-composite.R")),
  "4"  = list(desc = "HM nutrient distributions vs MILQ -> figures/figure4.jpeg",
              steps = c("figure-scripts/manuscript_figures/fig4-milq-boxplots.R")),
  "5"  = list(desc = "Tertiary targeted metabolites -> figures/figure5.png",
              # NB: the real generator is fig5-tertiary-composite.R (+ the fig5C-* Panel C
              # builder); the old fig5.R is a STALE purl with blank Panel B/C placeholders.
              # run-tertiary-msea-dual.R added 2026-08-25: it's Panel B's actual data
              # source (tertiary_msea_dual.csv) and was previously missing from this chain,
              # so a clean-tree rebuild would have silently plotted a stale/absent Panel B
              # instead of regenerating it -- see Manuscript/CODE_AUDIT_2026-08-25.md #4.
              steps = c("src/metaboanalyst/run-tertiary-msea-dual.R",             # -> tertiary_msea_dual.csv (Panel B)
                        "figure-scripts/manuscript_figures/fig5C-triglyceride.R",  # Panel C
                        "figure-scripts/manuscript_figures/fig5-tertiary-composite.R")),                        # A + B + C -> figure5.png
  "6"  = list(desc = "Exploratory MSEA (A) + Mummichog (B) + proteome GO (C) + cross-compartment pathways (D) -> figures/figure6.png",
              # All four panels are WEB-INDEPENDENT. 6A: local reference-background ORA on
              # the NOMINAL (P<0.05) milk foreground (script 47 -> run-untargeted-msea.R).
              # 6D: Trenton's pathway-annotated cross-compartment figure, with the web
              # pathway_results.csv + hand-coded tribble replaced by a LOCAL KEGG pathway
              # run (script 54 supplement status -> run-compartment-pathway.R).
              steps = c("src/2 analysis/47-annotate-milk-features.R",         # -> milk_nominal_putative_annotation.csv (6A foreground)
                        "src/metaboanalyst/run-untargeted-msea.R",            # -> untargeted_msea_combined.csv (6A local ORA)
                        "src/2 analysis/54-supplement-detection-status.R",    # -> supplement_status_fdr_features.csv (6D supplement labels)
                        "src/metaboanalyst/run-compartment-pathway.R",        # -> compartment_tracking local KEGG pathway (6D)
                        "src/2 analysis/55-proteomics-go-uniprot.R",          # -> results/proteomics_go_uniprot*.csv (Panel C data)
                        "figure-scripts/manuscript_figures/fig6A-untargeted-msea.R",        # Panel A
                        "figure-scripts/manuscript_figures/fig6B-mummichog.R",       # Panel B
                        "figure-scripts/manuscript_figures/fig6C-proteomics.R",      # Panel C
                        "figure-scripts/manuscript_figures/fig6D-crosscompartment.R",# Panel D (pathway-annotated)
                        "figure-scripts/manuscript_figures/fig6-composite.R")),            # A+B+C+D -> figure6.png
  # RETIRED 2026-08-26: cut from the printed supplement and manuscript entirely (not
  # just renumbered) per author decision. Static image + its unused fig1.R rebuild
  # archived to archive/figures/; fig1.R/fig1.Rmd archived to archive/figure-scripts/.
  # No other figure/table depends on this key.
  "X-S1-study-design" = list(desc = "[RETIRED, was Fig S1] Study design / sampling schematic -- deleted entirely, not renumbered",
              steps = character(0)),
  # renumbered 2026-08-26: was Fig S2, now Fig S1 (S1-study-design deleted, not just renumbered)
  "S1" = list(desc = "Child growth outcomes -> figures/figureS2_growth_outcomes.png (reconstructed 2026-08-20; replaces an extracted legacy image; filename keeps its legacy S2 name)",
              steps = c("figure-scripts/manuscript_figures/figS2-growth-outcomes.R")),
  # renumbered 2026-08-26: was Fig S3, now Fig S2
  "S2" = list(desc = "Micronutrient-deficiency RR forest -> figures/fig-milq-deficiency-reduction-forest-plot.png",
              steps = c("figure-scripts/manuscript_figures/figS3-milq-deficiency.R")),
  # renumbered 2026-08-26: was Fig S4, now Fig S3
  "S3" = list(desc = "Secondary (HMO + bioactive) forest -> figures/figureS4.png (filename keeps its legacy S4 name)",
              steps = c("figure-scripts/manuscript_figures/supporting/figS4-secondary-forest.R")),
  # renumbered 2026-08-26: was Fig S5, now Fig S4
  "S4" = list(desc = "Triglyceride treatment-specific means -> figures/figureS5_triglyceride_means.png (reconstructed 2026-08-20; filename keeps its legacy S5 name)",
              steps = c("figure-scripts/manuscript_figures/figS5-triglyceride-means.R")),
  # renumbered 2026-08-26: was Fig S6, now Fig S5
  "S5" = list(desc = "Milk microbiome alpha diversity -> figures/figureS6_microbiome_diversity.png (reconstructed 2026-08-20; filename keeps its legacy S6 name)",
              steps = c("figure-scripts/manuscript_figures/figS6-microbiome-diversity.R")),
  # RETIRED 2026-08-26: cut from the printed supplement and manuscript entirely (not
  # just renumbered) per author decision. Script and PNG archived (archive-never-delete)
  # to archive/src/2 analysis/ and archive/figures/cross_compartment/.
  "X-S7-blood-transfer" = list(desc = "[RETIRED, was Fig S7] Blood cross-compartment transfer -- deleted entirely, not renumbered",
              steps = character(0)),
  # RETIRED 2026-08-16: cut from the printed supplement because it re-presented
  # Table S9. Script and PNG are kept (archive-never-delete) and the key is retained
  # so the artifact can still be rebuilt; it is no longer a numbered exhibit.
  "X-blood-class" = list(desc = "[RETIRED, was Fig S8] Blood chemical-class enrichment -> figures/figureS8_blood_class_enrichment.png; the result is printed as Table S9",
              steps = c("figure-scripts/manuscript_figures/figS8-blood-class-enrichment.R")),
  # renumbered: was Fig S9, then Fig S8 (2026-08-25), now Fig S6 (2026-08-26, after S7's deletion)
  "S6" = list(desc = "Direction-split ORA -> figures/figureS_{primary,tertiary}_ora_by_direction.png",
              steps = c("figure-scripts/manuscript_figures/figS9-ora-by-direction.R")),
  # S10/S11 previously had NO generator: both were orphan PNGs from Trenton
  # (trenton_xcompartment_*.png) whose input CSV was regenerated by a script while the
  # images were not, so they drifted apart silently. Both are now scripted from the
  # in-repo artifacts; see Manuscript/PROVENANCE_AUDIT.md.
  # RETIRED 2026-08-16: cut from the printed supplement as a near-duplicate of main
  # Fig 6D. The script still runs and is STILL REQUIRED, because it produces the
  # Table S11 compound list under Trenton's grouping rule; only its PNG is unused.
  "X-compartment-tracking" = list(desc = "[RETIRED as a figure, was Fig S10] Cross-compartment tracking -> figures/cross_compartment/figureS10_compartment_tracking.png; STILL NEEDED for results/compartment_tracking/trenton_pathway_compound_list.csv (Table S11)",
               steps = c("src/2 analysis/54-supplement-detection-status.R",    # -> supplement_status_fdr_features.csv
                         "figure-scripts/manuscript_figures/figS10-compartment-tracking.R",
                         "src/metaboanalyst/run-compartment-pathway-trenton.R")),
  # renumbered: was Fig S11, then Fig S9 (2026-08-25), now Fig S7 (2026-08-26, after S7's deletion)
  "S7" = list(desc = "Cross-compartment volcano grid -> figures/cross_compartment/figureS11_crosscompartment_volcanoes.png (filename keeps its legacy S11 name)",
               steps = c("figure-scripts/manuscript_figures/figS11-crosscompartment-volcanoes.R")),
  # RETIRED 2026-08-26: cut from the printed supplement and manuscript entirely (not
  # just renumbered) per author decision. Script and PNG archived (archive-never-delete)
  # to archive/figure-scripts/manuscript_figures/ and archive/figures/.
  "X-S10-proteomics-individual" = list(desc = "[RETIRED, was Fig S10] Individual-protein effects, untargeted proteomics -- deleted entirely, not renumbered",
               steps = character(0))
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
