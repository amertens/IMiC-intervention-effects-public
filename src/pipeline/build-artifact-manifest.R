# =============================================================================
# build-artifact-manifest.R
#
# Emit results/ARTIFACT_MANIFEST.csv -- the machine-readable contract between this
# repo and anything that consumes its outputs (chiefly the online supplement's
# port_results.R, which today hardcodes a hand-maintained path list and needs a new
# entry every time an artifact lands somewhere new).
#
# The manifest answers, for every manuscript exhibit:
#   exhibit  - "Fig 6D", "Table S11", ... (the label the manuscript uses)
#   kind     - figure | table | data
#   path     - repo-relative artifact path
#   generator- the ONE script that writes it ("" = static/no generator, which is
#              itself a finding: see Manuscript/PROVENANCE_AUDIT.md)
#   depends  - semicolon-separated upstream artifacts it is derived from; this is
#              what makes the staleness check possible (a figure must never be
#              older than the data it plots -- the Fig. S10 failure mode)
#   exists / bytes / mtime - stamped live at build time
#
# The REGISTRY below is the single source of truth. Adding an exhibit means adding
# one row here -- not editing a downstream repo.
#
# Run from repo root: Rscript src/pipeline/build-artifact-manifest.R
# =============================================================================
suppressMessages({ library(data.table) })
root <- paste0(here::here(), "/")
OUT  <- paste0(root, "results/ARTIFACT_MANIFEST.csv")

FS  <- "figure-scripts/manuscript_figures/"
A2  <- "src/2 analysis/"
MA  <- "src/metaboanalyst/"

r <- function(exhibit, kind, path, generator = "", depends = "", notes = "")
  data.table(exhibit, kind, path, generator, depends, notes)

REGISTRY <- rbindlist(list(
  # ---- main figures -------------------------------------------------------
  r("Fig 1",  "figure", "figures/figure1.png", paste0(FS, "fig1-ml-vim-classifier.R"),
    "figure-data/SL_vim_plot_data.RDS"),
  r("Fig 2",  "figure", "figures/figure2.png", paste0(FS, "fig2-primary-forest.R"),
    "figure-data/primary_forest_plots.RDS"),
  r("Fig 3",  "figure", "figures/figure3.png", paste0(FS, "fig3-primary-volcano-composite.R"),
    "figures/figure4_panelB_msea.png"),
  r("Fig 3B", "figure", "figures/figure4_panelB_msea.png", paste0(FS, "fig3B-pathway.R"),
    "results/metaboanalyst/primary_pathway_local/primary_pathway_all_cells.csv",
    "legacy output name: fig3B writes figure4_panelB_msea.png"),
  r("Fig 4",  "figure", "figures/figure4.jpeg", paste0(FS, "fig4-milq-boxplots.R")),
  r("Fig 5",  "figure", "figures/figure5.png", paste0(FS, "fig5-tertiary-composite.R"),
    "figures/figure4_panelC_tg_composition.png"),
  r("Fig 5C", "figure", "figures/figure4_panelC_tg_composition.png", paste0(FS, "fig5C-triglyceride.R"),
    "results/metaboanalyst/triglyceride_fa/triglyceride_fa_composition_combined.csv",
    "legacy output name: fig5C writes figure4_panelC_*.png"),
  r("Fig 6",  "figure", "figures/figure6.png", paste0(FS, "fig6-composite.R"),
    paste("figures/figure6_panelA_untargeted_msea.png",
          "figures/figure6_panelB_mummichog.png",
          "figures/figure6_panelC_proteomics_go.png",
          "figures/figure6_panelD_crosscompartment.png", sep = ";")),
  r("Fig 6A", "figure", "figures/figure6_panelA_untargeted_msea.png", paste0(FS, "fig6A-untargeted-msea.R"),
    "results/metaboanalyst/untargeted_msea/untargeted_msea_combined.csv"),
  r("Fig 6B", "figure", "figures/figure6_panelB_mummichog.png", paste0(FS, "fig6B-mummichog.R"),
    "results/metaboanalyst/mummichog_s5/milk_mummichog_tableS5.csv"),
  r("Fig 6C", "figure", "figures/figure6_panelC_proteomics_go.png", paste0(FS, "fig6C-proteomics.R"),
    "results/proteomics_go_uniprot.csv"),
  r("Fig 6D", "figure", "figures/figure6_panelD_crosscompartment.png", paste0(FS, "fig6D-crosscompartment.R"),
    paste("results/compartment_tracking/trenton_linked_crosscompartment.csv",
          "results/compartment_tracking/metabolite_pathways_trenton.csv", sep = ";"),
    "Trenton's name-based grouping (src/trenton-ports/compartment-tracking.R::ct_trenton_track()), same rule as Table S11 and Fig S10 -- see Manuscript/PROVENANCE_AUDIT.md 6.4 (resolved)"),

  # ---- supplementary figures ---------------------------------------------
  # RETIRED 2026-08-26: deleted entirely (not renumbered) per author decision. Static
  # image + fig1.R's unused rebuild archived to archive/figures/; fig1.R/fig1.Rmd
  # archived to archive/figure-scripts/manuscript_figures/.
  r("(retired) Fig S1", "figure", "archive/figures/figureS1.png", "", "",
    "DELETED 2026-08-26, not renumbered. Was STATIC: the exact submitted schematic; deliberately not rebuilt"),
  # renumbered 2026-08-26: was Fig S2, now Fig S1 (Fig S1-study-design deleted above)
  r("Fig S1", "figure", "figures/figureS2_growth_outcomes.png",
    paste0(FS, "figS2-growth-outcomes.R"), "results/growth_intervention_effects_results.RDS",
    "Filename keeps its legacy S2 name. PROMOTED 2026-08-25 (was STATIC, the submitted image at Manuscript/qmd/extracted/media_supplement/media/image1.png): supplement_v2.qmd already embedded this reconstructed figure since 2026-08-20 -- this row previously disagreed with that and with build_all_manuscript_figures.R, which already treated it as live; corrected to match. Static original vs. reconstruction is compared side-by-side in Manuscript/figure_comparison/build_figure_comparison.R / figures_old_vs_new.html -- see Manuscript/CODE_AUDIT_2026-08-25.md #5"),
  # renumbered 2026-08-26: was Fig S3, now Fig S2
  r("Fig S2", "figure", "figures/fig-milq-deficiency-reduction-forest-plot.png",
    paste0(FS, "figS3-milq-deficiency.R")),
  # renumbered 2026-08-26: was Fig S4, now Fig S3
  r("Fig S3", "figure", "figures/figureS4.png", paste0(FS, "supporting/figS4-secondary-forest.R"),
    "", "Filename keeps its legacy S4 name."),
  # renumbered 2026-08-26: was Fig S5, now Fig S4
  r("Fig S4", "figure", "figures/figureS5_triglyceride_means.png",
    paste0(FS, "figS5-triglyceride-means.R"),
    paste0("data/merged_analysis_datasets.RDS;",
           "results/adjusted_intervention_effects_results_clean.RDS;",
           "results/fat_adjusted_metabolomics_intervention_effects_results.RDS"),
    "Filename keeps its legacy S5 name. PROMOTED 2026-08-25 (was STATIC, the submitted image at Manuscript/qmd/extracted/media_supplement/media/image5.png). Same situation as Fig S1 above -- comparison preserved in Manuscript/figure_comparison/"),
  # renumbered 2026-08-26: was Fig S6, now Fig S5
  r("Fig S5", "figure", "figures/figureS6_microbiome_diversity.png",
    paste0(FS, "figS6-microbiome-diversity.R"), "results/microbiome_diversity_intervention_effects_results.RDS",
    "Filename keeps its legacy S6 name. PROMOTED 2026-08-25 (was STATIC, the submitted image at Manuscript/qmd/extracted/media_supplement/media/image6.png). Same situation as Fig S1 above -- comparison preserved in Manuscript/figure_comparison/"),
  # RETIRED 2026-08-26: deleted entirely (not renumbered) per author decision. Script +
  # PNG archived to archive/src/2 analysis/ and archive/figures/cross_compartment/.
  r("(retired) Fig S7", "figure", "archive/figures/cross_compartment/fig_blood_transfer_supplement.png",
    paste0(A2, "45-figure-blood-transfer-supplement.R"),
    "results/blood_compartment_all_FDRsig_ATE.csv",
    "DELETED 2026-08-26, not renumbered."),
  r("(retired) Fig S8", "figure", "figures/figureS8_blood_class_enrichment.png",
    paste0(FS, "figS8-blood-class-enrichment.R"),
    "results/blood_chemical_class_enrichment_directional.csv"),
  # renumbered: was Fig S8a/S9a, now Fig S6a (2026-08-26, after Fig S7's deletion)
  r("Fig S6a", "figure", "figures/figureS_primary_ora_by_direction.png",
    paste0(FS, "figS9-ora-by-direction.R"),
    "results/metaboanalyst/primary_combined/primary_ora_upregulated.csv"),
  r("Fig S6b", "figure", "figures/figureS_tertiary_ora_by_direction.png",
    paste0(FS, "figS9-ora-by-direction.R"),
    "results/metaboanalyst/tertiary_combined/tertiary_ora_upregulated.csv"),
  r("(retired) Fig S10", "figure", "figures/cross_compartment/figureS10_compartment_tracking.png",
    paste0(FS, "figS10-compartment-tracking.R"),
    paste("results/supplement_status_fdr_features.csv",
          "results/compartment_tracking/metabolite_pathways_trenton.csv", sep = ";"),
    "Not cited in manuscript_v4.qmd or supplement_v2.qmd -- fully subsumed by Table S8 (feature list) + Table S11 (pathways) + figS10_tracking_counts.csv (Results-text counts). Still built (two-pass, see script header) because it emits the Table S11 query list; kept out of the manuscript's numbered figure set on purpose. Still shown as an online-only bonus figure in IMiC-intervention-effects-supplement's 10_cross_compartment_blood.Rmd (§10.5)"),
  # renumbered: was Fig S11, then Fig S9 (2026-08-25), now Fig S7 (2026-08-26, after old Fig S7's deletion)
  r("Fig S7", "figure", "figures/cross_compartment/figureS11_crosscompartment_volcanoes.png",
    paste0(FS, "figS11-crosscompartment-volcanoes.R"),
    paste("results/blood_compartment_adjusted_combined_arms_intervention_effects_results_clean.RDS",
          "results/adjusted_combined_arms_intervention_effects_untargeted_results_clean_ATE.RDS",
          "results/supplement_status_fdr_features.csv", sep = ";"),
    "Filename keeps its legacy S11 name."),
  # RETIRED 2026-08-26: deleted entirely (not renumbered) per author decision. Script +
  # PNG archived to archive/figure-scripts/manuscript_figures/ and archive/figures/.
  r("(retired) Fig S10 individual-protein", "figure", "archive/figures/figureS10_proteomics_individual.png",
    paste0(FS, "figS10-proteomics-individual.R"),
    "results/adjusted_combined_arms_intervention_effects_proteomics_results_clean.RDS",
    "DELETED 2026-08-26, not renumbered. Was the individual-protein plot (distinct from the retired compartment-tracking exhibit above, which also used to carry the 'Fig S10' slot at an earlier point)."),

  # ---- supplementary tables ----------------------------------------------
  r("Table S1", "table", "results/metaboanalyst/primary_pathway_local/primary_pathway_all_cells.csv",
    paste0(MA, "run-primary-pathway-local.R"),
    "results/combined_intervention_effects_results_combined_arms.RDS",
    "CORRECTED 2026-08-26: previously pointed at primary_combined/primary_combined_supplementary_table.csv (run-primary.R), an ORA-module file (SMPDB names, no Impact column, ~438 rows) that is not the printed Table S1. The printed Table S1 is the KEGG pathway-topology (pathora) run -- the SAME run that produces Fig. 3B -- and was a frozen hand-pasted transcription of Trenton's submitted MetaboAnalyst web exports. Per author decision it was regenerated 2026-08-26 from the current pipeline (run-primary-pathway-local.R -> primary_pathway_all_cells.csv, with metabolite hit membership in the sibling primary_pathway_hits_all_cells.csv); a few borderline pathways shift across FDR 0.05 vs the submitted version (current KEGG metpa library). supplement_v2.qmd Table S1 is now a script-regenerable pipe table, not frozen text."),
  r("Table S2", "table", "results/metaboanalyst/tertiary_msea/tertiary_msea_dual.csv",
    paste0(MA, "run-tertiary-msea-dual.R"), "",
    paste0("EXTERNALIZED 2026-08-26: the 463-row table was removed from supplement_v2.qmd (Editor SS3.14, same treatment as Table S7) and replaced with a caption + pointer to results/tables/table_s2_tertiary_msea_full.csv (written by src/pipeline/rebuild-submission-tables.R) for the Zenodo Online Appendix. ",
           "CORRECTED 2026-08-26: this row previously pointed at tertiary_combined_supplementary_table.csv (run-secondary-tertiary.R), which is neither the printed table nor Fig 5B's source and was never actually verified. The printed Table S2 was a stale, undocumented web export (Trenton's original MetaboAnalyst downloads via run-tertiary-msea-from-tables.R) that disagreed with Fig 5B's own local re-implementation with no note explaining why. Per author decision, Table S2 was regenerated from tertiary_msea_dual.csv (the same source Fig 5B uses) so table and figure agree -- see src/pipeline/rebuild-submission-tables.R for the exact recipe (463 rows, unfiltered, sorted by P, matching Table S1's convention).")),
  r("Table S3", "table",
    paste("results/metaboanalyst/triglyceride_fa/triglyceride_fa_composition_combined.csv",
          "results/metaboanalyst/triglyceride_fa/triglyceride_fa_composition_stratified.csv", sep = ";"),
    paste0(FS, "fig5C-triglyceride.R"), "",
    "CORRECTED 2026-08-26: the printed table mixes combined-arm ELICIT/Mumta-LW rows with stratified-arm MISAME-III rows, but only the combined-arm half had ever been (re)built from this script -- the MISAME-III rows had drifted from current data (up to 272x P-value drift, 14 rows absent from a current re-run) and 16 rows had alpha-/gamma-linolenic acid names swapped (this script's own fatty_acid_name output was and is correct; only the printed table had it backwards). fig5C-triglyceride.R now supports both arm framings via IMIC_TG_ARM_FRAMING=combined|stratified (env var, default combined) and writes a companion CSV named accordingly; Table S3 is the P<0.05 union of both, sorted by P -- see src/pipeline/rebuild-submission-tables.R for the exact recipe (44 rows)."),
  r("Table S4", "table", "results/metaboanalyst/untargeted_msea/untargeted_msea_combined_perContrastFDRsig_tableS4.csv",
    paste0(MA, "run-untargeted-msea.R"), "",
    "VERIFIED 2026-08-26: re-running this script reproduces the printed Table S4 row-for-row (330-row combined output, 68-row FDR-significant artifact, all.equal() TRUE) and Fig 6A re-renders byte-identical. The old note here claiming Table S4 was a web-export-only, non-byte-matching table is superseded -- this is now the local scripted source of truth, same as the other supplementary tables."),
  r("Table S5", "table", "results/metaboanalyst/mummichog_s5/milk_mummichog_tableS5.csv",
    paste0(MA, "run-milk-mummichog-s5.R"),
    paste("results/adjusted_intervention_effects_res_untargeted_metabolomics_clean_ATE.RDS",
          "results/adjusted_combined_arms_intervention_effects_untargeted_results_clean_ATE.RDS", sep = ";"),
    "EXTERNALIZED 2026-08-26: the 505-row printed table (raw P < 0.05 of this file) was removed from supplement_v2.qmd (Editor SS3.14, same treatment as Table S7) and replaced with a caption + pointer to results/tables/table_s5_mummichog_full.csv (written by src/pipeline/rebuild-submission-tables.R) for the Zenodo Online Appendix."),
  r("Table S6", "table", "results/metaboanalyst/proteomics_go/proteomics_go_tableS6.csv",
    paste0(MA, "run-proteomics-go.R"), "results/proteomics_go_uniprot.csv",
    paste0("EXTERNALIZED 2026-08-26: the inline table was removed from supplement_v2.qmd (Editor SS3.14) and replaced with a caption + pointer to results/tables/table_s6_proteomics_go_full.csv (this 182-row scripted output; written by src/pipeline/rebuild-submission-tables.R) for the Zenodo Online Appendix. NOTE the previously-inline qmd table had ~510 rows and predated the 2026-08-12 proteomics regeneration -- externalizing to this canonical artifact also drops that stale row set; the FDR-significant down GO-BP counts driving Fig 6C (27 MISAME-III / 58 Mumta-LW) are unchanged. ",
           "RESOLVED 2026-08-26: manuscript now matches this scripted output (31->27, 65->58 down GO-BP terms), superseding the frozen Manuscript/Table_S6_proteomics_GO.xlsx (Trenton's 2026-08-10 approved run, 31/65) after the underlying proteomics effect estimates were regenerated 2026-08-12 (background grew ~1,196/1,555 -> 1,320/1,681 proteins, matching the individual-protein-analysis counts already cited elsewhere in the manuscript -- see the pre-submission audit for the full diagnosis). The scripted output IS now the source of truth for Table S6 and Fig 6C; re-running this two-hop chain (55-proteomics-go-uniprot.R writes the depends file above, then run-proteomics-go.R writes the Table S6 CSV) should reproduce the manuscript. If it drifts again, it means the upstream proteomics effect estimates changed again -- diagnose there, not here.")),
  r("Table S7", "table", "results/tables/table_s7_primary_secondary_native_units.csv",
    "src/3 visualizations/build_table_s7.R",
    "results/adjusted_combined_arms_intervention_effects_unscaled_results_clean.RDS",
    paste0("Corrected 2026-08-25: previously misattributed to ", A2, "clean_results.R, which is a genuine ",
           "upstream dependency (it writes the depends RDS above) but not the table's actual writer -- ",
           "see Manuscript/CODE_AUDIT_2026-08-25.md #2")),
  r("Table S8", "table", "results/table_s8_cross_compartment.csv",
    paste0(A2, "46-build-table-s8-cross-compartment.R"),
    "results/cross_compartment_proteome_overlap.csv",
    paste("TOP-LEVEL SPRAWL: should move to results/tables/table_s8_cross_compartment.csv.",
          "FIXED 2026-08-26: Panel C (GPX3/SELENOP milk values) had drifted stale -- its input,",
          "results/cross_compartment_proteome_overlap.csv, was dated 2026-08-04, before the",
          "2026-08-12 proteomics-RDS regeneration (same root cause as the Table S6/Mumta-LW fixes).",
          "Panel C read +0.64/+0.63 SD where the manuscript's already-correct printed table (and the",
          "individual-protein Results paragraph) says +0.69/+0.70. Re-ran 13-cross-compartment-",
          "proteome-overlap.R then this script; both now reproduce the printed values exactly. The",
          "printed qmd table itself was never wrong -- only this artifact chain had gone stale.")),
  r("Table S9", "table", "results/blood_chemical_class_enrichment_directional.csv",
    paste0(A2, "52-blood-class-enrichment-direction-sensitivity.R"), "",
    "TOP-LEVEL SPRAWL: should move to results/tables/table_s9_blood_class_enrichment.csv"),
  r("Table S10", "table", "results/tables/table_s10_temporal_persistence.csv",
    paste0(A2, "57-table-s10-temporal-persistence.R"),
    "results/blood_compartment_all_FDRsig_ATE.csv"),
  r("Table S11", "table", "results/tables/table_s11_compartment_pathway.csv",
    paste0(MA, "run-compartment-pathway-trenton.R"),
    "results/compartment_tracking/trenton_pathway_compound_list.csv"),

  # ---- supporting data artifacts the supplement site also serves ----------
  r("Table S11 support", "data", "results/tables/figS10_tracking_counts.csv",
    paste0(FS, "figS10-compartment-tracking.R")),
  r("Compartment tracking", "data", "results/compartment_tracking/trenton_linked_crosscompartment.csv",
    paste0(MA, "run-compartment-pathway-trenton.R"), "results/supplement_status_fdr_features.csv",
    "Trenton's name-based grouping, >=2 compartments (Fig 6D's rule); the earlier union-find version (run-compartment-pathway.R, linked_upregulated_plot_data.csv) is superseded and kept only as figS10-compartment-tracking.R's GROUPING='union-find' sensitivity view, not a manuscript artifact"),
  r("Supplement detection", "data", "results/supplement_status_fdr_features.csv",
    paste0(A2, "54-supplement-detection-status.R"), "",
    "TOP-LEVEL SPRAWL: should move to results/compartment_tracking/")
))

stamp <- function(paths) {
  # a path entry may list several files joined by ";" (e.g. Table S3 = combined +
  # stratified CSV): exists = ALL parts exist; bytes = sum; mtime = oldest part.
  rbindlist(lapply(paths, function(p) {
    parts <- trimws(unlist(strsplit(p, ";", fixed = TRUE)))
    parts <- parts[nzchar(parts)]
    full  <- file.path(root, parts)
    ok    <- file.exists(full)
    if (!length(ok) || !all(ok))
      return(data.table(exists = FALSE, bytes = NA_integer_, mtime = ""))
    info <- file.info(full)
    data.table(exists = TRUE,
               bytes  = as.integer(sum(info$size)),
               mtime  = format(min(info$mtime), "%Y-%m-%dT%H:%M:%S"))
  }))
}

man <- cbind(REGISTRY, stamp(REGISTRY$path))
setcolorder(man, c("exhibit", "kind", "path", "generator", "depends",
                   "exists", "bytes", "mtime", "notes"))
fwrite(man, OUT)

cat("wrote", OUT, "\n")
cat("  artifacts:", nrow(man), "| missing:", sum(!man$exists),
    "| without a generator:", sum(man$generator == ""), "\n")
if (any(!man$exists))
  cat("  MISSING:\n", paste0("   - ", man[exists == FALSE]$path, collapse = "\n"), "\n")
