# manuscript_figures/

Build scripts for **every figure of the IMiC *Science* manuscript and its supplement**, consolidated in this one folder. These scripts used to be split across `figure-scripts/`, `figure-scripts/manuscript_figures/`, and `src/3 visualizations/`. They now all live here.

## Filename convention

The published image filename encodes the figure number as the **submitted manuscript** prints it: **Fig N** goes to `figures/figureN.*`, **Fig SN** goes to `figures/figureSN_<slug>.*`, and the generator script name matches (`figN-*.R`, `figSN-*.R`). Artifacts that are *not* numbered exhibits take an `appendix_` prefix (online appendix) or `retired_` instead. `results/ARTIFACT_MANIFEST.csv` is authoritative; regenerate it with `Rscript src/pipeline/build-artifact-manifest.R` after any rename.

Most scripts export through `save_figure_3way()`, which writes PDF, EPS and PNG for the typesetter. All paths resolve with `here::here()` or are relative to the **repo root**, so run every script **from the repo root** (`Rscript figure-scripts/manuscript_figures/<script>.R`).

## Main figures

| Fig | Embedded file (in `manuscript_v4_combined.qmd`) | Script | Output call |
| --- | --- | --- | --- |
| **Fig 1** | `figures/figure1.png` | `fig1-ml-vim-classifier.R` | `save_figure_3way(name="figure1")` |
| **Fig 2** | `figures/figure2.png` | `fig2-primary-forest.R` | `save_figure_3way(name="figure2")` |
| **Fig 3** | `figures/figure3.png` (combined MAIN) plus `figures/figure3_stratified_supplement.png` (arm-stratified SUPPLEMENT) | `fig3-primary-volcano-composite.R` (Panel A volcanoes, emits both figures), `fig3B-pathway.R` (Panel B pathway-impact: combined to `figures/figure3_panelB_msea.png`, stratified to `figures/figure3_panelB_stratified.png`), and `src/metaboanalyst/run-primary-pathway-stratified.R` (stratified Panel B data) | `save_figure_3way(name="figure3"` and `name="figure3_stratified_supplement")` |
| **Fig 4** | `figures/figure4.jpeg` | `fig4-milq-boxplots.R` | `ggsave(here("figures/figure4.jpeg"))` |
| **Fig 5** | `figures/figure5.png` | `fig5-tertiary-composite.R`, plus `fig5C-triglyceride.R` for Panel C | `save_figure_3way(name="figure5")` |
| **Fig 6** | `figures/figure6.png` | `fig6-composite.R` (stitches A+B+C) | panels: `fig6A-untargeted-msea.R` (A), `fig6B-mummichog.R` (B), `fig6C-proteomics.R` (C, with data from `src/2 analysis/55-proteomics-go-uniprot.R`) |

Fig 6's Panel C is the UniProt-native GO over-representation from `src/2 analysis/55-proteomics-go-uniprot.R`, following Trenton's finished method: each protein is counted once, and only down-regulated terms survive FDR.

> **Naming.** As of 2026-09-08 generator *and* output filenames both follow the submitted figure numbers, and the legacy off-by-one names are gone. The renumbering map is below. Note that the old `fig5.R` is a stale `knitr::purl` output with blank Panel B/C placeholders; `fig5-tertiary-composite.R` is the real generator.

## Supplementary figures

The submitted manuscript (aee9284, Version 2) lists **`Figs. S1 to S5` and `Tables S1 to S11`** on the supplement title page. That is the whole numbered supplementary figure set. There is no Fig. S6 or beyond. Since 2026-09-08 the script and output filenames carry those submitted numbers directly.

| Fig | Embedded file | Script | Output call |
| --- | --- | --- | --- |
| **Fig S1** | `figures/figureS1_growth_outcomes.png` | `figS1-growth-outcomes.R` | `save_figure_3way(name="figureS1_growth_outcomes")`; child growth outcomes |
| **Fig S2** | `figures/figureS2_milq_deficiency.png` | `figS2-milq-deficiency.R` | `ggsave(...)`; micronutrient-deficiency RR forest |
| **Fig S3** | `figures/figureS3_hmo_bioactives.png` | `figS3-hmo-bioactives.R` | `save_figure_3way(name="figureS3_hmo_bioactives")`; HMO and bioactive forest |
| **Fig S4** | `figures/figureS4_triglyceride_means.png` | `figS4-triglyceride-means.R` | `save_figure_3way(name="figureS4_triglyceride_means")`; triglyceride treatment means |
| **Fig S5** | `figures/figureS5_microbiome_diversity.png` | `figS5-microbiome-diversity.R` | `save_figure_3way(name="figureS5_microbiome_diversity")`; milk microbiota alpha diversity |

### Renumbering map (2026-09-08)

Every supplementary script and output used to be off by one against the submitted paper. This was a residue of the 2026-08-26 renumbering: the study-design Fig. S1 was cut, which shifted S2 through S6 down to S1 through S5, but only the *labels* moved at the time. The filenames did not. Both now agree.

| Submitted | Old script | Old output | New script | New output |
| --- | --- | --- | --- | --- |
| Fig S1 | `figS2-growth-outcomes.R` | `figureS2_growth_outcomes.*` | `figS1-growth-outcomes.R` | `figureS1_growth_outcomes.*` |
| Fig S2 | `figS3-milq-deficiency.R` | `fig-milq-deficiency-reduction-forest-plot.png` | `figS2-milq-deficiency.R` | `figureS2_milq_deficiency.png` |
| Fig S3 | `supporting/figS4-secondary-forest.R` | `figureS4.*` | `figS3-hmo-bioactives.R` | `figureS3_hmo_bioactives.*` |
| Fig S4 | `figS5-triglyceride-means.R` | `figureS5_triglyceride_means.*` | `figS4-triglyceride-means.R` | `figureS4_triglyceride_means.*` |
| Fig S5 | `figS6-microbiome-diversity.R` | `figureS6_microbiome_diversity.*` | `figS5-microbiome-diversity.R` | `figureS5_microbiome_diversity.*` |
| Fig 3B (panel) | `fig3B-pathway.R` | `figure4_panelB_msea.png` | *(unchanged)* | `figure3_panelB_msea.png` |
| Fig 5C (panel) | `fig5C-triglyceride.R` | `figure4_panelC_tg_composition*.png` | *(unchanged)* | `figure5_panelC_tg_composition*.png` |

Fig. S3's generator also moved up out of `supporting/`, so all five now sit together.

## Not in the submitted supplement

These artifacts are still built and still useful, but they are **not numbered supplementary figures** and must not be cited as such. Each carries an `appendix_` (online appendix) or `retired_` prefix instead of a false S-number.

| Artifact | Script | Status |
| --- | --- | --- |
| `figures/appendix_primary_ora_by_direction.png`, `figures/appendix_tertiary_ora_by_direction.png` | `appendix-ora-by-direction.R` | Online appendix. These are the "up- and down-regulated pathways in the Online Appendix" the Results point to. Was mislabelled **Fig S6a/S6b**. |
| `figures/cross_compartment/appendix_crosscompartment_volcanoes.png` | `appendix-crosscompartment-volcanoes.R` | Online appendix. Was mislabelled **Fig S7**. The printed cross-compartment result is Table S8 (features) and Table S11 (pathways), not a figure. |
| `figures/appendix_blood_class_enrichment.png` | `appendix-blood-class-enrichment.R` | Retired as a figure 2026-08-16 (was Fig S8); the result is printed as **Table S9**. See the writer mismatch below. |
| `figures/cross_compartment/appendix_compartment_tracking.png` | `appendix-compartment-tracking.R` | Retired as a figure 2026-08-16 (was Fig S10), but the **script is still required**. It emits `results/compartment_tracking/trenton_pathway_compound_list.csv` (the Table S11 query) and `results/tables/compartment_tracking_counts.csv` (Results-text counts). |
| `figures/retired_figureS1_study_design.{eps,pdf}` | *(static, no generator)* | Study-design schematic, cut from the manuscript 2026-08-26. The `.png` lives in `archive/figures/`. |
| `figures/retired_figureS10_proteomics_individual.{eps,pdf}` | *(archived)* | Individual-protein plot, cut 2026-08-26. |

> **Blood chemical-class writer mismatch.** `appendix-blood-class-enrichment.R` writes `results/blood_class_enrichment_figure.png`. Note the directory (`results/`, not `figures/`) and the different basename. The published `figures/appendix_blood_class_enrichment.png` is a hand-copied, renamed version of that output, so no script step writes the published filename directly. Its input CSV comes from `src/2 analysis/52-blood-class-enrichment-direction-sensitivity.R`.

## `supporting/`

Auxiliary panel builders and drafts, kept for reference and not wired into the main render:

- `fig1c_subpanel.R`: standalone build of the age-at-collection violin sub-panel.
- `figSX_untargeted_metabolomics_heatmap.R`: supplementary triglyceride heatmap.
- `fig2_andrew_draft.Rmd`: exploratory Fig 1/2 alternative, not used in the submission.

## Shared helpers and inputs

- `../0_figure-functions.R`, that is `figure-scripts/0_figure-functions.R`, holds the `tableau10` palette, `theme_imic()`, `save_figure_3way()`, and the `science_dims` width presets. `src/0-config.R` sources it through `here()`, and most scripts load that; the `fig{1,2,4,5}` scripts source it directly. It stays in `figure-scripts/` **on purpose**, because many scripts reference it by that path.
- Inputs are read from the repo's `results/`, `data/`, `figure-data/`, and `figures/` directories through `here::here()`. The key RDS files are `figure-data/SL_vim_plot_data.RDS` (Panel A of `fig1-ml-vim-classifier.R`, which includes the ELICIT infant-azithromycin **negative control**, built by `src/3 visualizations/5-SL_VIM_plots.R` from `results/SL_individual_lab_vim_res.RDS`), `figures/figure-data/pca_intervention_effects_results.RDS` (Panel B of the same figure), `results/combined_intervention_effects_results_{stratified,combined}_arms.RDS` (the volcanoes in `fig3-primary-volcano-composite.R` and `fig5-tertiary-composite.R`), and `data/merged_analysis_datasets.RDS`.
- `input data/` holds cached RDS copies used by some of the `.Rmd` twins.

## How to render (from the repo root)

The easiest route is the single driver. It rebuilds everything in figure-number order with the full data, panel and composite chain per figure, and it is the source of truth for the mapping between filenames and figure numbers:

```bash
Rscript figure-scripts/manuscript_figures/build_all_manuscript_figures.R          # all figures
Rscript figure-scripts/manuscript_figures/build_all_manuscript_figures.R 1 3 6    # only Figs 1, 3, 6
```

Or run individual scripts:

```bash
Rscript figure-scripts/manuscript_figures/fig1-ml-vim-classifier.R          # Fig 1
Rscript figure-scripts/manuscript_figures/fig2-primary-forest.R             # Fig 2
Rscript figure-scripts/manuscript_figures/fig3B-pathway.R                   # Fig 3 panel B
Rscript figure-scripts/manuscript_figures/fig3-primary-volcano-composite.R  # Fig 3
Rscript figure-scripts/manuscript_figures/fig4-milq-boxplots.R              # Fig 4
Rscript figure-scripts/manuscript_figures/fig5C-triglyceride.R              # Fig 5 panel C
Rscript figure-scripts/manuscript_figures/fig5-tertiary-composite.R         # Fig 5
Rscript figure-scripts/manuscript_figures/fig6A-untargeted-msea.R           # Fig 6 panel A
Rscript figure-scripts/manuscript_figures/fig6-composite.R                  # Fig 6
Rscript figure-scripts/manuscript_figures/figS1-growth-outcomes.R           # Fig S1
Rscript figure-scripts/manuscript_figures/figS2-milq-deficiency.R           # Fig S2
Rscript figure-scripts/manuscript_figures/figS3-hmo-bioactives.R            # Fig S3
Rscript figure-scripts/manuscript_figures/figS4-triglyceride-means.R        # Fig S4
Rscript figure-scripts/manuscript_figures/figS5-microbiome-diversity.R      # Fig S5
```

Online-appendix and retired artifacts, which are not numbered exhibits:

```bash
Rscript figure-scripts/manuscript_figures/appendix-ora-by-direction.R           # ORA by direction
Rscript figure-scripts/manuscript_figures/appendix-crosscompartment-volcanoes.R # cross-compartment volcanoes
Rscript figure-scripts/manuscript_figures/appendix-blood-class-enrichment.R     # then copy to figures/appendix_blood_class_enrichment.png
Rscript figure-scripts/manuscript_figures/appendix-compartment-tracking.R       # still required for the Table S11 query list
```

Label placement is seeded (`seed = 123` on the ggrepel layers), so re-running any of these reproduces its output byte for byte.

## Science-journal compliance

All four `save_figure_3way()` main figures render at the **Science full-page** standard (`science_dims$full_page` = 7.25 in / 184 mm wide, height <= 9.5 in). They meet the Reviewer 2 section 2.6 readability floors (axis-tick >= 7 pt, axis-title >= 8 pt, sub-panel and strip text >= 8 pt, legend >= 7 pt) and the Science >= 5 pt minimum. Verify after any theme change:

```r
source("figure-scripts/0_figure-functions.R")
str(science_dims)              # dimension presets
formals(save_figure_3way)      # defaults: width = 7.25, height = 9.5
```
