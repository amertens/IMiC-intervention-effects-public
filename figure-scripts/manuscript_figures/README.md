# manuscript_figures/

Build scripts for **every figure of the IMiC *Science* manuscript and its supplement**, consolidated in this one folder. Historically these scripts were split across `figure-scripts/`, `figure-scripts/manuscript_figures/`, and `src/3 visualizations/`; they now all live here.

## Filename convention

The published image filename encodes the figure number: **Fig N -> `figures/figureN.*`**, **Fig SN -> `figures/figureSN.*`** (a few supplement panels keep descriptive names). The generator script name does **not** always match the figure number (historical), so always go by the `save_figure_3way(name=...)` / `ggsave(filename=...)` target, documented below.

Most scripts export via `save_figure_3way()` (PDF + EPS + PNG for the typesetter). All paths are resolved with `here::here()` or are relative to the **repo root**, so run every script **from the repo root** (`Rscript figure-scripts/manuscript_figures/<script>.R`).

## Main figures

| Fig | Embedded file (in `manuscript_v2_combined.qmd`) | Script | Output call |
| --- | --- | --- | --- |
| **Fig 1** | `figures/figure1.png` | `fig1-ml-vim-classifier.R` | `save_figure_3way(name="figure1")` |
| **Fig 2** | `figures/figure2.png` | `fig2-primary-forest.R` | `save_figure_3way(name="figure2")` |
| **Fig 3** | `figures/figure3.png` (combined MAIN) + `figures/figure3_stratified_supplement.png` (arm-stratified SUPPLEMENT) | `fig3-primary-volcano-composite.R` (Panel A volcanoes — emits both figures) + `fig3B-pathway.R` (Panel B pathway-impact: combined -> `figures/figure4_panelB_msea.png`, stratified -> `figures/figure3_panelB_stratified.png`) + `src/metaboanalyst/run-primary-pathway-stratified.R` (stratified Panel B data) | `save_figure_3way(name="figure3"` and `name="figure3_stratified_supplement")` |
| **Fig 4** | `figures/figure4.jpeg` | `fig4-milq-boxplots.R` | `ggsave(here("figures/figure4.jpeg"))` |
| **Fig 5** | `figures/figure5.png` | `fig5-tertiary-composite.R` (+ `fig5C-triglyceride.R` for Panel C) | `save_figure_3way(name="figure5")` |
| **Fig 6** | `figures/figure6.png` | `fig6-composite.R` (stitches A+B+C) | panels: `fig6A-untargeted-msea.R` (A), `fig6B-mummichog.R` (B), `fig6C-proteomics.R` (C — data from `src/2 analysis/55-proteomics-go-uniprot.R`) |

> **Naming caveat.** Generator filenames follow the manuscript figure numbers, but several *output* filenames keep legacy numbers -- e.g. `fig3B-pathway.R` writes `figures/figure4_panelB_msea.png`. `results/ARTIFACT_MANIFEST.csv` is authoritative for which script writes which exhibit. The old `fig5.R` was a stale `knitr::purl` output with blank Panel B/C placeholders and has been removed; `fig5-tertiary-composite.R` is the real generator.

> **Fig 6 status — BUILT.** `fig6-composite.R` stitches Panel A (untargeted MSEA), Panel B (Mummichog), and Panel C (milk proteome GO) into `figures/figure6.png`. Panel C is the UniProt-native GO over-representation from `src/2 analysis/55-proteomics-go-uniprot.R` (Trenton's finished method: each protein counted once; only down-regulated terms survive FDR). No placeholder remains.

## Supplementary figures

| Fig | Embedded file | Script | Notes |
| --- | --- | --- | --- |
| **Fig S1** | `figures/figureS1.png` | **STATIC — the exact submitted image** (do not rebuild) | `fig1.R` rebuilds a re-rendered version to `figures/figureS1_rebuilt_fig1R.png`, kept for reference but NOT embedded. Per author request, S1 is the verbatim submitted study-design figure (only the number changed, Fig 1 → Fig S1). A higher-res original PNG can replace it if located. |
| **Fig S2** | `extracted/media_supplement/media/image1.png` | *(none — pre-rendered static image embedded in the qmd)* | child-growth outcomes |
| **Fig S3** | `figures/fig-milq-deficiency-reduction-forest-plot.png` | `figS3-milq-deficiency.R` | micronutrient-deficiency RR forest |
| **Fig S4** | `figures/figureS4.png` | `supporting/figS4-secondary-forest.R` | `save_figure_3way(name="figureS4")`; HMO + bioactive forest |
| **Fig S5** | `extracted/media_supplement/media/image5.png` | *(none — pre-rendered static image)* | triglyceride treatment means |
| **Fig S6** | `extracted/media_supplement/media/image6.png` | *(none — pre-rendered static image)* | milk microbiota alpha diversity |
| **Fig S7** | `figures/cross_compartment/fig_blood_transfer_supplement.png` | `src/2 analysis/45-figure-blood-transfer-supplement.R` | **not in this folder** — it is a numbered step of the blood cross-compartment analysis pipeline; left in place to preserve pipeline order |
| **Fig S8** | `figures/figureS8_blood_class_enrichment.png` | `figS8-blood-class-enrichment.R` | see mismatch note below |
| **Fig S9** | `figures/figureS_primary_ora_by_direction.png` + `figures/figureS_tertiary_ora_by_direction.png` | `figS9-ora-by-direction.R` | writes both panels (primary + tertiary ORA by direction) |

> **Fig S8 filename mismatch.** `figS8-blood-class-enrichment.R` writes `results/blood_class_enrichment_figure.png` (note: `results/`, not `figures/`, and a different basename). The published `figures/figureS8_blood_class_enrichment.png` is a hand-copied/renamed version of that output — there is no script step that writes the published filename directly. Its input CSV comes from `src/2 analysis/52-blood-class-enrichment-direction-sensitivity.R`.

## `supporting/`

Auxiliary panel builders and drafts kept for reference (not all wired into the main render):

- `figS4-secondary-forest.R` — **active**: builds Fig S4 (`figureS4.png`).
- `fig1c_subpanel.R` — standalone build of the Fig S1 Panel C age-at-collection violin (also inlined by `fig1.R`).
- `figSX_untargeted_metabolomics_heatmap.R` — supplementary triglyceride heatmap.
- `fig2_andrew_draft.Rmd` — exploratory Fig 1/2 alternative (not used in submission).

## Shared helpers and inputs

- `../0_figure-functions.R` (i.e. `figure-scripts/0_figure-functions.R`) — `tableau10` palette, `theme_imic()`, `save_figure_3way()`, and the `science_dims` width presets. Sourced by `src/0-config.R` (via `here()`), which most scripts load; the `fig{1,2,4,5}` scripts source it directly. **Left in `figure-scripts/` on purpose** — many scripts reference it by that path.
- Inputs are read from the repo `results/`, `data/`, `figure-data/`, and `figures/` directories via `here::here()`. Key RDS: `figure-data/SL_vim_plot_data.RDS` (fig1-ml-vim-classifier.R Panel A — includes the ELICIT infant-azithromycin **negative control**; built by `src/3 visualizations/5-SL_VIM_plots.R` from `results/SL_individual_lab_vim_res.RDS`), `figures/figure-data/pca_intervention_effects_results.RDS` (fig1-ml-vim-classifier.R Panel B), `results/combined_intervention_effects_results_{stratified,combined}_arms.RDS` (fig3-primary-volcano-composite.R / fig5-tertiary-composite.R volcanoes), `data/merged_analysis_datasets.RDS`.
- `input data/` — cached RDS copies used by some `.Rmd` twins.

## How to render (from the repo root)

**Easiest: one driver rebuilds everything in figure-number order, with the full
data -> panel -> composite chain per figure (this driver is the source of truth
for the filename<->figure-number mapping):**

```bash
Rscript figure-scripts/manuscript_figures/build_all_manuscript_figures.R          # all figures
Rscript figure-scripts/manuscript_figures/build_all_manuscript_figures.R 1 3 6    # only Figs 1, 3, 6
```

Or run individual scripts (note the filename mismatch — see the tables above):

```bash
Rscript figure-scripts/manuscript_figures/fig1-ml-vim-classifier.R                                   # Fig 1
Rscript figure-scripts/manuscript_figures/fig2-primary-forest.R                   # Fig 2
Rscript figure-scripts/manuscript_figures/fig3-primary-volcano-composite.R                                   # Fig 3
Rscript figure-scripts/manuscript_figures/fig4-milq-boxplots.R                    # Fig 4
Rscript figure-scripts/manuscript_figures/fig5-tertiary-composite.R                # Fig 5
Rscript figure-scripts/manuscript_figures/fig6A-untargeted-msea.R               # Fig 6 panel A
Rscript figure-scripts/manuscript_figures/fig5C-triglyceride.R
Rscript figure-scripts/manuscript_figures/fig1.R                                   # Fig S1
Rscript figure-scripts/manuscript_figures/figS3-milq-deficiency.R        # Fig S3
Rscript figure-scripts/manuscript_figures/supporting/figS4-secondary-forest.R  # Fig S4
Rscript figure-scripts/manuscript_figures/figS8-blood-class-enrichment.R          # Fig S8 (then copy to figures/figureS8_blood_class_enrichment.png)
Rscript figure-scripts/manuscript_figures/figS9-ora-by-direction.R     # Fig S9
```

Fig S7 is produced by the blood pipeline: `Rscript "src/2 analysis/45-figure-blood-transfer-supplement.R"`.

## Science-journal compliance

All four `save_figure_3way()` main figures render at the **Science full-page** standard (`science_dims$full_page` = 7.25 in / 184 mm wide, height <= 9.5 in) and meet the Reviewer 2 §2.6 readability floors (axis-tick >= 7 pt, axis-title >= 8 pt, sub-panel / strip text >= 8 pt, legend >= 7 pt) and the Science >= 5 pt minimum. Verify after theme changes:

```r
source("figure-scripts/0_figure-functions.R")
str(science_dims)              # dimension presets
formals(save_figure_3way)      # defaults: width = 7.25, height = 9.5
```
