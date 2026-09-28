# IMiC: intervention effects on human milk composition

Analysis code for:

> Dailey-Chwalibóg, Mertens, et al. **"Nutritional interventions' impacts on human milk:
> three trials in low-resource settings."** *Science* (submitted).

The pipeline estimates the effects of three maternal-nutrition randomized trials on
multi-omics human-milk outcomes. The trials are MISAME-III (Burkina Faso), Mumta-LW
(Pakistan) and ELICIT (Tanzania). Outcomes span macronutrients, micronutrients,
B-vitamins, HMOs, bioactive proteins, targeted and untargeted metabolomics, proteomics,
and the milk microbiome. A cross-compartment extension (MISAME-III) compares the same
signals in maternal blood, milk, and infant blood, by direction and significance.

The paper's interactive online resource, which it cites by section number, is a separate
site: https://amertens.github.io/IMiC-intervention-effects-supplement/ (source:
https://github.com/amertens/IMiC-intervention-effects-supplement).

---

## 1. What is and is not in this repository

The repository holds the analysis code, aggregate results, and reference inputs. It holds
no participant-level data.

**Included**

| Path | Contents |
|---|---|
| `src/`, `functions/`, `figure-scripts/` | The analysis and figure-building code of record |
| `figure-data/SL_vim_plot_data.RDS` | Fig. 1A's classifier summary: cross-validated AUCs by study, visit and arm |
| `metadata/` | Milk/blood component groupings, the component specification, two trial codebooks the covariate cleaning reads, and WHO growth-velocity standards. No participant records. |
| `data/untargeted_annotation/`, `data/milq_age_specific_cutoffs_clean.RDS` | The manual annotation of untargeted milk features, and the MILQ reference cutoffs |
| `results/*_clean.RDS`, `*_clean.csv`, `results/subsetted results/` | Per-outcome effect estimates (point estimate, CI, p, FDR), one row per outcome × study × visit, plus a few aggregate intermediates. No participant records. See `results/README.md`. |
| `results/metaboanalyst/`, `results/compartment_tracking/`, `results/tables/` | The enrichment and pathway results the figures read, and the machine-readable supplementary tables. No participant records. |

**Deliberately excluded**

| Excluded | Why |
|---|---|
| `data/` participant files, `merged_analysis_datasets.RDS` | Individual-level trial data, governed by the IMiC consortium and the three trials' data-access agreements. Not publicly redistributable. |
| Feature-level result files (untargeted metabolomics, blood compartments, the combined-results RDS; 15–250 MB each) and the cross-platform feature-alignment table (`data/additional datasets/IMiC_alignment.csv`) | Size. They carry no participant rows, but they are too large for a code archive. Available on request. |
| `figure-data/`, except Fig. 1A's classifier summary | Saved `ggplot` objects embed their plot data, which can include per-participant rows. |
| Rendered figures, manuscript sources, slide decks | The code regenerates the figures, and the manuscript is published separately. |
| Code that only builds the online resource, exploratory and superseded scripts, development notes and tests | Not needed to reproduce the paper (section 4d). |

**Data availability.** The individual-level datasets underlying these analyses are held
by the IMiC consortium and the MISAME-III, Mumta-LW, and ELICIT trial teams, and are
available under a data-access agreement. See the paper's *Data and materials
availability* statement.

The estimation stage (`src/1 data prep/`, `src/2 analysis/`) needs the restricted data,
so it cannot be re-run from this repository alone. Section 4b lists which figures,
tables, and enrichment steps do rebuild from the shipped files.

---

## 2. Layout

```
src/
  0-config.R                    packages, confounder list (Wvars), outcome groupings,
                                shared aesthetics. Source this first; everything
                                assumes it has run.
  1 data prep/                  per-trial cleaning + merge into the harmonized
                                analytic dataset  [needs restricted data]
  2 analysis/                   TMLE estimation (combined-arms and arm-stratified, plus
                                the unscaled, unadjusted and trajectory variants that
                                clean_results.R assembles), FDR correction, and the
                                MISAME-III blood and cross-compartment analyses
                                [needs restricted data]
  3 visualizations/             Fig. 1A's classifier summary and the Table S1 builder
  metaboanalyst/                scripted MetaboAnalystR enrichment and pathway runs
                                (replaces the manual metaboanalyst.ca workflow)
  run_*.R                       orchestrators for the blood/cross-compartment runs
  README_blood_cross_compartment_pipeline.md   read before touching the blood code

functions/                      analytic helpers: run_bioTMLE() (a drtmle wrapper),
                                extract_bioTMLE_results(), SL learners, plot helpers

figure-scripts/
  0_figure-functions.R          theme_imic(), save_figure_3way(), and label helpers
  manuscript_figures/
    build_all_manuscript_figures.R   entry point: builds every figure
    figN-*.R, figSN-*.R              one generator per figure or panel
    study_colors.R, fig6_layout.R,   shared palette, Fig. 6 geometry, and the renderer
    render_msea_panelB.R             behind Figs. 5B and 6A
    README.md                        figure-by-figure map

metadata/                       component groupings, codebooks, growth standards
data/                           untargeted annotation, MILQ cutoffs
results/                        aggregate effect estimates (see results/README.md)
```

---

## 3. Requirements

- R 4.4.2 (the version used for the submitted analyses).
- Bioconductor + the `tlverse` stack. The full `library()` list is the top of
  `src/0-config.R`; the load-bearing ones are `drtmle`, `SuperLearner`, `sl3`,
  `origami`, `biotmle`, `SummarizedExperiment`, `BiocParallel`, `tidyverse`,
  `data.table`, `here`.
- MetaboAnalystR for the enrichment/pathway stage. The paper used version 4.3.0
  (`src/metaboanalyst/env/sessionInfo.txt`); `src/metaboanalyst/00-setup-metaboanalystr.R`
  installs the current GitHub version. The analyses run locally through the
  package, not the metaboanalyst.ca web tool, but the package downloads its compound
  and pathway libraries from metaboanalyst.ca (section 4b).
- `clusterProfiler` + an org annotation package for the proteome GO analysis.
- Figure scripts: `ggplot2`, `cowplot`, `ggrepel`, `ggforce`, `ggtext`, `patchwork`,
  `magick`, `ragg`, `readxl`, `writexl`, and `RColorBrewer`. Most of them also source
  `src/0-config.R`, so they need its packages too.
- The Mummichog steps shell out to a conda environment named `mummichog`. They call
  `conda` from the `PATH`; point them at another executable without editing code:
  ```bash
  export IMIC_CONDA_CMD=/path/to/conda      # or set IMIC_CONDA_CMD on Windows
  ```

All in-repo paths are resolved with `here::here()` from the repository root. Open
`imic_intervention_effects.Rproj`, or set the working directory to the repo root
before running anything.

**Windows users:** a few generated result filenames are long (up to 108 characters). If
you clone into an already-deep folder you may hit the legacy 260-character path limit
(`Filename too long`). Either clone somewhere shallow, or enable long paths once:

```bash
git config --global core.longpaths true
```

---

## 4. Reproducing the analysis

### 4a. From the restricted data (full pipeline)

```r
source("src/0-config.R")
```
then run, in order: `src/1 data prep/` (1 through 7); the milk analyses in
`src/2 analysis/` (scripts 1 to 10, in numeric order); `src/2 analysis/clean_results.R`,
which scripts 12 onward read; the rest of `src/2 analysis/`; then the figure driver
(section 4b). `src/1 data prep/4-pca-reductions.R` reads the microbiome merge written by
`src/2 analysis/3_adjusted_analysis_microbiome.R`, so run that script first. The blood and
cross-compartment extension has its own orchestrators. Read
`src/README_blood_cross_compartment_pipeline.md` first, then run `src/run_blood_full_adjusted.R`
followed by `src/run_blood_adjusted_downstream.R`.

### 4b. From this archive (no restricted data)

The estimation stage has already been run; its outputs are in `results/`.

```bash
# every figure, in figure-number order (writes to figures/)
Rscript figure-scripts/manuscript_figures/build_all_manuscript_figures.R
Rscript figure-scripts/manuscript_figures/build_all_manuscript_figures.R 2 S3   # a subset

# enrichment steps whose inputs ship here
Rscript src/metaboanalyst/run-untargeted-msea.R            # Table S5 and Fig. 6A's input
Rscript src/metaboanalyst/run-proteomics-go.R              # Table S7
Rscript src/metaboanalyst/build-compartment-query-list.R   # Table S11 query list, then
Rscript src/metaboanalyst/run-compartment-pathway.R        # Table S11 and Fig. 6D's input
```

`build_all_manuscript_figures.R` decides which script builds which figure. It runs each
step in a fresh `Rscript` process, because several generators clear the workspace with
`rm(list=ls())` or change the working directory. Figure key `6` reruns
`run-untargeted-msea.R` before drawing Fig. 6A, so it overwrites the shipped MSEA results
in `results/metaboanalyst/untargeted_msea/` with a rerun against the current
MetaboAnalyst libraries (see below).

Each generator was run on its own in a fresh copy of the release on 2026-09-28; every
figure and table marked as rebuilding came out byte-identical to the original:

| Status | Figures | Tables and other outputs |
|---|---|---|
| Rebuilds as shipped | Figs. 1 and 2; Fig. 3 Panel B; Fig. 6 Panels A, B and C; Figs. S1, S2, S3 and S5 | Tables S7 and S11 |
| Rebuilds, with small differences from library updates (see below) | | Table S5 (about 25 minutes) |
| Needs the restricted participant data | Fig. 4, Fig. S4 | Table S1 |
| Needs an on-request feature-level result file | Fig. 3 (Panel A and the composite), Fig. 5 (Panels A and B and the composite), Fig. 6D (so the Fig. 6 composite builds without Panel D) | Tables S2, S3, S6, S8, S9 and S10 |
| Needs the Biocrates Quant 500 "BioIDs" structure file (a vendor annotation file, not redistributable), as well as the on-request combined-results RDS | Fig. 5 Panel C | Table S4 |

The MetaboAnalystR steps (Tables S5 and S11, and the inputs to Figs. 3B, 5B and 6A)
download MetaboAnalyst's compound and pathway libraries from metaboanalyst.ca, so they
need internet access. MetaboAnalyst updates these libraries without version numbers, so a
rerun uses today's libraries rather than the ones the paper used. On 2026-09-27, a rerun
of the untargeted MSEA reproduced 272 of the 330 rows of `untargeted_msea_combined.csv`
(the full result behind Table S5) exactly. In the other 58 (MISAME-III at 14–21 days and
Mumta-LW at 2 months, both down-regulated), one compound mapped differently, which moved
six pathways across the FDR threshold. The shipped tables are the analyzed results. A
failed download skips the affected cell with a message in the log, so check the row
count against the shipped file.

Figure map, as declared by the driver:

| Fig. | Content |
|---|---|
| 1 | Machine-learning classification of arm from each component class (CV-AUC) (A) + effects on each class's first principal component (B) |
| 2 | Primary-outcome forest |
| 3 | Primary volcano (A) + KEGG pathway impact (B) |
| 4 | HM nutrient distributions vs MILQ reference |
| 5 | Tertiary targeted metabolites (incl. triglyceride panel C) |
| 6 | Exploratory MSEA (A) + Mummichog (B) + proteome GO (C) + cross-compartment comparison (D) |
| S1 to S5 | Growth outcomes, micronutrient-deficiency RR, secondary (HMO/bioactive) forest, triglyceride means, microbiome alpha diversity |

The supplement has Figs. S1 to S5 and Tables S1 to S11; there is no Fig. S6 or beyond.
Generator and output filenames follow the printed figure numbers, so
`figS1-growth-outcomes.R` writes `figures/figureS1_growth_outcomes.png`, which the paper
prints as Fig. S1. `figure-scripts/manuscript_figures/README.md` maps every figure.

### 4c. Tracing a figure or table back to its code

`results/ARTIFACT_MANIFEST.csv` has one row per printed figure, table, and separately
built panel (28 rows). Its columns:

| Column | Meaning |
|---|---|
| `exhibit` | The label the manuscript uses, such as "Fig 6D" or "Table S11" |
| `kind` | `figure` or `table` |
| `path` | Where the artifact is written |
| `generator` | The script that writes it |
| `inputs` | The files that script reads |
| `availability` | `shipped` (rebuilds from this archive), `on request` (needs a feature-level result file), or `restricted` (needs participant data) |

### 4d. What was left out of this release, and why

This repository holds what is needed to reproduce the printed paper: every script on the
path from the raw data to a printed figure or table, or to a number stated in the text,
plus the analyses whose outputs `clean_results.R` assembles into the shared result files.
Everything else stays in the authors' working repository:

- Code that only builds the online resource's pages (subgroup, trajectory, PCA and
  cross-compartment figures, direction-split enrichment panels). The online resource has
  its own repository.
- Exploratory, pilot, and superseded analyses, including earlier versions of the
  cross-compartment matching and the MetaboAnalyst runs that were replaced by the ones
  used in the paper.
- Development material: test suites, validation notes, investigation logs, and
  intermediate result files that no retained script reads.

## 5. Naming conventions

- **`study == "Vital"` is the Mumta-LW trial (Pakistan).** "Vital" is the legacy
  codebase name for the Mumta-LW lactation cohort, and it is the value used throughout
  the analytic datasets, results files, and `studytime` strings; `"Vital-40"` is the
  Mumta-LW 1.5-month visit and `"Vital-56"` the 2-month visit (the numbers are the
  target infant age in days). The manuscript uses Mumta-LW in all public-facing
  labels. The code was not renamed, to avoid breaking downstream RDS keys.
- **`study == "Misame"` is MISAME-III.** Same legacy-short-name convention.
- Analyses come in pooled/combined-arm and arm-stratified variants
  (`*_combined_arms` vs `*_arm_strat`), and in scaled and `*_unscaled` forms. These are
  parallel analyses, so check which one a number came from before comparing across
  files.

---

## 6. Methods and reproducibility notes

- **Estimation.** TMLE through `drtmle::drtmle()` (Benkeser & van der Laan), wrapped as
  `run_bioTMLE()` in `functions/bioTMLE_functions.R`, which is a light adaptation of the
  bioTMLE package interface. `cv_folds = 1` (no cross-fitting) in the main runs.
- **Determinism.** `run_bioTMLE()` sets a fixed seed (`12345`) for each study × visit
  run, and runs execute serially, so estimates reproduce across machines. Figure label
  placement is seeded too; `figure-scripts/manuscript_figures/README.md` notes the one
  caveat (ggrepel's time limit on a much slower machine).
- **Multiple testing.** Benjamini-Hochberg applied within each
  (outcome group × study × visit). A more conservative correction pooled across studies
  and visits within an outcome group is retained as `pval_adj_global` for sensitivity.
- **Confounders.** The adjustment set is `Wvars` in `src/0-config.R`. Gestational age at
  birth, maternal BMI, and mid-upper-arm circumference are deliberately *not* adjusted
  for, since prenatal interventions may affect them.
- **Enrichment.** `src/metaboanalyst/` runs the pathway and enrichment analyses locally
  in MetaboAnalystR in place of the metaboanalyst.ca web tool. For the one cell whose web
  results were saved, the scripted pathway analysis reproduced them exactly and the ORA
  reproduced the set size and hit count (its README has the numbers). Those saved web
  results are not included here.
- **Figure styling.** Every figure exports through `save_figure_3way()` or `ggsave()`
  with the helpers in `figure-scripts/0_figure-functions.R`. Panels have no y-axis
  ticks or horizontal gridlines: most use `theme_imic()`, and the forest plots (Figs. 1,
  2 and S3) use `theme_bw()` with those elements removed.

---

## 7. Citing

Please cite the paper. To cite the code itself, also cite this archive, deposited on
Zenodo at https://doi.org/10.5281/zenodo.22104233.

## 8. License

Released under the MIT License, see [`LICENSE`](LICENSE). You may use, modify, and
redistribute this code, including commercially, provided the copyright notice and
permission notice are retained.

This covers the code and documentation in this repository only. It does not grant any
rights to the underlying IMiC trial data, which is not distributed here and remains
governed by the consortium's data-access agreements (section 1).
