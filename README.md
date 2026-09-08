# IMiC: intervention effects on human milk composition

Analysis code for:

> Dailey-Chwalibóg, Mertens, et al. **"Nutritional interventions' impacts on human milk:
> three trials in low-resource settings."** *Science* (submitted).

The pipeline estimates the effects of three maternal-nutrition randomized trials on
multi-omics human-milk outcomes. The trials are **MISAME-III** (Burkina Faso),
**Mumta-LW** (Pakistan) and **ELICIT** (Tanzania). Outcomes span macronutrients,
micronutrients, B-vitamins, HMOs, bioactive proteins, targeted and untargeted
metabolomics, proteomics, and the milk microbiome. A cross-compartment extension
follows the same signals along maternal diet, maternal blood, milk, and infant blood.

---

## 1. What is and is not in this repository

This is a **code archive**, not a data archive.

**Included**

| | |
|---|---|
| `src/`, `functions/`, `figure-scripts/` | The complete analysis and figure-building code |
| `metadata/` | Data dictionaries, milk/blood component groupings, WHO growth-velocity standards. No participant records. |
| `data/msea/`, `data/untargeted_annotation/` | Aggregate enrichment (ORA/MSEA) exports and metabolite name maps |
| `results/*_clean.RDS`, `*_clean.csv` | Per-outcome **effect estimates** (point estimate, CI, p, FDR) for the named-outcome analyses, one row per outcome x study x visit. No participant records. |
| `results/metaboanalyst/`, `results/compartment_tracking/` | Pathway / ORA / mummichog / GO enrichment output, at pathway and compound level. No participant records. |

**Deliberately excluded**

| | Why |
|---|---|
| `data/` participant files, `merged_analysis_datasets.RDS` | **Individual-level trial data**, governed by the IMiC consortium and the three trials' data-access agreements. Not publicly redistributable. |
| Feature-level untargeted metabolomics and blood-compartment result files (~670 MB) | Size. These carry no participant rows, but they are too large for a code archive. Available on request. |
| `figure-data/*.RDS` | Saved `ggplot` objects, which embed their plot data; three of them carry per-participant rows. Regenerable from the included `src/3 visualizations/` scripts. |
| Rendered figures, tables, manuscript sources, slide decks | Outputs rather than inputs. The code regenerates them, and the manuscript is published separately. |
| Collaborator working trees, exploratory side analyses, superseded scripts | Not part of the analysis of record. |

**Data availability.** The individual-level datasets underlying these analyses are held
by the IMiC consortium and the MISAME-III, Mumta-LW, and ELICIT trial teams, and are
available under a data-access agreement. See the paper's *Data and materials
availability* statement.

**Consequence for reproduction:** the estimation stage (`src/1 data prep/`,
`src/2 analysis/`) requires the restricted data, so it cannot be re-run from this
repository alone. The enrichment and pathway stage does run from what is shipped here,
along with the figure panels that read only its output. Section 4b sets out which
inputs are missing, and how to check any individual figure for yourself.

---

## 2. Layout

```
src/
  0-config.R                    packages, confounder list (Wvars), outcome groupings,
                                shared aesthetics. Source this first; everything
                                assumes it has run.
  1 data prep/                  per-trial cleaning + merge into the harmonized
                                analytic dataset  [needs restricted data]
  2 analysis/                   TMLE estimation (combined-arms and arm-stratified),
                                unscaled + trajectory variants, FDR correction,
                                cross-compartment and mummichog analyses,
                                results cleanup  [needs restricted data]
  3 visualizations/             forest, volcano, and outcome-panel plots
  4 secondary analyses/         secondary-contrast and baseline-VIM analyses
  metaboanalyst/                scripted MetaboAnalystR enrichment/pathway pipeline
                                (replaces the manual metaboanalyst.ca workflow)
  pipeline/                     artifact manifest + staleness checks + table rebuild
  trenton-ports/                local ports of the collaborator analyses that
  trenton-manuscript-analyses/  produce Figs 3B, 6A, 6C and the compartment tracking
  run_*.R                       orchestrators for the blood/cross-compartment runs
  README_blood_cross_compartment_pipeline.md   read before touching the blood code

functions/                      analytic helpers: run_bioTMLE() (a drtmle wrapper),
                                extract_bioTMLE_results(), SL learners, plot helpers

figure-scripts/
  0_figure-functions.R          shared theme_imic() + save_figure_3way() export helper
  manuscript_figures/
    build_all_manuscript_figures.R   >>> single entry point for every figure <<<
    fig1-*.R ... fig6-*.R            one generator per main-text figure/panel
    figS1-*.R ... figS5-*.R          one generator per supplementary figure
    appendix-*.R                     online-appendix panels (not numbered exhibits)

metadata/                       data dictionaries, component specs, growth standards
data/                           aggregate enrichment exports + annotation name maps
results/                        aggregate effect estimates (see results/README.md)
```

---

## 3. Requirements

- **R 4.4.2** (the version used for the submitted analyses).
- Bioconductor + the `tlverse` stack. The full `library()` list is the top of
  `src/0-config.R`; the load-bearing ones are `drtmle`, `SuperLearner`, `sl3`,
  `origami`, `biotmle`, `SummarizedExperiment`, `BiocParallel`, `tidyverse`,
  `data.table`, `here`.
- **MetaboAnalystR 4.3.0** for the enrichment/pathway stage. Install it with
  `src/metaboanalyst/00-setup-metaboanalystr.R`. The pipeline runs fully
  **local and offline**, with no metaboanalyst.ca calls.
- `clusterProfiler` + an org annotation package for the proteome GO analysis.
- The mummichog steps shell out to a **conda** environment named `mummichog`. The conda
  executable defaults to a Windows path; override it without editing code:
  ```bash
  export IMIC_CONDA_CMD=/path/to/conda      # or set IMIC_CONDA_CMD on Windows
  ```

All in-repo paths are resolved with `here::here()` from the **repository root**. Open
`imic_intervention_effects.Rproj`, or set the working directory to the repo root
before running anything.

Every R file in this release parses cleanly, checked with `parse()` across the tree.

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
then run, in order: `src/1 data prep/` (1 through 7), `src/2 analysis/` (numeric order,
with `clean_results.R` last), then `src/3 visualizations/` and the figure scripts.
The blood and cross-compartment extension has its own orchestrators. Read
`src/README_blood_cross_compartment_pipeline.md` first, then use `src/run_blood_full.R`
and `src/run_blood_full_adjusted.R`.

### 4b. From the shipped results (no restricted data)

`clean_results.R` and the enrichment stage have already been run; their outputs are in
`results/`, `results/metaboanalyst/`, and `results/compartment_tracking/`.

```bash
# scripted MetaboAnalystR enrichment / pathway analysis (re-derives results/metaboanalyst/)
Rscript src/metaboanalyst/run-primary.R

# every manuscript figure, in figure-number order
Rscript figure-scripts/manuscript_figures/build_all_manuscript_figures.R

# or a subset
Rscript figure-scripts/manuscript_figures/build_all_manuscript_figures.R 3 6
```

`build_all_manuscript_figures.R` is the **source of truth** for which script builds
which figure. It declares the full data, panel, and composite chain per figure, and it
runs each step in a fresh `Rscript` subprocess, because the generators call
`rm(list=ls())` and so must not share a session.

**Be precise about what rebuilds here.** Not every generator runs from this archive as
shipped. Some stop on an input that is deliberately not included, for one of three
reasons:

1. **They need the restricted participant data.** `fig4-milq-boxplots.R` and
`fig1-ml-vim-classifier.R` read `data/merged_analysis_datasets.RDS`, which cannot be
made public. See the data-availability note in section 1.
2. **They need an excluded large result file**, one of
`combined_intervention_effects_results_*`, `*_untargeted_results_clean_ATE.RDS`, or
`blood_compartment_*_clean.RDS` (41 to 103 MB each). These are aggregate, but too large
to archive. Available on request.
3. **They need a `figure-data/*.RDS` intermediate**, which this archive does not ship.
Those files are saved `ggplot` objects, and a ggplot embeds its plot data; three of them
carry per-participant rows, and `figure1c_subplot.RDS` includes `subjid` for all 1,545
participants. Rather than vet each embedded layer, the whole directory is excluded. They
are regenerable: the `src/3 visualizations/` scripts that write them are included, and
re-running those with the restricted data recreates them.

Broadly, the pathway and enrichment panels read only the aggregate enrichment outputs
that do ship here, while the participant-level and feature-level panels need the
restricted or on-request data. You do not have to take that summary on trust. Each row
of `results/ARTIFACT_MANIFEST.csv` names the artifacts its generator reads in the
`depends` column, so for any given figure you can check directly whether its inputs are
present in this archive.

Figure map, as declared by that driver:

| Fig | Content |
|---|---|
| 1 | ML variable-importance / CV-AUC classifier |
| 2 | Primary-outcome forest |
| 3 | Primary volcano (A) + KEGG pathway-impact (B) |
| 4 | HM nutrient distributions vs MILQ reference |
| 5 | Tertiary targeted metabolites (incl. triglyceride panel C) |
| 6 | Exploratory MSEA (A) + mummichog (B) + proteome GO (C) + cross-compartment pathways (D) |
| S1 to S5 | Growth outcomes, micronutrient-deficiency RR, secondary (HMO/bioactive) forest, triglyceride means, microbiome alpha diversity |

The submitted supplement lists **Figs. S1 to S5 and Tables S1 to S11**. There is no
Fig. S6 or beyond. Artifacts that are not numbered exhibits carry an `appendix_` prefix
(they appear in the online appendix at https://doi.org/10.5281/zenodo.22104233) or a
`retired_` prefix, and their generators are keyed `A-*` or `X-*` in the driver rather
than given a figure number. The direction-split ORA panels and the cross-compartment
volcano panel are `appendix_`.

Generator and output filenames both follow the submitted figure numbers, so
`figS1-growth-outcomes.R` writes `figures/figureS1_growth_outcomes.png`, which the
paper prints as **Fig. S1**. The off-by-one legacy names were removed on 2026-09-08.
The old-to-new table is in `figure-scripts/manuscript_figures/README.md`.

---

## 4c. Tracing a figure or table back to its code

`results/ARTIFACT_MANIFEST.csv` is the authoritative map. It holds one row per
manuscript exhibit (25 figures, 11 tables, 3 data artifacts), with:

| Column | Meaning |
|---|---|
| `exhibit` | The label the manuscript uses, such as "Fig 6D" or "Table S11" |
| `path` | Where the artifact is written |
| `generator` | The **one** script that writes it |
| `depends` | Upstream artifacts it is derived from |

Start there rather than guessing from filenames. Since 2026-09-08 the two agree, and
`fig3B-pathway.R` writes `figures/figure3_panelB_msea.png`. The manifest is still the
authority, because it also records which artifacts are online-appendix or retired
rather than numbered exhibits.

**Two known gaps.** The manifest names two generators that do not exist in the source
repository at all: `figure-scripts/manuscript_figures/figS10-proteomics-individual.R`
and `src/2 analysis/45-figure-blood-transfer-supplement.R`. Both are pre-existing
provenance gaps, not omissions from this release.

## 4d. What was left out of this release, and why

The repository ships the analysis of record. Removed: pilot analyses, figures made for
collaborator memos, data handoff scripts, unfinished drafts, a stale `knitr::purl`
output with blank panels, superseded gene-level proteomics (replaced by the UniProt
version used in the paper), and a separate secondary-contrast ("svn") thread belonging
to a different paper. Each removal was checked against that script's own header.

Scripts on the path from raw data to a manuscript exhibit were all kept, including
cases where a stricter cut might have looked possible. Much of the blood and
cross-compartment pipeline assembles its output filenames at runtime. For example,
`12-blood-compartment-intervention-effects.R` writes
`paste0(here::here(), "/results/blood_compartment_", .tag, "...")`, so no static
analysis can prove those scripts unused. They are retained deliberately.

Scripts that had no header comment carry a generated one listing the files they read
and write, recovered from their syntax tree. Those headers state only what the code
does, never why.

## 5. Naming conventions: read this before interpreting results

- **`study == "Vital"` is the Mumta-LW trial (Pakistan).** "Vital" is the legacy
  codebase name for the Mumta-LW lactation cohort, and it is the value used throughout
  the analytic datasets, results files, and `studytime` strings; `"Vital-40"` is the
  Mumta-LW 1.5-month visit and `"Vital-56"` the 2-month visit. The manuscript uses
  **Mumta-LW** in all public-facing labels. The code was not renamed, to avoid breaking
  downstream RDS keys. Treat `"Vital"` and `"Mumta-LW"` as the same study.
- **`study == "Misame"` is MISAME-III.** Same legacy-short-name convention.
- Analyses come in **pooled/combined-arm** and **arm-stratified** variants
  (`*_combined_arms` vs `*_arm_strat`), and in scaled and `*_unscaled` forms. These are
  deliberate parallel analyses rather than duplicates, so check which frame a number
  came from before comparing across files.

---

## 6. Methods and reproducibility notes

- **Estimation.** TMLE through `drtmle::drtmle()` (Benkeser & van der Laan), wrapped as
  `run_bioTMLE()` in `functions/bioTMLE_functions.R`, which is a light adaptation of the
  bioTMLE package interface. `cv_folds = 1` (no cross-fitting) in the main runs.
- **Determinism.** `run_bioTMLE()` sets a fixed seed (`12345`) for each study x visit
  run, and runs execute serially, so estimates reproduce across machines. Figure label
  placement is seeded as well, so re-running a figure script reproduces its output byte
  for byte.
- **Multiple testing.** Benjamini-Hochberg applied within each
  (outcome group x study x visit). A more conservative correction pooled across studies
  and visits within an outcome group is retained as `pval_adj_global` for sensitivity.
- **Confounders.** The adjustment set is `Wvars` in `src/0-config.R`. Gestational age at
  birth and maternal anthropometry are deliberately *not* adjusted for, since prenatal
  interventions may affect them.
- **Enrichment is web-independent.** `src/metaboanalyst/` reproduces the submitted
  pathway results with a local MetaboAnalystR run. `src/metaboanalyst/FIDELITY-REPORT.md`
  documents the numeric check against the one downloaded ground-truth
  metaboanalyst.ca result, which it matches exactly.
- **Figure styling.** All manuscript panels use `theme_imic()` in
  `figure-scripts/0_figure-functions.R` (Helvetica, minimum font sizes, 7.25 in
  full-page width) and export through `save_figure_3way()`.

---

## 7. Citing

Please cite the paper. This archive is deposited on Zenodo with its own DOI; cite it
alongside the paper when referring to the code specifically.

## 8. License

Released under the **MIT License**, see [`LICENSE`](LICENSE). You may use, modify, and
redistribute this code, including commercially, provided the copyright notice and
permission notice are retained.

This covers the **code and documentation in this repository only**. It does not grant
any rights to the underlying IMiC trial data, which is not distributed here and remains
governed by the consortium's data-access agreements (section 1).
