# Blood / cross-compartment / mummichog pipeline

This is the maternal/infant **blood** extension of the milk bioTMLE pipeline, built for
the IMiC *Science* revision: it estimates BEP-vs-control intervention effects (ATEs) in
maternal plasma, maternal/infant dried-blood (VAMS) metabolomics, and maternal blood
proteomics, then asks whether the milk effects are mirrored along the chain
**maternal diet → maternal blood → milk → infant blood**.

All scripts assume the working directory is the repo root and that `0-config.R` (the
tlverse stack + `extract_res`, `ci_to_pvalue`, `Wvars`) is sourced first. Rscript lives
at `C:/Program Files/R/R-4.4.2/bin/Rscript.exe` (not on PATH).

## Shared helpers — read this first

**`src/2 analysis/_blood_helpers.R`** holds everything that used to be copy-pasted across
scripts: the matching tolerances, the m/z–RT extractors, the cross-compartment matcher,
the threshold-free concordance methods, and the bioTMLE result tidier. Numbered scripts
`source()` it instead of redefining these. If you change a tolerance or a method, change
it **here**, once.

## Pipeline order

| # | Script | Does | Key output |
|---|--------|------|-----------|
| 7 | `1 data prep/7-blood-compartment-prep.R` | Links lab samples → trial arms (+ milk-ID covariate bridge); builds the 4-level `arm` and period-aware binary `class`; splits postnatal VAMS into maternal/infant | `data/blood/merged_blood_datasets.RDS` |
| 12 | `2 analysis/12-blood-compartment-intervention-effects.R` | bioTMLE ATEs, **arm-stratified** (4 cells vs Control) | `…intervention_effects_results.RDS` |
| 12b | `2 analysis/12-…_combined_arms.R` | bioTMLE ATEs, **combined-arms** (period-aware binary = primary) | `…combined_arms…_results.RDS` |
| — | `2 analysis/clean_blood_results.R` | Tidies raw bioTMLE → long; **FDR (BH) per visit × dataset** | `…_results_clean.RDS` |
| — | `2 analysis/combine_blood_results.R` | Master table + slim FDR-sig CSV | `blood_compartment_all_*` |
| 13 | `2 analysis/13-cross-compartment-proteome-overlap.R` | Milk ↔ maternal-blood proteins by **UniProt**; selenoproteins | `cross_compartment_proteome_*` |
| 14 | `2 analysis/14-cross-compartment-metabolite-mzrt-match.R` | Milk↔blood metabolites by **m/z–RT (25 ppm)**; maternal↔infant by shared ID | `cross_compartment_metab_*` |
| 15 | `2 analysis/15-blood-mummichog-pathway-analysis.R` | Mummichog pathway enrichment on blood, per compartment × ion mode | `blood_mummichog_pathways*` |
| 16 | `2 analysis/16-milk-mummichog-for-comparison.R` | Same for milk (comparison) | milk mummichog tables |
| 17 | `3 visualizations/17-cross-compartment-pathway-figure.R` | Milk vs maternal vs infant pathway heatmap | `figures/cross_compartment/pathway_*` |
| 18 | `2 analysis/18-directional-mummichog.R` | **Directional** mummichog (up/down feature subsets) → pathway directions | `directional_pathways.csv` |
| 19 | `2 analysis/19-cross-compartment-arrow-contrast-summary.R` | Trenton's **arrow × contrast** table: n-sig, pairs ≤25 ppm, % concordant (bulk + FDR-strong). **Demoted to supporting** (ppm-only) | `cross_compartment_arrow_contrast_summary.csv` |
| 19b | `2 analysis/19b-cross-compartment-fdr-first-lists.R` | **PRIMARY (FDR-first reorientation):** small per-arrow × per-contrast lists of features FDR-sig in BOTH, matched (25 ppm / id), **annotated** individually | `cross_compartment_fdr_first_lists.csv` |
| 20 | `2 analysis/20-cross-compartment-threshold-free-panel.R` | **Threshold-free** panel (RRHO / GSEA / weighted-r / anchored) over the same arrows. Supporting (ppm-only) | `cross_compartment_threshold_free_panel.csv` |
| 21 | `2 analysis/21-bep-supplement-to-infant-tracer.R` | **[cross-platform]** supplement (V1) → milk/plasma (V1↔V1 clean) + → infant (V1↔V3 isobaric) directional enrichment of supplement-abundant features among BEP-up features | `bep_supplement_tracer.csv` |
| — | `3 visualizations/blood_compartment_volcano_plots.R` | Per-compartment volcano plots | `figures/blood_volcano/` |

## Orchestrators (in `src/`)

Single-session runners that source `0-config.R` once, set the override flags, and source
the numbered scripts in order. They exist because the untargeted bioTMLE is the slow part
(~hours; 38k-feature VAMS dominate) and `SnowParam` workers must be capped to avoid OOM.

| Runner | What | Notes |
|--------|------|-------|
| `run_blood_full.R` | unadjusted full run (prep → stratified → combined → clean) | `BLOOD_WORKERS` caps PSOCK workers |
| `run_blood_full_adjusted.R` | adjusted full run, writes `adjusted_`-tagged outputs | don't run concurrently with the unadjusted one |
| `run_blood_adjusted_recover.R` | resumes the adjusted run after an OOM at the combined step (2 workers) | |
| `run_ever_prenatal_bep_sensitivity.R` | "ever-prenatal-BEP" sensitivity for maternal compartments | confirmed the combined-arms framing isn't attenuated |
| `archive/src/pilot_threshold_free_concordance.R` | one-off pilot of the 4 threshold-free methods on the identity-anchored pair | **superseded by script 20** (the productionized version); **archived** 2026-06-25 |

## Key conventions

- **`BLOOD_ADJUST`** (default FALSE) switches every script between unadjusted and the
  covariate-adjusted, `adjusted_`-tagged outputs. Adjusted (combined-arms) is the analysis of record.
- **Combined-arms vs stratified.** Combined-arms (period-aware binary) is the powered
  primary; the 4-cell stratified is a sensitivity check (small-n, read qualitatively).
- **Platform / matching.** Milk + prenatal plasma + prenatal VAMS are LC method **V1**;
  postnatal VAMS (maternal + infant) is **V3**. So milk↔blood is cross-platform mass-match
  (isobaric-limited); maternal↔infant VAMS share one V3 catalogue (exact-ID, the
  reviewer-proof comparison). See `Manuscript/cross_compartment_ANALYSIS_PLAN.md`.
- **Measurement constraint (Kim).** Absolute intensities are not comparable across
  datasets — every result is a within-dataset BEP-vs-control effect, compared by *direction*.
