# Blood and cross-compartment pipeline

In a MISAME-III subset with paired samples, these scripts estimate BEP-versus-control
average treatment effects (bioTMLE, as for milk) in maternal venous plasma, maternal
prenatal and postnatal dried-blood microsamples (VAMS), infant postnatal VAMS, and
maternal blood proteomics. They then compare the effects with milk by direction and
significance only, because intensities are not comparable across platforms. The
covariate-adjusted, combined-arms (pooled postnatal-BEP vs control) results are the
analysis of record; the arm-stratified results feed the threshold-free check (script 20).

Run everything from the repo root. `src/0-config.R` supplies the tlverse stack and
`extract_res`, `ci_to_pvalue` and `Wvars`. Shared tolerances, m/z-RT catalogues, the
accurate-mass matcher (±25 ppm within ionization mode) and the bioTMLE result tidier
are in `src/2 analysis/_blood_helpers.R`; change a tolerance there, once.

## Run order

1. `Rscript src/run_blood_full_adjusted.R`: prep (`1 data prep/7-blood-compartment-prep.R`),
   both script-12 variants, `clean_blood_results.R`, with `BLOOD_ADJUST = TRUE`.
2. `Rscript src/run_blood_adjusted_downstream.R`: every downstream script below, in
   dependency order.

The bioTMLE runs take hours (the 38k-feature postnatal VAMS dominate); `BLOOD_WORKERS`
caps the parallel workers to avoid running out of memory. `N_FEATURES_SUBSET` (default
`NULL`, all features) can be set before sourcing a script-12 file for a quick test run.

## Scripts (all in `src/2 analysis/`) and what they support

| Script | Produces (in `results/`) | Supports |
|---|---|---|
| `12-blood-compartment-intervention-effects.R`, `12-..._combined_arms.R` | `blood_compartment_[adjusted_][combined_arms_]intervention_effects_results.RDS` | Blood ATEs (Methods) |
| `clean_blood_results.R` | `..._results_clean.RDS` (BH FDR per dataset x visit) | Blood ATEs (Methods) |
| `combine_blood_results.R` | `blood_compartment_all_intervention_effects_results_clean.RDS`, `blood_compartment_all_FDRsig_ATE.csv` | Input to 57 |
| `57-table-s10-temporal-persistence.R` | `tables/table_s10_temporal_persistence.csv` | Table S10 |
| `13-cross-compartment-proteome-overlap.R` | `cross_compartment_proteome_overlap_adjusted.csv` | Input to 46 |
| `46-build-table-s8-cross-compartment.R` | `table_s8_cross_compartment.csv`, `table_s8_fragment.md` | Table S8 (incl. SELENOP and GPX3 raised in milk and blood) |
| `52-blood-class-enrichment-direction-sensitivity.R` | `blood_chemical_class_enrichment_directional.csv`, `blood_plasma_class_match_sensitivity.csv` | Table S9 |
| `15-blood-mummichog-pathway-analysis.R` | `mummichog_output_adjusted/`, `blood_mummichog_pathways_adjusted.csv` | Input to 23 and 32 |
| `32-signed-pathway-direction.R` | `signed_pathway_direction.csv` | Binomial sign test (Methods); fatty-acid synthesis up in maternal blood (Discussion) |
| `16-milk-mummichog-for-comparison.R` | `mummichog_output/` (milk runs) | Input to 33 |
| `33-signed-pathway-direction-milk.R` | `signed_pathway_direction_milk.csv` | Sign test for milk; input to 41 |
| `41-fat-synthesis-timepoint-table.R` | `fat_synthesis_timepoint_table.csv` | De novo fatty-acid synthesis down in milk at each visit (Discussion) |
| `18-directional-mummichog.R` | `mummichog_output_directional/`, `directional_pathways.csv` | Input to 47 |
| `47-annotate-milk-features.R` | `milk_nominal_putative_annotation.csv`, `milk_fdr_sig_putative_annotation.csv` | Table S5 and Fig. 6A (via `src/metaboanalyst/run-untargeted-msea.R`); input to 54 |
| `23-annotate-fdr-features.R` | `fdr_sig_putative_annotation.csv` | Input to 54 |
| `50-blood-mummichog-input.R` | `blood_mummichog_input.{csv,xlsx}`, `mummichog_input_blood/` | Input to 54 |
| `54-supplement-detection-status.R` | `supplement_status_fdr_features.csv` | "16 of the 20" linked features mass-matched to the BEP supplement; Fig. 6D and Table S11 (via `src/metaboanalyst/run-compartment-pathway.R`) |
| `20-cross-compartment-threshold-free-panel.R` | `cross_compartment_threshold_free_panel.csv` | Rank-based, threshold-free agreement check (Methods) |
| `29-figure-overall-shift-agreement.R` | console table; `figures/cross_compartment/fig_overall_shift_agreement.png` | "16 of 18 maternally significant features agreed in direction in the infant" |
| `55-proteomics-go-uniprot.R` | `proteomics_go_uniprot.csv` | Fig. 6C and Table S7 (via `src/metaboanalyst/run-proteomics-go.R`); milk proteome, independent of the blood results |

## Inputs

All of these scripts need data that are not in this archive: the participant-level
blood datasets (`data/blood/`, built by the prep script) [restricted data], the
feature-level milk and blood result files, and the feature catalogues in
`data/additional datasets/` (`ProcessedDataMISAME3_plasma.csv`,
`ProcessedDataMISAME3_VAMS.csv`, `metabolite_description_vam_with_global_id.csv`,
`IMiC_alignment.csv` and the per-study compound-ID tables), which are not shipped.
`metadata/blood_component.Rdata` and
`data/untargeted_annotation/untargeted_to_annotate_manually_annotated.xlsx` are included.

## Mummichog

Scripts 15, 16 and 18 call the `mummichog` Python package through
`conda run -n mummichog`. Create a conda environment named `mummichog` with the package
installed, and set the environment variable `IMIC_CONDA_CMD` to the conda executable if
`conda` is not on PATH (the default is `conda`). Settings match the milk analysis:
10 ppm, `human_mfn` network, significance cutoff p < 0.05, one run per compartment x
visit x ionization mode.

## Conventions

- `BLOOD_ADJUST` (default `FALSE`) switches scripts 12, `clean_blood_results.R`, 13 and 15
  between unadjusted outputs and `adjusted_`-tagged adjusted outputs; the other
  downstream scripts read the adjusted results directly.
- Milk, maternal plasma and prenatal VAMS were run on LC method V1; postnatal VAMS
  (maternal and infant) on V3. Milk-blood metabolite matches are therefore by accurate
  mass, and maternal-infant VAMS matches by shared feature id. Proteins are matched by
  UniProt accession.
- Metabolite identities, supplement matches and pathway assignments are accurate-mass
  based and provisional.
