# `results/`: aggregate effect estimates

These files are the output of the estimation stage (`src/2 analysis/`) and the enrichment
steps (`src/metaboanalyst/`). Each row is one estimated intervention effect, pathway, or
summary statistic. There are no participant records here: no IDs and no per-subject
measurements.

They are shipped so that the steps downstream of estimation (figures, tables, and some
enrichment runs) can be re-run without the restricted individual-level trial data.
`ARTIFACT_MANIFEST.csv` lists every printed figure and table with the script that makes
it, its inputs, and whether it rebuilds from this archive (`shipped`), needs an
on-request result file (`on request`), or needs the participant data (`restricted`).

## Effect-estimate files

| File | Contents | Used by |
|---|---|---|
| `adjusted_combined_arms_intervention_effects_results_clean.RDS` | Covariate-adjusted effects with intervention arms pooled into one contrast per trial (the paper's main analysis). 13,299 rows. | Fig. 4 |
| `adjusted_combined_arms_intervention_effects_results_clean.csv` | CSV copy of the file above, for reading without R. | |
| `adjusted_intervention_effects_results_clean.RDS` | The same models with each randomized arm compared with control. 29,320 rows. | Fig. S4 |
| `adjusted_combined_arms_intervention_effects_unscaled_results_clean.RDS` | Pooled-arm effects on the native measurement scale (not Z-scored). | Table S1 |
| `subsetted results/*.csv` | Pooled-arm estimates split by outcome group: primary macronutrients, micronutrients and B-vitamins; secondary HMOs and bioactive proteins; tertiary targeted metabolites. | Figs. 2 and S3 |

Fig. 5C and Table S4 default to the pooled-arm frame; Table S4 also reports the
arm-stratified MISAME-III results. The proteome GO analysis (Fig. 6C, Table S7) reads
both frames. Check which frame a number came from before comparing across files.

## Schema

The `*_results_clean` files share 24 columns:

| Column | Meaning |
|---|---|
| `study` | `Misame` = MISAME-III, `Vital` = Mumta-LW, `Elicit` = ELICIT. See the naming note in the top-level README. |
| `visit`, `studytime` | Visit label (`1 mo.`, `14-21 days`, ...) and study plus visit (`Elicit (1 mo.)`, `Vital (1.5 mo.)`, ...) |
| `contrast` | Treatment contrast estimated |
| `biomarker`, `label`, `label_f`, `description` | Outcome identity and display labels |
| `outcome_group`, `category`, `category_raw` | Primary / secondary / tertiary grouping and assay category |
| `measure` | Effect measure |
| `est`, `cil`, `ciu` | Point estimate and 95% confidence interval |
| `pval`, `chi_pval` | Wald and chi-square p-values |
| `pval_adj`, `chi_pval_adj` | Benjamini–Hochberg FDR within (outcome group × study × visit) |
| `pval_adj_global` | More conservative BH pooled across studies and visits within an outcome group (sensitivity) |
| `sig`, `sigFDR`, `chi_sig`, `chi_sigFDR` | Significance flags before / after FDR |

## Other files

| File | Used by |
|---|---|
| `growth_intervention_effects_results.RDS` | Fig. S1 |
| `milq_deficiency_reduction_analysis_results.RDS` | Fig. S2 |
| `fat_adjusted_metabolomics_intervention_effects_results.RDS` | Fig. S4 |
| `microbiome_diversity_intervention_effects_results.RDS` | Fig. S5 |
| `pca_intervention_effects_results.RDS` | Fig. 1B |
| `proteomics_go_uniprot.csv` | Fig. 6C and Table S7 |
| `milk_nominal_putative_annotation.csv` | Fig. 6A and Table S5 (putative names for nominally significant untargeted milk features) |
| `supplement_status_fdr_features.csv` | Fig. 6D and Table S11 (FDR-significant cross-compartment features and whether each was detected in the BEP supplement) |
| `table_s8_cross_compartment.csv`, `blood_chemical_class_enrichment_directional.csv` | Tables S8 and S9 |

## Subdirectories

| Directory | Contents |
|---|---|
| `metaboanalyst/` | Enrichment and pathway results behind Figs. 3B, 5B, 5C, 6A and 6B and Tables S2–S7 (see `src/metaboanalyst/README.md`). |
| `compartment_tracking/` | Cross-compartment linkage behind Fig. 6D and Table S11: the linked up-regulated features, their pathway labels, and the compound list sent to the pathway analysis. |
| `tables/` | Tables S1, S10 and S11, and the full-table downloads the supplement names. `table_s3_tertiary_msea_full.csv`, `table_s6_mummichog_full.csv` and `table_s7_proteomics_go_full.csv` are formatted copies of the runner outputs in `metaboanalyst/` (Table S3: same values as `tertiary_msea/tertiary_msea_dual.csv`; Table S6: the pathways in `mummichog/milk_mummichog_pathways.csv` with P < 0.05; Table S7: identical to `proteomics_go/proteomics_go_pathways.csv`). |

## Not included

Feature-level results (untargeted metabolomics, and the blood-compartment and
cross-compartment analyses) are the same kind of aggregate estimate, but they are too
large for a code archive (15–250 MB each). The steps that read them (Panel A of Figs. 3
and 5, Fig. 6D, and Tables S2, S3, S6, S8, S9 and S10) will not rebuild from this
repository alone. The files are available from the authors on request.
