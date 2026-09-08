# `results/`: aggregate effect estimates

These files are the output of `src/2 analysis/clean_results.R`. Each row is **one
estimated intervention effect** for one outcome, in one study, at one visit. There are
**no participant records here**, no IDs, no per-subject measurements.

They are shipped so that everything downstream of estimation (FDR correction,
enrichment/pathway analysis, figures, tables) can be re-run without the restricted
individual-level trial data.

## Files

| File | Contents |
|---|---|
| `adjusted_intervention_effects_results_clean.RDS` | **Main analysis.** Covariate-adjusted, arm-stratified effects for the named outcomes. |
| `adjusted_combined_arms_intervention_effects_results_clean.RDS` | Same, with treatment arms pooled (combined-arm frame). |
| `adjusted_combined_arms_intervention_effects_results_clean.csv` | CSV of the above, for inspection without R. |
| `unadjusted_intervention_effects_results_clean.RDS` | Unadjusted comparison. |
| `adjusted_intervention_effects_unscaled_results_clean.RDS`<br>`adjusted_combined_arms_intervention_effects_unscaled_results_clean.RDS` | Effects on the native measurement scale rather than standardized. |
| `adjusted_intervention_effects_traj_results_clean.RDS`<br>`adjusted_combined_arms_intervention_effects_traj_results_clean.RDS` | Trajectory (change-over-visits) analyses. |
| `adjusted_combined_arms_intervention_effects_proteomics_results_clean.RDS` | Proteomics, protein-level. |
| `adjusted_combined_arms_and_visits_intervention_effects_results_clean.RDS` | Arms *and* visits pooled. |
| `CHILD_nutrient_medians.csv`, `CHILD_nutrient_z_cutoffs.csv` | CHILD-cohort reference medians and z-score cutoffs used by the MILQ comparison figure. |

**Pick the right frame.** The pooled (`combined_arms`) and arm-stratified files answer
different questions and are both used in the paper, the proteomics results and the
Fig. 5C triglyceride analysis come from the arm-stratified frame, most other panels from
the pooled frame. Check which file a number came from before comparing across outputs.

## Schema

Shared by the `*_results_clean` files (24 columns in the main file, 29,320 rows):

| Column | Meaning |
|---|---|
| `study` | `Misame` = MISAME-III, `Vital` = **Mumta-LW**, `Elicit` = ELICIT. See the naming note in the top-level README. |
| `visit`, `studytime` | Visit label, e.g. `Elicit (1 mo.)`, `Misame (14-21 days)`, `Vital (1.5 mo.)` |
| `contrast` | Treatment contrast estimated |
| `biomarker`, `label`, `label_f`, `description` | Outcome identity and display labels |
| `outcome_group`, `category`, `category_raw` | Primary / secondary / tertiary grouping and assay category |
| `measure` | Effect measure |
| `est`, `cil`, `ciu` | Point estimate and 95% confidence interval |
| `pval`, `chi_pval` | Wald and chi-square p-values |
| `pval_adj`, `chi_pval_adj` | Benjamini–Hochberg FDR **within** (outcome group × study × visit) |
| `pval_adj_global` | More conservative BH pooled across studies and visits within an outcome group (sensitivity) |
| `sig`, `sigFDR`, `chi_sig`, `chi_sigFDR` | Significance flags before / after FDR |

## Subdirectories

| Directory | Contents |
|---|---|
| `metaboanalyst/` | Output of the scripted MetaboAnalystR pipeline (`src/metaboanalyst/`): ORA and pathway-impact cells per outcome group and direction, the mummichog run for Table S5, the proteome GO run for Table S6, and the assembled supplementary tables. Pathway- and compound-level. |
| `compartment_tracking/` | Cross-compartment linkage tables behind Fig. 6D and Fig. S10, matched features, linked up-regulated sets, and the pathway/compound lists. |

Together with the `*_clean` files above, these are what let the pathway and enrichment
figures rebuild without the restricted data. Also present at the top level: a handful of
aggregate intermediates the figure scripts read directly, 
`pathway_replication_matrix.csv`, `proteomics_go_uniprot.csv`,
`blood_chemical_class_enrichment_directional.csv`,
`microbiome_diversity_intervention_effects_results.RDS`,
`growth_intervention_effects_results.RDS`,
`pca_intervention_effects_results.RDS`,
`milq_deficiency_reduction_analysis_results.RDS`,
`fat_adjusted_metabolomics_intervention_effects_results.RDS`, and
`sapient_annotation_handoff.csv`.

## Not included

Feature-level results, untargeted metabolomics and the blood-compartment /
cross-compartment analyses, are the same kind of aggregate estimate but run to roughly
**670 MB** across eight files (1–6.8 million rows each), which is too large for a code
archive. Figures that depend on them (parts of Figs 3, 6, and S6–S7) will not rebuild
from this repository alone. They are available from the authors on request.
