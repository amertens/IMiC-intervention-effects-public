# manuscript_figures/

Build scripts for the printed figures of the IMiC *Science* manuscript (Figs 1–6) and its
supplement (Figs S1–S5). `build_all_manuscript_figures.R` runs them in the right order.

Script and output names carry the printed figure number: **Fig N** is built by `figN-*.R`
and written to `figures/figureN.*`; **Fig SN** is built by `figSN-*.R` and written to
`figures/figureSN_<slug>.*`. Most scripts export through `save_figure_3way()`, which
writes PDF, EPS and PNG. Run every script **from the repo root**
(`Rscript figure-scripts/manuscript_figures/<script>.R`). Each script's header lists its
inputs and outputs and says when it needs restricted or on-request data.

## Exhibits

"On request" = a feature-level result file that is not shipped because of its size;
"restricted" = participant-level data.

| Exhibit | Script(s), in run order | Inputs | Output |
| --- | --- | --- | --- |
| **Fig 1** | `src/3 visualizations/5-SL_VIM_plots.R`, `fig1-ml-vim-classifier.R` | `figure-data/SL_vim_plot_data.RDS` (shipped; rebuilt by `5-SL_VIM_plots.R` from restricted data), `results/pca_intervention_effects_results.RDS` | `figures/figure1.{png,pdf,eps}` |
| **Fig 2** | `fig2-primary-forest.R` | `results/subsetted results/primary_{macro,micro,bvit}.csv` | `figures/figure2.{png,pdf,eps}` |
| **Fig 3** | `src/metaboanalyst/run-primary-pathway-local.R`, `fig3B-pathway.R` (Panel B), `fig3-primary-volcano-composite.R` (Panel A + composite) | `results/metaboanalyst/primary_pathway_local/primary_pathway_all_cells.csv`; combined-results RDS (on request) | `figures/figure3_panelB_msea.png`, `figures/figure3.{png,pdf,eps}` |
| **Fig 4** | `fig4-milq-boxplots.R` | `data/merged_analysis_datasets.RDS` (restricted), `data/milq_age_specific_cutoffs_clean.RDS`, `results/adjusted_combined_arms_intervention_effects_results_clean.RDS` | `figures/figure4.jpeg` |
| **Fig 5** | `src/metaboanalyst/run-tertiary-msea-dual.R`, `fig5C-triglyceride.R` (Panel C), `fig5-tertiary-composite.R` (Panels A, B + composite) | `results/metaboanalyst/tertiary_msea/tertiary_msea_dual.csv`; combined-results RDS (on request); Biocrates Quant 500 "BioIDs" structure file for Panel C (vendor file, not redistributable) | `figures/figure5_panelB_msea.png`, `figures/figure5_panelC_tg_composition{,_nolegend}.png`, `figures/figure5.{png,pdf,eps}` |
| **Fig 6** | `fig6A-untargeted-msea.R`, `fig6B-mummichog.R`, `fig6C-proteomics.R`, `fig6D-crosscompartment.R`, then `fig6-composite.R` (upstream steps: see the driver) | `results/metaboanalyst/untargeted_msea/untargeted_msea_combined.csv`, `results/metaboanalyst/mummichog/milk_mummichog_pathways.csv`, `results/proteomics_go_uniprot.csv`, `results/compartment_tracking/{linked_crosscompartment,metabolite_pathways}.csv`; milk untargeted ATE RDS for Panel D (on request) | `figures/figure6_panel{A,B,C,D}_*.png`, `figures/figure6_legend_ABC.png`, `figures/figure6.{png,pdf}` |
| **Fig S1** | `figS1-growth-outcomes.R` | `results/growth_intervention_effects_results.RDS` | `figures/figureS1_growth_outcomes.{png,pdf,eps}` |
| **Fig S2** | `figS2-milq-deficiency.R` | `results/milq_deficiency_reduction_analysis_results.RDS` | `figures/figureS2_milq_deficiency.png` |
| **Fig S3** | `figS3-hmo-bioactives.R` | `results/subsetted results/secondary_{hmo,bioactives}.csv` | `figures/figureS3_hmo_bioactives.{png,pdf,eps}` |
| **Fig S4** | `figS4-triglyceride-means.R` | `results/adjusted_intervention_effects_results_clean.RDS`, `results/fat_adjusted_metabolomics_intervention_effects_results.RDS`, `data/merged_analysis_datasets.RDS` (restricted) | `figures/figureS4_triglyceride_means.{png,pdf,eps}` |
| **Fig S5** | `figS5-microbiome-diversity.R` | `results/microbiome_diversity_intervention_effects_results.RDS` | `figures/figureS5_microbiome_diversity.{png,pdf,eps}` |

`fig5C-triglyceride.R` also writes the Table S4 CSVs
(`results/metaboanalyst/triglyceride_fa/triglyceride_fa_composition_{combined,stratified}.csv`).
Fig 6C uses the UniProt-level GO over-representation from
`src/2 analysis/55-proteomics-go-uniprot.R`: each protein is counted once, and only
down-regulated terms reach FDR significance.

## Shared helpers

| File | Holds |
| --- | --- |
| `../0_figure-functions.R` | `tableau10`, `canonical_label()`, `plotmath_label()` (B-vitamin subscripts in PDF/EPS), `imic_logp_title` (the shared "–Log₁₀(*P*-value)" axis title), `theme_imic()`, `science_dims`, `save_figure_3way()`, `get_bh_cutoff()`. |
| `study_colors.R` | Study colours (Okabe-Ito: ELICIT blue, MISAME-III orange, Mumta-LW reddish purple), a symbol per study, and the category palette of the Fig 3A / 5A volcanoes. |
| `render_msea_panelB.R` | The renderer behind Figs 5B and 6A. |
| `fig6_layout.R` | Fig 6 panel sizes and print-scale factors. `fig6-composite.R` prints the actual scales on every run and warns if they drift from these constants. |

## How to render (from the repo root)

```bash
Rscript figure-scripts/manuscript_figures/build_all_manuscript_figures.R            # all figures
Rscript figure-scripts/manuscript_figures/build_all_manuscript_figures.R 2 S1 S3    # only these keys
```

The driver's keys are `1`–`6` and `S1`–`S5`. It runs each step in a fresh `Rscript`
process, creates `figures/` and `figure-data/` if they are missing, and prints a summary of
failed steps. To run one script by hand, run its upstream steps first (see the table).

## Reproducibility of label placement

ggrepel places labels at draw time, so every repel layer sets its own seed (`seed = 123`
in Figs 3A, 3B, 5A, 5B, 5C and 6A; `seed = 1` in Figs 6B and 6C). Re-running a script on
the same machine therefore reproduces its output byte for byte. Figs 3B, 5A, 5B and 6A
also set `max.time = Inf`, so placement stops on the iteration count alone and is the same
on any machine. Figs 3A, 5C, 6B and 6C keep ggrepel's default 0.5-second time limit, so a
much slower computer can place a crowded label slightly differently.

## Sizes and theme

- **Theme.** Panels have no y-axis ticks and no horizontal gridlines. Most use
  `theme_imic()`; the forest plots (Figs 1, 2 and S3) use `theme_bw()` with white facet
  strips and those elements removed.
- **Width.** Most figures are drawn 7.25 in wide (Science's two-column width). Fig 4
  (10 in) and Figs S2 (8 in), S4 and S5 (9 in) are drawn wider and are scaled down to fit.
- **Height.** Figs 1, 2, 3 and 5 are 10.6–10.7 in tall, more than the 9.5 in page, so they
  print at about 89% of their drawn size. Fig 6 is assembled from panel images and prints
  about 9.3 in tall at 7.25 in wide.

Scripts written before September 2026 used supplementary numbers one higher (Fig S2–S6
for today's S1–S5).
