# Pathway and enrichment analyses (`src/metaboanalyst/`)

These scripts run the paper's metabolite pathway and enrichment analyses with
MetaboAnalystR 4.3.0, replacing the point-and-click metaboanalyst.ca web workflow
used during the analysis with code that gives the same result on every run. The
same folder holds the Mummichog analysis of the untargeted milk metabolome, the
proteome Gene Ontology table, and the cross-compartment pathway analysis.

## Runners

Output paths are under `results/metaboanalyst/` unless a full path is given.

| Runner | Method | Exhibit | Output |
|---|---|---|---|
| `run-primary-pathway-local.R` | KEGG pathway analysis (hypergeometric test, relative-betweenness topology) of FDR-significant primary metabolites, one cell per study × time point | Fig 3B, Table S2 | `primary_pathway_local/primary_pathway_all_cells.csv` |
| `run-tertiary-msea-dual.R` | Over-representation analysis (ORA) against SMPDB metabolite sets; separate metabolite and lipid passes, 1,268-name reference metabolome | Fig 5B, Table S3 | `tertiary_msea/tertiary_msea_dual.csv` |
| `run-untargeted-msea.R` | ORA against SMPDB metabolite sets, nominal P < 0.05 foreground, 615-name reference metabolome | Fig 6A, Table S5 | `untargeted_msea/untargeted_msea_combined.csv` (Fig 6A), `untargeted_msea/untargeted_msea_combined_fdr_sig.csv` (Table S5) |
| `run-milk-mummichog.R` | Mummichog (human_mfn network, 10 ppm), directional, per ionization mode | Fig 6B, Table S6 | `mummichog/milk_mummichog_pathways.csv` |
| `run-proteomics-go.R` | Gene Ontology ORA of the milk proteome with `clusterProfiler::enricher` and a UniProt-to-GO map (computed by `src/2 analysis/55-proteomics-go-uniprot.R`; this runner builds the table) | Table S7 | `proteomics_go/proteomics_go_pathways.csv` |
| `build-compartment-query-list.R` | Links FDR-significant features across maternal blood, milk and infant blood (MISAME-III) and builds the pathway query list | Table S11 (query) | `results/compartment_tracking/linked_upregulated.csv`, `results/compartment_tracking/pathway_compound_list.csv` |
| `run-compartment-pathway.R` | KEGG pathway analysis of the cross-compartment query list | Fig 6D, Table S11 | `results/compartment_tracking/linked_crosscompartment.csv`, `results/compartment_tracking/metabolite_pathways.csv`, `results/tables/table_s11_compartment_pathway.csv` |
| `build-reference-metabolome.R` | Maps the study reference metabolome (1,593 HMDB IDs) to MetaboAnalyst compound names | Background for Fig 6A, Table S5 | `reference/reference_metabolome_1593_matched_names.txt` |
| `00-setup-metaboanalystr.R` | Installs MetaboAnalystR and its dependencies | n/a | `env/sessionInfo.txt` |

## Helper files (`R/`)

- `run-pathway.R`: `run_pathway()` (MetaboAnalystR pathway module) and `harvest_pathway()` (pathway table plus compound-to-pathway hits).
- `run-ora.R`: `run_ora()` (MetaboAnalystR ORA module), including the session patch for a MetaboAnalystR 4.3.0 crash that happens after the results are computed.
- `harvest.R`: reads ORA/pathway results into tables; `apply_pathway_size_floor()` for the reporting size filter.
- `build-cells.R`: splits an intervention-effects table into per-cell query lists.
- `label-map.R`: maps targeted-assay labels to MetaboAnalyst compound names.
- `lipid-name-map.R`: converts Quant 500 lipid names to LIPID MAPS abbreviations for the lipid pass.
- `build-reference.R`: scripted compound-name matching for building a reference metabolome.
- `compartment-tracking.R`: cross-compartment linkage (`ct_compartment_pairs()`) and name-based grouping (`ct_name_track()`).

`reference/` holds the reference metabolomes (`refMetabolomeForQER.csv`, the
1,593-ID HMDB list and its 615 matched names) and the KEGG pathway-name map.

## Inputs

- On request (feature-level, too large to ship): `results/combined_intervention_effects_results_{combined,stratified}_arms.RDS`
  (primary and tertiary runners). `run-milk-mummichog.R` also needs the untargeted
  milk effect estimates and `data/additional datasets/IMiC_alignment.csv`, which are
  likewise not shipped.
- Shipped: `results/milk_nominal_putative_annotation.csv` (untargeted ORA),
  `results/proteomics_go_uniprot.csv` (Table S7), `results/supplement_status_fdr_features.csv`
  (cross-compartment analysis), and the files in `reference/`. The untargeted ORA,
  proteome table, cross-compartment analysis and reference-metabolome build run from
  shipped files alone.

## How to run

From the repository root:

```bash
Rscript src/metaboanalyst/00-setup-metaboanalystr.R        # once
Rscript src/metaboanalyst/run-primary-pathway-local.R      # Fig 3B, Table S2   [needs on-request file]
Rscript src/metaboanalyst/run-tertiary-msea-dual.R         # Fig 5B, Table S3   [needs on-request file]
Rscript src/metaboanalyst/run-untargeted-msea.R            # Fig 6A, Table S5
Rscript src/metaboanalyst/run-milk-mummichog.R             # Fig 6B, Table S6   [needs on-request files and the mummichog conda env]
Rscript src/metaboanalyst/run-proteomics-go.R              # Table S7
Rscript src/metaboanalyst/build-compartment-query-list.R   # Table S11 query list
Rscript src/metaboanalyst/run-compartment-pathway.R        # Fig 6D data, Table S11
```

`build-compartment-query-list.R` must run before `run-compartment-pathway.R`.
`run-milk-mummichog.R` calls the `mummichog` command-line tool in a conda
environment named `mummichog`; set `IMIC_CONDA_CMD` if `conda` is not on the PATH.

## Library downloads

MetaboAnalystR downloads its compound, metabolite-set and pathway libraries from
metaboanalyst.ca when a run starts, so the MetaboAnalystR runners need internet
access; the analyses themselves run locally. Downloads are cached for the rest of
the R session. If the libraries on metaboanalyst.ca change, results can differ
slightly from the shipped outputs.

## Validation

The scripted engines were checked against the metaboanalyst.ca results saved for
the ELICIT 1-month up-regulated primary cell. The pathway engine reproduced the saved
web pathway results exactly for all six pathways in that cell (e.g. nicotinate and
nicotinamide metabolism, raw P = 1.2011 × 10⁻⁵, FDR = 1.19 × 10⁻³), and the ORA
engine reproduced the set size and hit count reported for that pathway in the web
ORA (7 and 5); no web ORA P-values were saved for that cell. The saved web results
are not included in this repository.
