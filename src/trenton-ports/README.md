# Trenton-ports: faithful R conversions of Trenton's Rmd analyses

Each script here is a **faithful, top-to-bottom R port of one of Trenton's `.Rmd`
notebooks** from `Dailey-Chwalibóg, Mertens et al. 2026 Science/2. Scripts/`. The
goal is fidelity: same input data, same filters, same cell definitions, same test
parameters, same order, *unless Trenton made a clear coding error*, which is
called out inline where it occurs.

The only thing that changes versus the Rmd is that the **manual metaboanalyst.ca
steps** (uploading each cell's compound list / the reference metabolome to the web
tool, then downloading the result) are replaced by the repo's scripted, offline
**MetaboAnalystR engine** (`src/metaboanalyst/R/run-ora.R`, `run-pathway.R`,
`harvest.R`, `build-reference.R`). That engine is validated to reproduce Trenton's
downloaded results **exactly** (see below), so the automation introduces no numeric
change.

## Fidelity gate (why we trust the automation)

The scripted ORA/pathway engine matches Trenton's downloaded MetaboAnalyst files
to the last digit on every cell for which he saved a `Download/`:

| Cell (downloaded ground truth) | module | agreement |
|---|---|---|
| Elicit · Up · 1 mo (Combined) | Pathway | 6/6 pathways, rel-diff 0 |
| Elicit · Up · 1 mo (Combined) | ORA | 60/60 pathways, raw-p & FDR Δ = 0 |
| Elicit · Up · 5 mo (Combined) | ORA | 60/60, Δ = 0 |
| Elicit · Down · 5 mo (Combined) | ORA | 9/9, Δ = 0 |
| Mumpta · Up · 1.5 mo (Combined) | ORA | 76/76, Δ = 0 |

(Validation drivers: `src/metaboanalyst/validate-against-trenton.R` and the
ground-truth `msea_ora_result.csv` files from the 2026-08-10 export.)

## Scripts

| Port | Source Rmd | Engine call | Output |
|---|---|---|---|
| `primary-msea-ora.R` | Primary Outcomes (MSEA ORA).Rmd | `run_ora` (msetora + reference) | `results/trenton-ports/primary_msea_ora/` |
| `primary-pathway-analysis.R` | Primary Outcomes (Pathway Analysis).Rmd | `run_pathway` (pathora) | `results/trenton-ports/primary_pathway/` |
| `exploratory-proteomics-genelevel.R` | Exploratory Outcomes (Proteomics).Rmd | clusterProfiler `enrichGO` | `results/trenton-ports/proteomics_genelevel/` |
| `exploratory-proteomics-uniprot.R` | Exploratory Outcomes (Proteomics - UniProt).Rmd | clusterProfiler `enricher` (UniProt TERM2GENE) | `results/trenton-ports/proteomics_uniprot/` |
| `compartment-tracking.R` | Compartment Tracking.Rmd | mummichog + his 25 ppm / name / feature-id linkage | `results/trenton-ports/compartment_tracking/` |

## Figure-code ports (recovered from the old `imicPaper*.Rmd`, 2025)

The submitted **figure-generating** code (which the 2026 analysis-only notebooks
lacked) was found in Trenton's older scripts and ported here:

| Panel | Source Rmd | Port | Status |
|---|---|---|---|
| **Fig 6C** proteome-GO | `imicPaperProteomics.Rmd` (`★Untargeted Proteomics Figure`) | `figure-6c-proteome-go-submitted.R` | Reproduces the **submitted** figure (global-universe gene-level enrichGO). **Fixes** his copy-paste bug (Mumpta-up filtered MISAME). This is the global-universe result that makes up-regulated terms survive FDR; the corrected per-cell/UniProt methods overturn it. |
| **Fig 6A** untargeted MSEA | `imicPaperUntargetedMetabolomics.Rmd` (`★Untargeted Metabolites Figure`) | `figure-6a-untargeted-msea-submitted.R` | Exact figure code + reproducible cell logic + `run_ora` wiring. **Gated** on Trenton's untargeted annotation keys (not in-repo); runs 1:1 once supplied. |
| **Fig 5C** fatty-acid enrichment | `imicPaperTriglycerides.Rmd` (`★ Fatty Acids Figure`) | already faithfully reproduced in-repo at `figure-scripts/manuscript_figures/fig5C-triglyceride.R` (needs the Biocrates Quant500 file, which IS present). No duplicate port needed. |
| **Fig 6B** mummichog | `imicUntargetedMetabolomicsMummichog.Rmd` | mummichog **run** code present (conda, `-u 10 -n human_mfn -c 0.05`); the pathway-volcano **plot** code is incomplete in his file ("START OVER"), still pending. |

Run any port from the repo root, e.g.:

```bash
Rscript src/trenton-ports/primary-msea-ora.R
Rscript src/trenton-ports/figure-6c-proteome-go-submitted.R
```
