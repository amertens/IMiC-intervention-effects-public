# Trenton's manuscript-analysis scripts

Trenton Dailey-Chwalibóg's original R Markdown scripts for the manuscript
enrichment/pathway analyses, copied verbatim from his Science submission package
(`Dailey-Chwalibóg, Mertens et al. 2026 Science/2. Scripts/`). Kept here so his
analyses are reproducible from this repo.

| File | Produces | Notes |
|---|---|---|
| `Primary Outcomes (MSEA ORA).Rmd` | Primary over-representation (up/down), SMPDB, 35-nutrient reference | Builds the reference metabolome, exports queries; enrichment run on the MetaboAnalyst web tool |
| `Primary Outcomes (Pathway Analysis).Rmd` | **Fig 3B** KEGG pathway-impact plot | Query = `sigFDR==1`, pos+neg combined; ellipse+labels for pathways significant in ≥2 cells (`pathway_count >= 2`). Reads `pathway_results.csv` exports (see `results/metaboanalyst/primary_pathway_trenton/`) |
| `Exploratory Outcomes (Proteomics - UniProt).Rmd` | Fig 6C proteome GO (UniProt-native, each protein once) | Analysis only, no figure code; the 6C figure was to be handed over separately |
| `Exploratory Outcomes (Proteomics) [OLD gene-level].Rmd` | superseded gene-level GSEA | Reference only, this is the inflated gene-level test that was replaced |
| `Compartment Tracking.Rmd` | Cross-compartment transfer (Fig 6D / Table S8; the old Fig S7 was cut 2026-08-26) | 87 plotting calls, his original cross-compartment figure code; our `fig6D-crosscompartment.R` was built from his DATA and should be reconciled against this |

## What is NOT here
His original figure scripts for the **targeted (Fig 5, 5B, 5C), untargeted MSEA (6A),
and mummichog (6B)** panels are not in the submission package; that package ships
*our* adapted scripts (its `Andrew scripts/` folder). Those panels are currently
matched from the submitted images, not reproduced from his source.

## How they run
These are MetaboAnalyst-web-driven: they prepare inputs, the enrichment/pathway
steps are run on metaboanalyst.ca, and the exported result CSVs are read back in.
They depend on his `1. Data/` RDS files and the exported result folders.
