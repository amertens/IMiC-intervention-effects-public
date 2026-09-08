# Re-deriving Table S2 (tertiary MSEA): findings

**Date:** 2026-07-22
**Goal:** reproduce the manuscript's Table S2 (tertiary metabolite-set enrichment) with a
scripted, reproducible pipeline and diff it against the published Q-values.

## Method discovery (confirmed)

Table S2 was produced with the **ORA / over-representation module + a name-matched reference
metabolome** (the same method as the primary analysis) **not** the pathway-analysis module our
Phase-3 first used (which was based on Trenton's verbal "tertiary = pathway analysis" description).

Evidence:
- Table S2's column is **"Enrichment Ratio"** (ORA), not "Impact" (pathway module).
- Its **Totals are restricted** (Glycine/Serine = 13, not the full SMPDB set ~50), i.e. a reference
  metabolome was applied.
- Re-running with ORA + a matched reference reproduces those restricted totals (Ammonia Recycling =
  7 exact, Glutamate = 8 exact, Glycine/Serine 11 vs 13).

Table S2 is also the **stratified (per-arm) analysis**: the same study/timepoint carries duplicate
pathway rows with different hit counts (different arms), collapsed under one study/timepoint label
(no arm column).

## Reproduction result

Tooling built: `R/build-reference.R` (`build_matched_reference()` automates Trenton's manual
MetaboAnalyst ID-conversion: 581 raw tertiary names → 106 SMPDB-matched names) and
`run-tertiary-msea.R` (ORA + reference, combined and stratified, signs the enrichment ratio by
direction). Outputs: `results/metaboanalyst/tertiary_msea/`.

**Structural reproduction is excellent:**
- **All 45 published FDR-significant pathway rows reproduced** by pathway × study × timepoint
  (stratified).
- **Total agreement 45/45** (|ours − published| ≤ 2); hits within ±1; directions match; the same
  pathway families (amino-acid/nitrogen + energy/fatty-acid-oxidation, all downregulated).

**But the Q-values do NOT match, ours are uniformly non-significant (FDR ≈ 1.0 vs published
~1e-6 to 1e-2), 0/45 significant.**

## Root cause (localized)

The gap is entirely in the ORA **background denominator**, visible in `Expected`:

| pathway (MISAME) | Total | published Expected | our Expected | ratio |
|---|---|---|---|---|
| Spermidine & Spermine Biosynth. | 5 | 0.0477 | 0.991 | ~21× |
| Glycine & Serine Metabolism | 13 | 0.155 | ~1.1 | ~7× |

`Expected = Total × (querySize / backgroundSize)`. Totals match, so **Trenton's effective background
was ~20× larger than our matched reference (106)**, his query was a tiny fraction of the background
(≈1%), ours is ≈20%. With a query that is 20% of the background, no pathway can be
over-represented after FDR, hence 0/45.

## Open question / what closes the gap

The exact background Trenton's MetaboAnalyst run used for Table S2 is unknown to us. Candidates:
- the full SMPDB/HMDB **library** as background (large denominator), but that usually also un-restricts
  the Totals, which Table S2 did not show, so this is a partial fit at best;
- a larger measured-metabolome reference (e.g. all measured metabolites incl. unmatched, or a
  differently-curated tertiary reference) than our 106-name matched set;
- a MetaboAnalyst version/setting difference in how `Expected` is computed when a reference is applied.

**This is also a methodological sensitivity worth the authors' attention:** the tertiary enrichment
*significance* hinges on the ORA background choice. Under a conservative measured-metabolome
background (ours), the tertiary pathway enrichment does not survive FDR; under a much larger
background (Table S2), it does. The pathway *identities and directions* are robust either way.

**To fully reconcile Table S2 numerically, we need Trenton's exact reference file + MetaboAnalyst ORA
settings (background definition).** The pathways, totals, hits, and directions are already
reproduced.
