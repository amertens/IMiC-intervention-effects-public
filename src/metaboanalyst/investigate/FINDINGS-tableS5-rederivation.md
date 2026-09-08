# Re-deriving Table S5 (untargeted milk Mummichog): findings

**Date:** 2026-07-23
**Script:** `src/metaboanalyst/run-milk-mummichog-s5.R`
**Output:** `results/metaboanalyst/mummichog_s5/milk_mummichog_tableS5.csv`

## Result: 196 / 198 published rows reproduced (99%)

Matched on pathway × study × timepoint × contrast × regulation × ionization mode.
This is by far the closest of the three supplementary-table reproductions
(cf. Table S2 structural-only, Table S6 not reproduced).

| check | result |
|---|---|
| published rows matched | **196 / 198** |
| overlap size within ±2 | 158 / 196 (81%) |
| our `pathway_size` > published | 173 / 196 (88%) |
| median `pathway_size` ratio (ours / published) | **1.39** |
| median p-value ratio (ours / published) | **2.17** |

## Inputs and conventions established

- **Input:** `results/adjusted_intervention_effects_res_untargeted_metabolomics_clean_ATE.RDS`
  (STRATIFIED arms), its contrasts (MISAME BEP/BEP, BEP/IFA, IFA/BEP; Vital
  BEP+ExBf, BEP+ExBf+AZT; Elicit Az., Nico, Nico+Az.) match Table S5's Contrast column.
- **m/z–RT:** `data/additional datasets/IMiC_alignment.csv`, per study
  (`*_MISAME3` for MISAME; `*_CHILD_ELICIT_VITAL` for Vital/Elicit).
- **Enrichment Ratio = overlap_size / pathway_size**, signed by direction: verified
  against published rows (31/68 = 0.456; 30/77 = 0.390; 9/12 = 0.75). This is *not*
  an observed/expected ratio.
- **FDR:** Benjamini–Hochberg across pathways within each run.
- **Mummichog config:** human_mfn, 10 ppm, p-cutoff 0.05, per ionization mode.
- **Direction** (CORRECTED 2026-08: see note below): each directional run submits
  **only the in-direction features**, with their real p-values, so Mummichog's
  background is the in-direction feature set. This matches Trenton's own directional
  runs (his `_neg` output file contains only the down-regulated features).

  > **Note, method superseded.** An earlier version of this note (and
  > `18-directional-mummichog.R`) described submitting the FULL feature list each run
  > and neutralising the opposite direction to p = 1. That approach left the
  > opposite-direction features in the background and roughly doubled it, shifting
  > every permutation p-value even though the significant set was unchanged. The
  > current `run-milk-mummichog-s5.R` submits in-direction features only (see its
  > header and `run_cell()`); that is the authoritative Table S5 method. Cite the
  > script, not the older full-list/p=1 description.

## Mixed arm frameworks (explains the last 8 rows)

The first pass matched 190/198. All 8 misses were **Mumta-LW 1.5 mo, contrast "BEP"**, 
a *combined*-arm contrast, while every other row is stratified. This matches the
Methods statement that "combined-arm analyses [were] used for time points preceding
antibiotic administration". Adding that one combined-arm cell
(`COMBINED_CELLS` in the script) took the match to **196/198**.

So Table S5 legitimately mixes combined- and stratified-arm results, keyed to whether
the timepoint precedes azithromycin.

## Remaining gap: one upstream feature filter

> **Caveat:** the reproduction statistics in this section were computed under the
> earlier full-list/p=1 background. The corrected in-direction-only method (above)
> roughly halves the per-direction background, so it likely shrinks this
> `pathway_size` gap; re-run the comparison before treating the numbers below as
> current.

The residual difference is systematic and has a single cause: **our submitted feature
list is larger than Trenton's** (median `pathway_size` ratio 1.39). A larger mapped
background inflates `pathway_size` and makes p-values less extreme, which explains
both remaining discrepancies with one mechanism (median p ratio 2.17).

**Open question for Trenton:** was a feature filter applied before submitting to
Mummichog (a detection/quality threshold, or a restricted subset of the alignment
key)? Our MISAME 3–4 mo negative-mode run submits ~27,700 features and yields
`pathway_size` 94 where the published table shows 68.

Note this does *not* affect which pathways are found, the pathway identities,
directions, contrasts and ranking reproduce; only the exact p/FDR values shift.
