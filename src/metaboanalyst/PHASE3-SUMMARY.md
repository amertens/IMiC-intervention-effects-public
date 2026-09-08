# Phase 3 summary: secondary & tertiary milk pathway analysis

**Date:** 2026-07-22
**Spec:** `docs/superpowers/specs/2026-07-22-metaboanalystr-phase3-secondary-tertiary-design.md`

Extended the pipeline to the **secondary** and **tertiary** outcome groups across
**combined and stratified** arms, using the **Pathway module** (`pathora`, filter OFF), 
Trenton's rule. The Phase-0/1/2 machinery was generalized (`run_outcome_group`,
`build_supplementary_table`) rather than duplicated; `run-primary.R` is now a thin wrapper
over the general driver, and the primary nicotinate **7/5** result still reproduces
(regression, test-09/test-11/test-13).

## What ran

| group / arms | cells with results | skipped | significant pathways |
|---|---|---|---|
| secondary / combined | 0 | 7 | n/a |
| secondary / stratified | 0 | 26 | n/a |
| tertiary / combined | 8 | 6 | 6 |
| tertiary / stratified | 22 | 16 | 7 |

Outputs: `results/metaboanalyst/<group>_<arm>/` (per-cell `results.csv`/`membership.csv`,
`<group>_<arm>_all_cells.csv`, `_supplementary_table.csv`, `_significant.csv`,
`skipped_cells.csv`).

## Findings

- **Secondary produces no pathway results**: its features are HMOs and bioactives, which
  are not in SMPDB metabolic-pathway sets. Every cell is skipped with "too few mappable
  metabolites for enrichment". This is a genuine data/database limitation, not a bug.
- **Many tertiary cells skip for the same reason**: tertiary is dominated by lipids
  (triglycerides, phosphatidylcholines, ceramides, ...) that don't map to SMPDB pathways
  (even a 199-compound cell can leave < 3 mappable metabolites). The cells that DO map carry
  amino-acid / acylcarnitine / bile-acid signal.
- **Significant tertiary pathways (fdr_native < 0.05) are coherent and directional**: all
  in DOWNregulated arms, concentrated in amino-acid and fatty-acid-oxidation metabolism:
  - combined: Homocysteine Degradation, Oxidation of Branched-Chain Fatty Acids, β-Oxidation
    of Very-Long-Chain Fatty Acids, Carnitine Synthesis, Methionine Metabolism (Misame BEP /
    Elicit Nico, down).
  - stratified: Spermidine & Spermine Biosynthesis, Methionine, Glutamate, Glycine & Serine,
    Homocysteine Degradation (Misame BEP/BEP, IFA/BEP; Elicit Nico+Az., down).

## Robustness fix (important)

The first batch exposed a MetaboAnalystR 4.3.0 crash: when `current.msg` is seeded (as
`run_ora` leaves it) and a pathway cell hits the too-few-metabolites path, execution
proceeds into C code that **segfaults the entire process** (uncatchable). `run_pathway` now
un-seeds `current.msg` and wraps `CalculateOraScore`, converting the failure into a clean,
catchable "too few mappable metabolites" skip reason (commit `26c149f`, test-13).

## Out of scope / next

- **Exploratory / untargeted (NA outcome_group, ~99k mz/rt features) + proteomics:** not
  name-mappable → handled by the existing mummichog pipeline, not this one.
- **Blood compartments:** untargeted mz/rt; parked as a separate annotation project
  (`docs/superpowers/specs/2026-07-22-metaboanalystr-blood-annotation-followon.md`).
