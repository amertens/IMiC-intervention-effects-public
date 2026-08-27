# Pathway / enrichment pipeline (`src/metaboanalyst/`)

A scripted, push-button **MetaboAnalystR** pipeline that reproduces the milk metabolome/proteome
enrichment analyses for the IMiC intervention-effects paper — replacing the earlier manual
metaboanalyst.ca point-and-click workflow with code that runs the same analyses identically every
time.

It mirrors the original methodological rule:

| Outcome group | Module | MetaboAnalyst function |
|---|---|---|
| **Primary** (targeted B-vitamins / micronutrients) | Over-Representation Analysis (ORA) | `msetora`, with a name-matched reference metabolome |
| **Secondary / tertiary** (HMOs, bioactives, targeted lipidome) | Pathway / impact analysis | `pathora` |
| **Untargeted milk metabolome** | Mummichog | (directional, per ion mode) |
| **Untargeted milk proteome** | GO Biological Process | `clusterProfiler::enrichGO` |

## Architecture

```
 driver scripts                 core engine (R/)                        outputs
 --------------                 ----------------                        -------
 run-primary.R            ┐
 run-primary-pathway-...  │     run-outcome-group.R  ── build-cells.R   results/metaboanalyst/
 run-secondary-tertiary.R ├──►  (dispatch a group    ── run-ora.R          <group>_<arm>/
 run-tertiary-msea.R      │      through a module)   ── run-pathway.R         <cell>/results.csv
 run-milk-mummichog-s5.R  ┘                          ── harvest.R             <cell>/membership.csv
 run-proteomics-go.R  (standalone)                                           <group>_<arm>_all_cells.csv
                                                                             skipped_cells.csv
 build-supplementary-table.R  ── consolidates a group's cells into one tidy table
```

`run_outcome_group()` is the heart: it builds one "cell" per study × visit × contrast × direction
(`build-cells.R`), runs each cell through the chosen module (`run-ora.R` / `run-pathway.R`), tags
and harvests the results (`harvest.R`), and writes per-cell + consolidated CSVs. Cells that fail
legitimately (a pathway needs ≥ 3 mappable metabolites) are caught and logged to
`skipped_cells.csv`, never fatal.

## Scripts, at a glance

| Script | What it produces | Manuscript output |
|---|---|---|
| `run-primary.R` | Primary ORA, **combined** and **stratified** arms | Table S1 (primary MSEA), Fig 3 enrichment |
| `run-primary-pathway-compare.R` | Combined-arm primary via **pathway/impact** (for ORA-vs-pathway comparison) | — (comparison only) |
| `run-secondary-tertiary.R` | Secondary + tertiary pathway analysis, both arm sets | — |
| `run-tertiary-msea.R` | Tertiary targeted-metabolome ORA/MSEA, both arm sets | Table S2 |
| `run-milk-mummichog-s5.R` | Directional Mummichog of the untargeted milk metabolome | Table S5 |
| `run-proteomics-go.R` | GO-BP enrichment of the untargeted milk proteome | Table S6 |
| `build-supplementary-table.R` | Consolidates a group×arm's cells into one supplementary table | — |
| `validate-against-trenton.R` | Diffs our numbers against Trenton's downloaded ground-truth cell | — |
| `R/build-cells.R` | Splits a result frame into per-cell query + reference lists | — |
| `R/build-reference.R` | Builds the name-matched reference metabolome (MetaboAnalyst name-mapping) | — |
| `R/run-ora.R` | `msetora` wrapper (+ documented 4.3.0 bug workarounds) | — |
| `R/run-pathway.R` | `pathora` wrapper | — |
| `R/harvest.R` | Pulls the results + membership tables out of an `mSet` | — |
| `R/config-primary.R`, `R/config-pathway-groups.R`, `R/label-map.R` | Per-group configuration + label maps | — |

## How to run (from the repo root)

```bash
Rscript "src/metaboanalyst/run-primary.R"                 # primary ORA (combined + stratified)
Rscript "src/metaboanalyst/run-primary-pathway-compare.R" # combined-arm pathway variant
Rscript "src/metaboanalyst/run-secondary-tertiary.R"      # secondary + tertiary pathway
Rscript "src/metaboanalyst/run-tertiary-msea.R"           # tertiary MSEA (Table S2)
Rscript "src/metaboanalyst/run-milk-mummichog-s5.R"       # milk Mummichog (Table S5)  [needs conda mummichog env]
Rscript "src/metaboanalyst/run-proteomics-go.R"           # proteome GO (Table S6)
```

Outputs land in `results/metaboanalyst/<group>_<arm>/` (gitignored). `build-supplementary-table.R`
turns any of those into a clean table.

`results/metaboanalyst/untargeted_msea/` and `results/metaboanalyst/triglyceride_fa/` are the reproducible CSV exports of manuscript Table S4 (untargeted MSEA) and Table S3 (triglyceride fatty-acid composition) for the online supplement's §9.

Inputs are the combined-arm / stratified-arm intervention-effect result RDS
(`trenton scripts/1. Data/combined_intervention_effects_results_{combined,stratified}_arms.RDS`).

## Reproductions verified (so you can trust the swap)

- **Primary ORA**, ELICIT · up · 1 mo: Nicotinate & Nicotinamide Metabolism **total 7 / hits 5**,
  raw p = 3.74 × 10⁻⁴ — matches the manual run.
- **Primary pathway/impact**, same cell: raw p = **1.2011 × 10⁻⁵**, impact 0.205 — matches Trenton's
  downloaded pathway CSV exactly. (The pathway module uses the full SMPDB background — total 32 — so
  it is more significant than ORA's restricted-reference background — total 7. This is the
  ORA-vs-pathway difference to discuss, not a discrepancy.)

## MetaboAnalystR gotchas (all patched in-session, none change the numbers)

The messy parts are isolated in `R/run-ora.R` / `R/run-pathway.R`, each with a full explanatory
header. In brief (MetaboAnalystR **4.3.0**):
- `CalculateHyperScore()` ends in two off-public-web reporting side effects
  (`ExportOraMembershipJson`, `PlotORAMembership`) that crash after the numbers are already
  computed → replaced with no-ops.
- `AddErrMsg()` reads an uninitialized `current.msg`; a seeded value lets the too-few-metabolites
  path **segfault**, so we un-seed it right before scoring to force a clean, catchable R error.
- `InitDataObjects()` needs an explicit `dpi` argument (self-referential default bug).

## Environment

- MetaboAnalystR **4.3.0** (installed via `00-setup-metaboanalystr.R`; needs `qs2`, Rtools44).
- `Rscript` at `C:/Program Files/R/R-4.4.2/bin/Rscript.exe` (not on PATH).
- Milk Mummichog additionally needs the `mummichog` conda env.

## Related docs

- `FIDELITY-REPORT.md` — cell-by-cell fidelity of the replication.
- `PHASE3-SUMMARY.md` — secondary/tertiary milk pathway summary.
- `investigate/FINDINGS-*.md` — deep dives (reference metabolome, Table S2/S5 re-derivation).
