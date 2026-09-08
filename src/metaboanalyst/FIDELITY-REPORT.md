# MetaboAnalystR pipeline: Phase 0/1 fidelity report

**Date:** 2026-07-22
**Scope:** Reproducing Trenton's manual metaboanalyst.ca enrichment/pathway analysis
for the **primary** outcome group, with a scripted MetaboAnalystR pipeline
(`src/metaboanalyst/`). MetaboAnalystR 4.3.0, R 4.4.2, fully local/offline.

## Headline

The scripted pipeline reproduces Trenton's one downloaded MetaboAnalyst result
**exactly**, and reproduces the structure he described for the enrichment
(ORA) analysis. Every step is now code, no manual clicking on metaboanalyst.ca.

## Hard numeric gate: the golden pathway cell (PASSES)

The only cell for which Trenton has a downloaded ground-truth file is the
**Pathway Analysis** of *Elicit · Upregulated · 1 month (Combined)*
(`.../Pathway Analysis/Download/pathway_results.csv`). Our engine vs his file,
all six shared pathways:

| Pathway | Our raw p | Trenton raw p | rel. diff |
|---|---|---|---|
| Nicotinate and Nicotinamide Metabolism | 1.2011e-05 | 1.2011e-05 | 0 |
| Androstenedione Metabolism | 0.1295 | 0.1295 | 0 |
| Androgen and Estrogen Metabolism | 0.1609 | 0.1609 | 0 |
| Steroidogenesis | 0.2257 | 0.2257 | 0 |
| Glycine and Serine Metabolism | 0.2634 | 0.2634 | 0 |
| Tryptophan Metabolism | 0.2862 | 0.2862 | 0 |

Reproduced with `run_pathway()` (SMPDB / Homo sapiens, hypergeometric,
relative-betweenness topology, metabolome filter OFF), matching the exact call
sequence in the downloaded `Rhistory.R`. (`src/metaboanalyst/validate-against-trenton.R`,
`test-10`.)

## Primary enrichment (ORA): the actual deliverable

Trenton's rule: primary uses **ORA (over-representation)**. The driver
`run-primary.R` runs all **11 non-empty primary cells** (Misame×BEP, Vital×BEP,
Elicit×Nico; up/down), applying the reference metabolome. Results are written to
`results/metaboanalyst/primary/` (per-cell `ora_results.csv` +
`ora_membership.csv`, plus a combined `primary_enrichment_all_cells.csv`).

- **All 11 cells succeed; none skipped** (including the 1-compound Elicit·5mo·down cell).
- **Elicit · Upregulated · 1 mo · Nicotinate & Nicotinamide Metabolism:**
  total = **7**, hits = **5**, raw p = 3.67e-04, FDR = 0.036.
  The **total = 7 / hits = 5** structure matches exactly what Trenton described
  on the call ("in our subset there's seven ... we have five hits"). This is the
  strongest available corroboration of the ORA path, since no downloaded ORA
  result file exists to check the p-value against.

### Reference-metabolome fidelity (important)

The ORA background universe ("reference metabolome") materially changes the
numbers, and **which name list you use matters**:

| Reference used | Nicotinate total | hits | raw p |
|---|---|---|---|
| Trenton's curated file (35 MetaboAnalyst-matched names) | 7 | 5 | 3.67e-04 |
| Our data-derived labels (38 raw names from the RDS) | 4 | 3 | 7.33e-04 |
| No reference (full SMPDB library) | 35 | 5 | 2.28e-07 |

Trenton's curated file holds the **ID-conversion-matched** compound names
(e.g. `L-Tryptophan`, `alpha-Tocopherol`), which match the SMPDB metabolite-set
library far better than the raw data labels (`tryptophan`, `Alpha-tocopherol`).
The pipeline therefore uses his curated file
(`PRIMARY_CONFIG$reference_path`) to reproduce his structure. **Auto-generating
matched-name references for the other modalities is a Phase-3 task** (it means
running each outcome group's measured compounds through MetaboAnalyst's name
mapping, the scripted equivalent of Trenton's manual ID-conversion step).

## Reference-metabolome question: resolved

The apparent contradiction (his pathway `Rhistory.R` shows the filter OFF, yet
he described uploading a reference) is resolved (`investigate/FINDINGS-reference-metabolome.md`):

- **Pathway module: filter OFF.** Filter ON with the KEGG reference is not just
  different but *broken* on the SMPDB library, `CalculateOraScore` filters
  against KEGG IDs while SMPDB set members are HMDB IDs (zero overlap), so it
  collapses to the package's "too few sets" error. Filter OFF reproduces his
  downloaded numbers.
- **ORA module: reference ON**, via `Setup.HMDBReferenceMetabolome` +
  `SetMetabolomeFilter(TRUE)`, matching his web procedure and reproducing the
  7/5 structure.

## MetaboAnalystR 4.3.0 bugs found and worked around

1. **`InitDataObjects` self-referential default arg**, the 4th `dpi` argument
   must be passed explicitly (e.g. `150`), else "promise already under evaluation".
2. **`CalculateHyperScore` crashes locally**, it calls two reporting-only
   helpers (`ExportOraMembershipJson`, `PlotORAMembership`) with an NA `mSetObj`
   off the public web, crashing *after* results are computed. Patched by
   no-op'ing those two side effects (`run-ora.R`; verified numerically identical).
3. **`AddErrMsg` reads an uninitialized global** off-web, so error paths would
   hard-crash; seeded `current.msg` defensively.

None affect any numeric result, they are all environment/reporting bugs in the
package's off-web code path.

## Known limitation & request to Trenton

Our **hard numeric validation rests on a single downloaded cell** (the pathway
result above). The primary **ORA** cells (the actual deliverable) have only
query-list `.xlsx` files saved, no downloaded result folders, so they are
validated by structure (total/hits) and internal consistency, not against
downloaded numbers.

**Request:** please share the downloaded `Download.zip` (or the
`ora_results` / membership files) for a few primary **enrichment** cells, 
ideally *Elicit · Upregulated · 1 mo* and *5 mo*, so we can lock the ORA path
to exact numbers the same way we did for the pathway cell.

## What runs, and how

```
# one-time environment setup
Rscript src/metaboanalyst/00-setup-metaboanalystr.R
# reproduce the primary enrichment analysis end to end
Rscript src/metaboanalyst/run-primary.R          # -> results/metaboanalyst/primary/
# full test suite (golden gate + wrappers + build + validate)
Rscript -e "testthat::test_dir('src/metaboanalyst/tests')"
```

## Next phases (separate plans)

- **Phase 2: critical evaluation:** native FDR vs Trenton's own BH on raw p;
  ORA-vs-pathway choice; running both modules per cell; auto matched-name
  reference generation.
- **Phase 3: expansion:** secondary/tertiary/exploratory outcome groups and
  other milk/blood modalities via config rows (`build_cells` + a config like
  `config-primary.R`), plus the stratified-arms cells (Misame BEP/BEP, BEP/IFA, ...).
