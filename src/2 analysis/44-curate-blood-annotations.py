#!/usr/bin/env python
"""
44-curate-blood-annotations.py

Curate the putative blood-feature annotations (from 23-annotate-fdr-features.R,
which reads mummichog userInput_to_EmpiricalCompounds.tsv) into a clean,
shareable table of consistently-affected named features.

Input : results/fdr_sig_putative_annotation.csv   (FDR-sig blood features + mummichog empirical-compound names)
        results/four_compartment_bep_shortlist.csv (adds milk direction, by feature id)
Output: results/blood_consistent_annotated_features.csv  (one row per feature)
        results/blood_consistent_annotated_features_consistent.csv (cross-compartment subset)

Confidence tiers:
  high      - Sapient-provided single name (sapient_name filled)  e.g. Octenoylcarnitine
  medium    - single unambiguous mummichog candidate (n_candidates<=1, real name)
  ambiguous - multiple mummichog candidates for the same m/z (isobaric set)
  unresolved- junk/partial annotation (e.g. bare "Acid")

"Consistent" = FDR-up in infant blood AND in >=1 maternal compartment (plasma or maternal VAMS).
NOTE: milk annotations are computed separately (Trenton's mummichog outputs); milk direction
here is carried from the four-compartment screen where available.
"""
import csv, os, io, sys
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
RES  = os.path.join(ROOT, "results")

JUNK = {"acid", "", "na", "none", "- none -", "-"}

def fnum(x):
    try: return float(x)
    except Exception: return None

def clean_candidates(row):
    # Prefer the pipe-separated mummichog candidate list; if it is empty, fall back to
    # the Sapient name, then to best_annotation. Finally drop junk/placeholder tokens.
    raw = (row.get("putative_candidates") or "").strip()
    cands = [c.strip() for c in raw.split("|") if c.strip()] if raw else []
    if not cands:
        sapient = (row.get("sapient_name") or "").strip()
        best_ann = (row.get("best_annotation") or "").strip()
        cands = [sapient] if sapient else ([best_ann] if best_ann else [])
    cands = [c for c in cands if c.lower() not in JUNK]
    return cands

def main():
    ann_path = os.path.join(RES, "fdr_sig_putative_annotation.csv")
    rows = list(csv.DictReader(open(ann_path, encoding="utf-8")))

    # milk direction from four-compartment shortlist (keyed by feature id)
    milk = {}
    short_path = os.path.join(RES, "four_compartment_bep_shortlist.csv")
    if os.path.exists(short_path):
        for r in csv.DictReader(open(short_path, encoding="utf-8")):
            milk[r["feature"]] = r

    DS = {"VamsPostnatalInfant": "infant", "MaternalPlasma": "plasma",
          "VamsPostnatalMaternal": "matVAMS"}

    feats = {}
    for r in rows:
        fid = r["feature"]
        est = fnum(r["est"]); q = fnum(r["q"])
        comp = DS.get(r["dataset"], r["dataset"])
        cands = clean_candidates(r)
        sn = (r.get("sapient_name") or "").strip()
        f = feats.setdefault(fid, {
            "feature": fid, "mz": r.get("mz", ""), "mode": r.get("mode", ""),
            "sapient_name": sn, "candidates": set(), "hits": {}})
        if sn: f["sapient_name"] = sn
        for c in cands: f["candidates"].add(c)
        # record per-compartment best (min q) up-hit
        if est is not None and est > 0:
            prev = f["hits"].get(comp)
            if prev is None or (q is not None and q < prev[1]):
                f["hits"][comp] = (est, q if q is not None else 9.9,
                                   r.get("visit", ""))

    out = []
    for fid, f in feats.items():
        cands = sorted(f["candidates"])
        sapient = f["sapient_name"]
        # confidence tier: a Sapient name is best; a single mummichog candidate is
        # unambiguous; several candidates share the m/z (isobaric); none = unresolved.
        if sapient:
            conf, name = "high", sapient
        elif len(cands) == 1:
            conf, name = "medium", cands[0]
        elif len(cands) > 1:
            conf, name = "ambiguous", cands[0]
        else:
            conf, name = "unresolved", "(unresolved)"
        comps = f["hits"]
        # "consistent" requires an FDR-strict (q<0.05) UP hit in infant blood AND in at
        # least one maternal compartment (plasma or maternal VAMS).
        strict = {c for c, (e, q, v) in comps.items() if q is not None and q < 0.05}
        consistent = ("infant" in strict) and bool(strict & {"plasma", "matVAMS"})
        m = milk.get(fid, {})
        milk_up = ""
        if m:
            milk_est = fnum(m.get("milk_est"))
            if milk_est is not None:
                milk_up = "up" if milk_est > 0 else "down"
        # render one compartment's effect cell as "+est (q=..., visit)", or "" if absent
        def fmt(c):
            if c in comps:
                e, q, v = comps[c]
                return f"{e:+.2f} (q={q:.1e}, {v})"
            return ""
        out.append({
            "feature": fid, "mz": f["mz"], "mode": f["mode"],
            "name": name, "confidence": conf,
            "n_candidates": len(cands),
            "all_candidates": " | ".join(cands),
            "infant": fmt("infant"), "plasma": fmt("plasma"),
            "matVAMS": fmt("matVAMS"), "milk_direction": milk_up,
            "compartments_up_FDR": ",".join(sorted(strict)),
            "cross_compartment_consistent": "yes" if consistent else "",
        })

    # sort: cross-compartment-consistent rows first, then by confidence tier, then by
    # how many compartments the feature is FDR-up in (more first). Each key ascends, so
    # booleans/tiers are negated where "more desirable" should come first.
    corder = {"high": 0, "medium": 1, "ambiguous": 2, "unresolved": 3}
    def sort_key(r):
        not_consistent = r["cross_compartment_consistent"] != "yes"   # False (0) sorts first
        tier = corder[r["confidence"]]
        n_up = len(r["compartments_up_FDR"].split(",")) if r["compartments_up_FDR"] else 0
        return (not_consistent, tier, -n_up)
    out.sort(key=sort_key)

    cols = ["feature", "name", "confidence", "n_candidates", "all_candidates",
            "mz", "mode", "infant", "plasma", "matVAMS", "milk_direction",
            "compartments_up_FDR", "cross_compartment_consistent"]
    op = os.path.join(RES, "blood_consistent_annotated_features.csv")
    with open(op, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols); w.writeheader()
        for r in out: w.writerow(r)

    cons = [r for r in out if r["cross_compartment_consistent"] == "yes"]
    op2 = os.path.join(RES, "blood_consistent_annotated_features_consistent.csv")
    with open(op2, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols); w.writeheader()
        for r in cons: w.writerow(r)

    print(f"wrote {op}  ({len(out)} features)")
    print(f"wrote {op2}  ({len(cons)} cross-compartment-consistent)")
    print("\n== cross-compartment-consistent named features (infant + maternal, FDR<0.05 up) ==")
    for r in cons:
        print(f"  [{r['confidence']:9}] {r['name']:22} mz={r['mz']:>9} "
              f"| infant {r['infant']:22} | plasma {r['plasma']:22} | matVAMS {r['matVAMS']:20}"
              + (f" | alt: {r['all_candidates']}" if r['confidence'] == 'ambiguous' else ""))
    print("\n== confidence breakdown (all FDR-sig blood features) ==")
    from collections import Counter
    for k, v in Counter(r["confidence"] for r in out).most_common():
        print(f"  {k:11} {v}")

if __name__ == "__main__":
    main()
