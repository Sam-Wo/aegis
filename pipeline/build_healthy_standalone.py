#!/usr/bin/env python3
"""
Build a HEALTHY-ONLY standalone Aegis HTML (tumor groups stripped from the
embedded data). The app auto-hides all tumor UI when meta.tumors is empty.
Output: app/aegis_healthy_standalone.html
"""
import json
from pathlib import Path

AEGIS = Path(__file__).resolve().parent.parent
d = json.loads((AEGIS / "data" / "aegis_data.json").read_text(encoding="utf-8"))

groups = d["groups"]
keep = [i for i, g in enumerate(groups) if g.get("kind", "normal") != "tumor"]
remap = {old: new for new, old in enumerate(keep)}
newgroups = [groups[i] for i in keep]

newexpr = {}
for gene, entries in d["expr"].items():
    ne = [[remap[gi], me, pp] for gi, me, pp in entries if gi in remap]
    if ne:
        newexpr[gene] = ne
genes = sorted(newexpr.keys())
ann = {g: d["ann"][g] for g in genes if g in d.get("ann", {})}

meta = dict(d["meta"])
meta["n_groups"] = len(newgroups)
meta["n_genes"] = len(genes)
meta["total_cells"] = sum(g["n"] for g in newgroups)
meta["tumors"] = []
# meta["tissues"] already lists only normal tissues

out = dict(meta=meta, groups=newgroups, genes=genes, expr=newexpr, ann=ann)

tpl = (AEGIS / "app" / "index.html").read_text(encoding="utf-8")
html = tpl.replace("/*__AEGIS_DATA__*/ null", json.dumps(out, separators=(",", ":")))
dest = AEGIS / "app" / "aegis_healthy_standalone.html"
dest.write_text(html, encoding="utf-8")
kb = dest.stat().st_size / 1024
print(f"Wrote {dest.name} ({kb:.0f} KB) — {meta['n_groups']} groups, "
      f"{meta['n_genes']} genes, {meta['total_cells']:,} cells, tumors={meta['tumors']}")
