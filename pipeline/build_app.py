#!/usr/bin/env python3
"""
Aegis — App Builder
Injects data/aegis_data.json into app/template.html to produce a single
self-contained app/aegis.html (double-click to open; no server needed).
"""
import json
from pathlib import Path

AEGIS = Path(__file__).resolve().parent.parent
data = (AEGIS / "data" / "aegis_data.json").read_text(encoding="utf-8")
tpl = (AEGIS / "app" / "index.html").read_text(encoding="utf-8")

out = tpl.replace("/*__AEGIS_DATA__*/ null", json.dumps(json.loads(data), separators=(",", ":")))
dest = AEGIS / "app" / "aegis_standalone.html"
dest.write_text(out, encoding="utf-8")
kb = dest.stat().st_size / 1024
print(f"Wrote {dest.relative_to(AEGIS)} ({kb:.0f} KB, data embedded)")
