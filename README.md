# 🛡️ Aegis — Single-Cell Target Safety Atlas

**Query a cell-surface protein and instantly see where it is expressed across
healthy-tissue cell types — and, critically, whether those cells are dividing.**

Aegis is a lightweight, fully static web app for **off-tumor / on-target safety
assessment** of antibody–drug-conjugate (ADC) and other surface-targeted
therapeutics. For any surface gene it plots every healthy-tissue cell type as a
bubble:

| Encoding | Meaning |
|---|---|
| **X** | mean expression (log-normalised) |
| **Y** | **proliferation** (Tirosh S + G2M signature) |
| **Size** | % of cells expressing (raw count > 0) |
| **Colour** | cell type (or tissue / safety tier) |
| **Red ring** | Tier 1–2 (dose-limiting) tissue |

### Why the proliferation axis?
Most ADC payloads (topoisomerase-I inhibitors, auristatins, maytansinoids) kill
**dividing** cells. A target that sits on a *proliferating* normal cell type of a
*dose-limiting* organ — e.g. hematopoietic stem cells in bone marrow — is a
far bigger safety liability than the same target on quiescent cells. Aegis puts
that trade-off on a single chart: **upper-right + red ring = danger.**

---

## Quick start

The app is a static site — nothing to install to *view* it:

```bash
# from the aegis/ directory
python -m http.server 8787
# open http://127.0.0.1:8787/app/aegis.html
```

Or just open the self-contained `app/aegis.html` in a browser (data is embedded).

## Hosting

Aegis is **100% client-side** — no backend, no database, no runtime compute.
Serving it is just handing over static files:

- **GitHub Pages** (free): point Pages at `/app`.
- **Homelab**: `nginx`/`Caddy`, `python -m http.server`, or a
  `nginx:alpine` Docker image + the files (~10 MB image).
- Total footprint: app HTML (~20 KB) + Plotly (CDN, or ~3.5 MB self-hosted) +
  the data file (see below). A Raspberry Pi serves it comfortably.

The multi-GB single-cell `.h5ad` files are **only needed at build time** — they
never touch the server.

---

## How it's built

```
raw scBaseCount / CELLxGene h5ads  ──►  pipeline/build_summary.py  ──►  data/aegis_data.json
                                                                          │
                                                    pipeline/build_app.py ▼
                                                                    app/aegis.html
```

1. **`pipeline/build_summary.py`** — a *lightweight* pass (no PCA/UMAP) over a
   manifest of h5ads. Normalises, scores proliferation, and aggregates, per
   `(tissue, cell_type)`, the mean expression and % positive of every surface
   protein. Emits a compact `data/aegis_data.json` (a few MB).
2. **`pipeline/build_app.py`** — embeds that JSON into the app template to
   produce a single-file `app/aegis.html`.

```bash
python pipeline/build_summary.py --manifest manifests/local.csv --min-cells 20
python pipeline/build_app.py
```

Manifest format (`path,tissue`):

```csv
path,tissue
results/normal/liver/SRX9627330_annotated.h5ad,liver
results/normal/bone_marrow/SRX28908287_annotated.h5ad,bone_marrow
```

---

## Data

- **Source**: [scBaseCount](https://arcinstitute.org/) (Arc Institute Virtual
  Cell Atlas) / [CELLxGene](https://cellxgene.cziscience.com/) healthy-tissue
  single-cell RNA-seq. Cell-type labels are the native ontology annotations.
- **Surface proteins**: curated cell-surface list (CSPA / UniProt derived).
- **Proliferation**: Tirosh et al. 2016 S + G2M gene signatures via
  `scanpy.tl.score_genes_cell_cycle`.

## Tissue safety tiers

| Tier | Label | Example tissues |
|---|---|---|
| 1 | CRITICAL | bone marrow, heart, liver |
| 2 | SERIOUS | kidney, lung, colon, small intestine, pancreas |
| 3 | MANAGEABLE | skin, spleen, thymus |
| 4 | LOW | muscle, adipose, prostate |

---

## Disclaimer
Aegis is a **research tool** built on public data for hypothesis generation.
It is **not** a substitute for experimental validation or IND-enabling
toxicology. Expression ≠ protein-level accessibility, and single-cell dropout
can under-report true positivity.

## License
MIT — see [LICENSE](LICENSE).
