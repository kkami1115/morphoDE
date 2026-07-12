# morphoDE — Bulge vs. Dent: a cross-species test of interface geometry → gene expression

*morphoDE = **morpho**logy × **d**ifferential **e**xpression: does the shape of a cell-type interface (its bulges and dents) correspond to a difference in gene expression?*

**Do the convex ("bulge") and concave ("dent") parts of a cell-type interface
express different genes?** We tested this deliberately naive hypothesis across
**11 spatial-transcriptomics datasets and 126 interface–strata**, from a human
tumor to a plant leaf, and found the correspondence holds almost everywhere — it
is a general property of tissue architecture, not a tumor-specific phenomenon,
and it is reproducible across developmental stages, regeneration time points, and
individual specimens of the same organism.

> Built with **Claude Science** for the 2026 *Built with Claude: Life Sciences*
> hackathon (Researcher Track). The biological question and the public datasets
> predate the event; **the analysis in this repository was carried out during the
> hackathon.**

![Cross-species geometry-breaking null survival](figures/crossspecies_tree_of_life.png)

*Each dot is one interface–stratum: the number of top differential-expression
genes (out of 12) whose bulge-vs-dent signal survives a geometry-breaking null.
Rows are ordered phylogenetically (vertebrates → invertebrates → plant). Red =
human tumor, blue = animal, green = plant; short vertical line = per-species
median; hollow = zero survivors. Signal is present across the whole tree of life.*

---

## The headline result

| Metric | Value |
|---|---|
| Datasets tested | **11** |
| Interface–strata tested | **126** |
| Datasets with ≥1 interface carrying signal | **11 / 11** |
| Interface–strata with ≥1 surviving gene | **125 / 126** |
| Interface–strata with ≥6/12 surviving | **97** |
| Interface–strata at a perfect 12/12 | **16** |
| Interface–strata with zero survivors | **1** |

Species/organs spanned: **human** (HNSCC metastatic lymph node, Open-ST),
**macaque** (cerebellum), **axolotl** (regenerating brain), **zebrafish**
(embryo), **mouse** (three independent embryo atlases), **fly** (larva + pupa,
Stereo-seq), **planarian** (regeneration), and **Arabidopsis** (leaf).

An **interface–stratum** is one cell-type interface measured within one
developmental stage / regeneration time point / specimen. Multi-stage or
multi-specimen datasets contribute one row per stratum, so the same interface is
tested repeatedly across biological conditions — signal that holds across all of
them is far stronger evidence than a single snapshot. The single zero-survivor
(mouse digital-embryo *Paraxial mesoderm*, one specimen) is the null *working*:
its expression differences track another spatial axis, not the bulge/dent
geometry, so the geometry-breaking null correctly rejects it.

---

## Why this is not trivial (the design that makes the result trustworthy)

A naive bulge-vs-dent comparison is dominated by three confounds. Each is
removed by construction:

1. **Cell-type composition.** Bulges and dents contain different cell types.
   → We measure differential expression **within a single cell type**, and
   separately decompose the total difference into a composition term and a
   within-cell-type term (Kitagawa/Oaxaca-style). On the human tumor's genuine
   interface the split is **93 % within-cell-type / 7 % composition** — the
   bulge/dent difference is overwhelmingly a change in what a single tumor cell
   type expresses, not a reshuffling of cell types.

2. **Spatial autocorrelation & "it's just a different location".** Neighboring
   cells are similar, so any two regions differ. The standard fix — a
   spatial-autoregressive (SAR/GM_Lag) model — **is the wrong test here**: if
   geometry drives expression, that effect is itself spatially smooth and nearly
   collinear with the SAR lag term, so SAR absorbs it and reports
   non-significance. We instead use a **geometry-breaking null**: per-section
   azimuthal circular shift of the bulge/dent score. This preserves spatial
   smoothness but destroys the score↔position correspondence, asking exactly the
   right question — *"is the geometry-label ↔ expression correspondence stronger
   than a spatially-smooth random relabeling?"*

3. **Developmental stage / specimen (atlases only).** In a multi-stage atlas the
   interface geometry partly tracks stage/body size. → We analyze **each stage /
   specimen as its own interface–stratum** and report all of them, so a claim only
   stands if it holds across biological conditions, not because we cherry-picked a
   stratum.

The null is validated two ways on the human tumor (panel E of the proof figure):
**false-positive rate is nominal** (p<0.05 in ~5 % when expression carries no true
geometric signal) and **power is 100 % at log-fold 0.02** — far below the observed
effect sizes.

**The interface is an open surface, not a closed solid.** An earlier version used
the convex hull of the domain, which closes the tissue block into a solid and
mislabels the physical cut-faces of the acquired slab as "bulges." We now extract
the **genuine open interface sheet** — keeping only surface faces that have
sampled tissue on *both* sides (the real domain/non-domain boundary), which
discards the block's cut-faces. The sheet's voxel resolution is **density-adaptive**
(sized to the local cell density, no per-dataset hand-tuning). Bulge/dent are the
**normal displacement** of this open sheet; a second **mean-curvature** definition
agrees gene-for-gene (ρ≈0.86, 12/12 same sign). The sign is fixed to the physical
surface orientation: **positive = convex = the domain protrudes into its neighbor
(bulge)**, negative = concave = the domain recedes (dent).

**Causal direction is not identifiable** from a single-time-point snapshot and is
never claimed. We show only that a bulge/dent ↔ expression correspondence
survives after removing composition, spatial structure, location, and stage.

![Six-panel proof on the human tumor interface](figures/tumor_interface_proof.png)

*The full evidence on the human tumor (Open-ST). **A** the genuine open interface
sheet from two angles (red = bulge/convex, blue = dent/concave). **B** the
bulge–dent difference is 93 % within cell type, only 7 % composition. **C**
within-tumor DE: bulge/convex = IGKC·IGHG3·DHCR7 (immunoglobulin / cholesterol),
dent/concave = PTHLH·KRT17·FTH1 (hypoxia / invasion / keratinization). **D** all
12 genes clear the geometry-breaking null (colored = observed, grey = null). **E**
false-positive rate is nominal and power reaches 100 % by log-fold 0.02. **F** the
primary (normal-displacement) definition and a mean-curvature cross-check agree
gene-for-gene (ρ≈0.86).*

The interface itself — one open curved sheet, Z-truncated by the acquired tissue
block, not a closed solid:

![The tumor interface as one open curved sheet](figures/concept_3d_still.png)

---

## What Claude Science did here

- **Reframed a failing analysis.** The continuous-curvature → expression analysis
  collapsed under SAR. Claude diagnosed *why* (SAR is unfair to a smooth
  geometric effect) and replaced it with the geometry-breaking null that puts the
  hypothesis on trial correctly.
- **Caught and fixed a structural flaw mid-project.** The reviewer (the user)
  noticed the human-tumor interface was being rendered as a closed blob when the
  tissue was never acquired as a closed solid. Claude verified the problem
  quantitatively (bulge cells' 60 µm neighborhoods were only ~4 % non-tumor — they
  faced *unsampled space*, not the real boundary), rebuilt the method on the
  genuine open interface sheet, and re-ran all 126 interface–strata. The
  headline finding survived and sharpened (within-cell-type share 39 % → 93 %).
- **Verified the geometry sign against data, not intuition.** When a red/blue
  color looked backwards, Claude did not guess — it probed the sheet normal
  (tumor-fraction 0.62 vs 0.33 across the surface) to fix the bulge/dent sign from
  the physical surface orientation alone, then reported whatever biology followed.
- **Fanned the analysis across 11 datasets** as independent parallel sub-agents
  (fresh kernels), each handling one dataset's interface extraction, stratum
  detection, within-cell-type DE, and null test.
- **Packaged the method as a reusable skill** (`skill/interface-bulge-dent-de`)
  with a `kernel.py` of tested helpers, so the whole pipeline runs on a new
  dataset in one call.

---

## Reproducing

The whole method is three functions in `src/open_interface.py`, plus a skill of
scoring/decomposition helpers.

```
src/open_interface.py            # the method, self-contained:
  extract_open_interface(...)    #   labelled point cloud → genuine open sheet
                                 #   (both-sided faces only, density-adaptive pitch)
  band_scores(...)              #   normal displacement (convex=bulge) on the
                                 #   nearest-1/3 distance-quantile interface band
  de_and_null(...)              #   within-cell-type DE + geometry-breaking null
skill/interface-bulge-dent-de/   # SKILL.md + kernel.py — scoring/decomposition
  kernel.py                      #   tercile, de_meandiff, decompose_de,
                                 #   null_calibration, power_curve
src/legacy_closed_mesh/          # the superseded closed-hull driver (provenance)
```

**Inputs required per dataset:** per-cell 3D coordinates, a per-cell domain /
non-domain label (the interface is extracted *from the points*, no pre-built mesh
needed), per-cell **cell-type annotations**, per-cell section id, and a
log-normalized expression matrix.

See [`DATA.md`](DATA.md) for the public sources of every dataset. Raw `.h5ad`
files are **not** redistributed here (size + upstream licensing); `DATA.md` gives
the accession/URL for each.

Minimal usage sketch:

```python
from open_interface import extract_open_interface, band_scores, de_and_null
# XYZ: (n,3) coords; is_domain: (n,) bool; ann/section/Xn per cell
V, F, genuine_frac = extract_open_interface(XYZ, is_domain)         # open sheet
band, score, in_band = band_scores(XYZ, is_domain, V, F)            # bulge/dent
de, null, meta = de_and_null(Xn[band], genes, score, XYZ[band], section[band])
# de:   within-cell-type bulge-vs-dent log-fold per gene
# null: geometry-breaking null p-value per top gene
# meta: survivors, effect size, band cell counts
```

---

## Results in this repo

| Path | Contents |
|---|---|
| `results/bd_final_headline.csv` | **Primary result:** per-interface–stratum survivors / tested / band cells / genuine-fraction — 126 rows (open-surface method) |
| `results/bd_final_full.json` | Per-interface–stratum gene list, log-fold, and null p-value for the 12 tested genes each |
| `results/legacy_closed_mesh/` | The earlier closed-convex-hull results (38 interfaces) kept for provenance — **superseded**; see the method note in `DATA.md` |
| `report/` | Full write-up (Japanese original + English translation) |

`bd_final_headline.csv` / `bd_final_full.json` are the results the paper and
figures are built on. The `legacy_closed_mesh/` tables reproduce the earlier
closed-surface analysis and are retained only so the correction is auditable.

---

## Figures

- `crossspecies_tree_of_life.png` — the cross-species headline (above): every
  interface–stratum plotted by phylogeny, so the breadth across the tree of life
  is read directly.
- `tumor_interface_proof.png` — six-panel proof on the human tumor: (A) the
  genuine **open** interface sheet from two angles, colored bulge→dent;
  (B) composition decomposition (93 % within cell type); (C) within-tumor DE
  (bulge/convex = IGKC·IGHG3·DHCR7 immunoglobulin/cholesterol; dent/concave =
  PTHLH·KRT17·FTH1 hypoxia/invasion/keratinization); (D) geometry-breaking null,
  12/12; (E) power curve; (F) the two shape definitions — normal displacement and
  mean curvature — agree gene-for-gene (ρ≈0.86).
- `concept_3d_still.png` / `concept_3d_rotate.mp4` — the interface as one open
  curved sheet, Z-truncated by the acquired block (demo hook).

All figure labels are in English. The proof figure is self-contained; earlier
standalone panels (`convexity_definitions.png`, `null_calibration_power.png`,
`de_decomposition.png`) were folded into it.

---

## License

MIT — see [`LICENSE`](LICENSE). Applies to the code, skill, and derived result
tables/figures in this repository. Underlying raw datasets retain their own
upstream licenses (see `DATA.md`).
