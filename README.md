# morphoDE — Bulge vs. Dent: a cross-species test of interface geometry → gene expression

*morphoDE = **morpho**logy × **d**ifferential **e**xpression: does the shape of a cell-type interface (its bulges and dents) correspond to a difference in gene expression?*

**Do the convex ("bulge") and concave ("dent") parts of a cell-type interface
express different genes?** We tested this deliberately naive hypothesis across
**11 spatial-transcriptomics datasets and 38 tissue interfaces**, from a human
tumor to a plant root tip, and found the correspondence holds broadly — it is a
general property of tissue architecture, not a tumor-specific phenomenon.

> Built with **Claude Science** for the 2026 *Built with Claude: Life Sciences*
> hackathon (Researcher Track). The biological question and the public datasets
> predate the event; **the analysis in this repository was carried out during the
> hackathon.**

![Cross-species geometry-breaking null survival](figures/crossspecies_tree_of_life.png)

*Each bar is one interface: the number of top differential-expression genes (out
of 12) whose bulge-vs-dent signal survives a geometry-breaking null. Red = human
tumor, blue = animal (embryo/brain/larva/regeneration), green = plant. Four
phylogenetically unrelated interfaces (human tumor, axolotl brain VLMC,
zebrafish notochord, planarian blastema `l4`) reach a perfect 12/12.*

---

## The headline result

| Metric | Value |
|---|---|
| Datasets tested | **11** |
| Interfaces tested | **38** |
| Datasets with ≥1 interface carrying signal | **11 / 11** |
| Interfaces with ≥1 surviving gene | **35 / 38** |
| Interfaces with ≥6/12 surviving | **19** |
| Interfaces at a perfect 12/12 | **4** |
| Interfaces with zero survivors | **3** |

Species/organs spanned: **human** (HNSCC metastatic lymph node, Open-ST),
**macaque** (cerebellum), **axolotl** (regenerating brain), **zebrafish**
(embryo), **mouse** (three independent embryo atlases), **fly** (larva + embryo,
Stereo-seq), **planarian** (regeneration), and **Arabidopsis** (root tip).

The 3 zero-survivor interfaces are not failures — they are the null *working*.
They have real expression differences (e.g. mouse Rostral neurectoderm,
max|lfc|=0.47) but those differences track **another spatial axis** (body-axis
position), not the bulge/dent geometry, so the geometry-breaking null correctly
rejects them.

---

## Why this is not trivial (the design that makes the result trustworthy)

A naive bulge-vs-dent comparison is dominated by three confounds. Each is
removed by construction:

1. **Cell-type composition.** Bulges and dents contain different cell types.
   → We measure differential expression **within a single cell type**, and
   separately decompose the total difference into a composition term and a
   within-cell-type term (Kitagawa/Oaxaca-style). In the human tumor, the split
   is **61 % composition / 39 % within-cell-type**; the within-cell-type part is
   the interesting residue.

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
   convex-hull score partly tracks stage/body size. → We always analyze **within
   the largest single stage/specimen stratum**.

The null is validated two ways on the human tumor (`figures/null_calibration_power.png`):
**false-positive rate is nominal** (p<0.05 in 3 % when expression is shuffled to
zero true signal) and **power is 100 % at log-fold 0.02** — far below the observed
effect sizes.

We never use curvature. Bulge/dent are defined purely by **convex-hull deviation**
(global protrusion) and, as a sensitivity check, **smoothed-surface displacement**
(local bumpiness).

**Causal direction is not identifiable** from a single-time-point snapshot and is
never claimed. We show only that a bulge/dent ↔ expression correspondence
survives after removing composition, spatial structure, location, and stage.

---

## What Claude Science did here

- **Reframed a failing analysis.** The continuous-curvature → expression analysis
  collapsed under SAR. Claude diagnosed *why* (SAR is unfair to a smooth
  geometric effect) and replaced it with the geometry-breaking null that puts the
  hypothesis on trial correctly.
- **Fanned the analysis across 11 datasets** as independent parallel sub-agents
  (fresh kernels), each handling one dataset's interface extraction, stratum
  detection, within-cell-type DE, and null test.
- **Packaged the method as a reusable skill** (`skill/interface-bulge-dent-de`)
  with a `kernel.py` of tested helpers, so the whole pipeline runs on a new
  dataset in one call.

---

## Reproducing

The pipeline is provided as an Agent Skill plus a batch driver.

```
skill/interface-bulge-dent-de/   # SKILL.md + kernel.py — the method
  SKILL.md                       # workflow, design rationale, when-to-use
  kernel.py                      # compute_bulge_dent_scores, tercile,
                                 # de_meandiff, decompose_de,
                                 # azimuthal_shift_null, null_calibration, power_curve
src/bulge_dent_driver.py         # per-dataset: interface → hull deviation →
                                 # target-domain cells → stratum → within-type DE + null
```

**Inputs required per dataset:** a reconstructed 3D interface mesh (vertices +
faces), per-cell signed distance to the interface, per-cell 3D coordinates,
per-cell **cell-type annotations**, and a log-normalized expression matrix. Mesh
construction and interface extraction (`st3d_interface.compute_interface`) come
from the upstream 3DCurv project; this repo consumes its output.

See [`DATA.md`](DATA.md) for the public sources of every dataset. Raw `.h5ad`
files are **not** redistributed here (size + upstream licensing); `DATA.md` gives
the accession/URL for each.

Minimal usage sketch:

```python
# after loading the skill (defines the helpers in your kernel)
from bulge_dent_driver import run_interface
de, null, meta = run_interface(adata, interface_spec, skillns=globals())
# de:   within-cell-type bulge-vs-dent log-fold per gene
# null: geometry-breaking null p-value per top gene
# meta: survivors, stratum, effect size, cell counts
```

---

## Results in this repo

| Path | Contents |
|---|---|
| `results/bd_all_meta.csv` | Per-interface summary (survivors, stratum, effect size, ρ between the two definitions, cell counts) — 38 rows |
| `results/bd_all_null.csv` | Per-interface × per-gene geometry-breaking null (456 rows) |
| `results/bd_all_de.csv` | All-interface within-cell-type DE (161,407 rows) |
| `results/per_dataset/bd_{de,null,meta}_<dataset>.csv` | Same, split by dataset (11 × 3) |
| `results/openst_detail/` | Deep evidence on the human tumor: within-tumor DE, composition decomposition, geometry-null table, power curve, two-definition sensitivity, GO ORA |
| `figures/` | The five figures (cross-species summary + the four-part openst proof + supporting panels) |
| `report/` | Full write-up (Japanese original + English translation) |

---

## Figures

- `crossspecies_tree_of_life.png` — the cross-species headline (above): every
  interface plotted against a cladogram, so the breadth across the tree of life
  is read directly.
- `tumor_interface_proof.png` — six-panel proof on the human tumor: (A) the
  reconstructed 3D interface mesh shown from two angles, colored bulge→dent;
  (B) composition decomposition; (C) within-tumor DE (bulge = hypoxia/invasion
  SLC2A1·EGLN3·PTHLH·KRT17; dent = immunoglobulin); (D) geometry-breaking null;
  (E) power curve; (F) two-definition sensitivity.
- `convexity_definitions.png` — the two curvature-free bulge/dent definitions and
  their agreement (ρ≈0.36 — complementary scales).
- `null_calibration_power.png` — null calibration (nominal false-positive rate)
  + spike-in power curve.
- `de_decomposition.png` — composition vs within-cell-type split of the total
  difference.

All figure labels are in English.

---

## License

MIT — see [`LICENSE`](LICENSE). Applies to the code, skill, and derived result
tables/figures in this repository. Underlying raw datasets retain their own
upstream licenses (see `DATA.md`).
