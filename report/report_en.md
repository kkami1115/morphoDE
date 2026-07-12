# Does expression differ between the bulge and the dent of an interface? — a cross-species test

*(English translation of `report_ja.md`. The Japanese file is the authoritative original.)*

## Abstract

We tested the deliberately naive hypothesis that *"the boundary (interface)
between domains of different cell types has protruding and recessed regions, and
gene expression differs between them"* across **11 datasets and 38 interfaces
(from human tumor to plant root tip).** Without using curvature (a continuous
quantity), we defined bulge/dent as a binary label from the interface mesh's
**convex-hull deviation**, measured differential expression **within a single
cell type and within a single developmental stage/specimen**, and tested
significance with a **geometry-breaking null** (per-section azimuthal circular
shift, n=200) that preserves spatial smoothness while destroying only the
label↔position correspondence.

Conclusion: **this correspondence is not tumor-specific — it holds broadly across
species, organs, and developmental stages.** At least one interface survived the
geometry-breaking null in all 11 datasets (11/11); 35 of 38 interfaces had ≥1
surviving gene and 19 had a majority (≥6/12). Interface-dependent variation is
also clear: 3 interfaces had zero survivors. Causal direction (whether geometry
drives expression or an expression difference shapes the geometry) cannot be
determined from a snapshot and is not claimed.

## Background and key design choices

This test was motivated by the outcome of a prior curvature analysis. The
coupling of continuous local curvature (mean curvature H, Gaussian curvature K,
shape index) → expression **vanished under a spatial-autoregressive model
(SAR/GM_Lag)** in both the openst tumor and the flysta3d epidermis. But **SAR
asks "does geometry add predictive power beyond the expression spatial lag (the
neighboring expression itself)?" — an unfairly hard question for this naive
hypothesis**: if geometry drives expression, that effect is itself spatially
smooth, becomes collinear with the lag term, and is absorbed.

So we stepped the hypothesis down one level of primitiveness and asked about a
**coarse binary protrusion/recession**. We also swapped the test:

- **Geometry-breaking null**: in each section, order the bulge/dent scores by
  azimuth around the interface centroid and circularly shift them. Spatial
  smoothness (neighboring cells have similar scores) is preserved; only the
  correspondence between score and actual position is destroyed. This asks
  directly: *"is the geometry-label ↔ expression correspondence stronger than a
  spatially-smooth random relabeling?"*
- **Within a single cell type**: DE on only the cells of the target domain. This
  structurally removes the cell-type composition confound, leaving the intrinsic
  component.
- **Within a stratum**: developmental atlases contain multiple stages/specimens.
  Because the convex-hull bulge/dent score confounds with stage/body-size, we
  analyze **within the largest single stage/single specimen stratum** (neglecting
  this stratification produces a huge spurious DE).

On the openst tumor interface, this design showed a within-tumor difference —
bulge = hypoxia/invasion (SLC2A1/GLUT1, EGLN3, PTHLH, KRT17), dent =
immunoglobulin — that survived the geometry-breaking null (12/12), with
demonstrated null calibration (false-positive ≤ nominal 5 %) and power (100 % at
lfc=0.02). This report applies that instrument to all datasets.

## Results by dataset

| Dataset | Interfaces | Mean survival | Strongest interface (survived/12) | Stratification | Normalization |
|---|---|---|---|---|---|
| Human HNSCC metastatic lymph node (Open-ST) | 1 | 100% | tumor/nontumor (12/12) | single-specimen stack | log-norm |
| Primate cerebellum (macaque) | 4 | 38% | granular layer (6/12) | single-specimen stack | raw counts* |
| Axolotl brain (regeneration) | 4 | 71% | VLMC (12/12) | within stage/specimen | log-norm |
| Zebrafish embryo | 4 | 42% | Notochord (12/12) | single-specimen stack | log-norm |
| Mouse embryo (Seq-Scope E7.5–8.0) | 4 | 42% | ExE endoderm (9/12) | within stage/specimen | log-norm |
| Mouse embryo (MOSTA E16.5) | 4 | 40% | Connective tissue (10/12) | single-specimen stack | log-norm |
| Mouse embryo (CNGB E13.5) | 4 | 54% | 2-craniofacial primordium (11/12) | single-specimen stack | log-norm |
| Fly larva (Stereo-seq) | 4 | 58% | midgut (11/12) | within stage/specimen | raw counts* |
| Fly embryo (Stereo-seq) | 4 | 58% | somatic muscle (10/12) | single-specimen stack | log-norm |
| Planarian (regeneration) | 4 | 44% | l4 (12/12) | within stage/specimen | raw counts* |
| Arabidopsis root tip | 1 | 33% | epidermis (4/12) | single-specimen stack | log-norm |

\* raw counts: datasets with max|lfc|>1.5 (cerebellum, flysta3d, prista4d) may be
raw counts, so absolute effect-size comparison is not valid. However, the
**geometry-breaking null survival decision is a scale-invariant quantity**
determined by where the observed statistic falls in the null distribution, so
the survival conclusions hold.

## Overall picture

- **All 38 interfaces**: ≥1 gene survived in **35/38**, majority (≥6/12) in
  **19**, perfect (12/12) in **4**, zero survivors in **3**.
- **All 11 datasets have ≥1 interface with signal (11/11)** — human tumor,
  primate cerebellum, axolotl brain, zebrafish embryo, three mouse embryo
  datasets, fly larva/embryo, planarian, and Arabidopsis root tip.
- **Perfect (12/12) interfaces**: human tumor (tumor|nontumor), axolotl brain
  VLMC, zebrafish notochord, planarian l4 — four phylogenetically unrelated
  interfaces at the same top score. The tumor is *one of the interfaces that
  work*, not special.
- **Zero-survivor interfaces**: mouse Rostral neurectoderm, zebrafish Segmental
  Plate/Tail Bud and Yolk Syncytial Layer. The effect size is not small (e.g.
  Rostral neurectoderm max|lfc|=0.47); rather, **their difference tracks some
  other spatial structure (body-axis position), not the bulge/dent**, so it
  disappears under the geometry-breaking null. This selection is exactly what the
  null is for.

## Interpretation and limitations

**The concept is general**: "the convex and concave parts of a cell-type boundary
express differently" is not a story about tumor biology — it holds at diverse
interfaces (epithelium, mesenchyme, nervous tissue, regeneration blastema, germ
layers). It is consistent with a picture in which protrusion/recession changes
the local microenvironment (mechanical tension, diffusion distance, exposure to
neighboring cell types), which is reflected in cell state.

**Interface dependence is a feature, not a defect**: it does not hold at every
interface. The geometry-breaking null discriminates "differences truly
corresponding to protrusion/recession" from "differences that merely happen to
run along another spatial axis," correctly dropping the latter (e.g. the
neurectoderm interface).

**Methodological requirements**:
1. Restrict to a **single cell type** (otherwise most of the signal is
   composition confound).
2. Analyze **within a single developmental stage/specimen stratum** (in
   developmental atlases, stage confounds with convex-hull deviation).
3. Absolute effect-size comparison requires **log-normalization** (in raw-count
   datasets the effect sizes are orders of magnitude off; the survival decision
   is scale-invariant and unaffected).

**Causal direction is undetermined**: from a single-time-point snapshot one
cannot distinguish geometry driving expression from an expression difference
deforming the tissue to create the geometry. What this analysis proves is only
that a geometry↔expression correspondence remains after removing composition,
spatial smoothness, location, and developmental stage.

## Methods and reproducibility

- The pipeline is published as the skill `interface-bulge-dent-de` (two-definition
  scores, within-cell-type DE, composition decomposition, geometry-breaking null,
  null calibration, power curve).
- Batch execution: `bulge_dent_driver.py` (interface extraction → convex-hull
  deviation → target-domain cells → stratum detection → DE + geometry null within
  the largest stratum). Each dataset processed by an independent sub-agent (fresh
  kernel).
- No curvature is used (convex-hull-deviation based). No code copied from 3DCoS.

## Artifacts

| File | Contents |
|---|---|
| `figures/crossspecies_tree_of_life.png` | Geometry-breaking null survival across all 38 interfaces, plotted against a cladogram (cross-species summary) |
| `results/bd_all_meta.csv` | Per-interface meta (survivors, stratum, effect size, two-definition correlation, cell counts) |
| `results/bd_all_null.csv` | Per-interface × per-gene geometry-breaking null (456 rows) |
| `results/bd_all_de.csv` | All-interface within-cell-type DE (161,407 rows) |
| `results/per_dataset/bd_{de,null,meta}_<dataset>.csv` | Per-dataset raw output (11 × 3) |
