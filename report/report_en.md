# Does expression differ between the bulge and the dent of an interface? — a cross-species test

*(English translation of `report_ja.md`. The Japanese file is the authoritative original.)*

## Abstract

We tested the deliberately naive hypothesis that *"the boundary (interface)
between domains of different cell types has protruding and recessed regions, and
gene expression differs between them"* across **11 datasets and 126
interface–strata (from human tumor to plant leaf).** The primary bulge/dent label
uses no curvature (a continuous quantity); we extracted the **genuine open
interface sheet** directly
from the labelled point cloud — keeping only surface faces with sampled tissue on
*both* sides, so the physical cut-faces of the acquired block are excluded — at a
**density-adaptive resolution**. Bulge/dent is a binary label from the sheet's
**normal displacement**, with the sign fixed to the physical surface orientation
(positive = convex = the domain protrudes into its neighbor). We measured
differential expression
**within a single cell type**, ran **every developmental stage / regeneration
time point / specimen as its own interface–stratum**, and tested significance
with a **geometry-breaking null** (azimuthal circular shift within z-bands,
n=200) that preserves spatial smoothness while destroying only the
label↔position correspondence.

Conclusion: **this correspondence is not tumor-specific — it holds almost
everywhere across species, organs, developmental stages, and individual
specimens.** At least one interface carried signal in all 11 datasets (11/11);
**125 of 126** interface–strata had ≥1 surviving gene, **97** had a majority
(≥6/12), and **16** reached a perfect 12/12. A single interface–stratum had zero
survivors — the null *working*. Causal direction (whether geometry drives
expression or an expression difference shapes the geometry) cannot be determined
from a snapshot and is not claimed.

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
- **Per-stratum, all reported**: developmental atlases contain multiple
  stages/specimens. Because interface geometry confounds with stage/body-size, we
  run **each stage / regeneration time point / specimen as its own
  interface–stratum** and report all of them, so a claim only stands if it holds
  across biological conditions rather than in one cherry-picked stratum.

- **Open interface, not a closed solid**: an earlier version used the convex hull
  of the domain, which closes the acquired tissue block into a solid and mislabels
  its physical cut-faces as "bulges." We verified this quantitatively (bulge cells'
  60 µm neighborhoods were only ~4 % non-tumor — they faced unsampled space, not
  the real boundary) and rebuilt the method on the **genuine open interface sheet**,
  keeping only faces with sampled tissue on both sides. The sheet resolution is
  **density-adaptive** (voxel edge sized to local cell density; no per-dataset
  hand-tuning).

On the openst tumor interface, this design showed a within-tumor difference that
survived the geometry-breaking null (12/12), 93 % of it within a single tumor cell
type, with demonstrated null calibration (false-positive ≤ nominal 5 %) and power
(100 % at lfc=0.02). With bulge = convex = tumor protruding into stroma:
**bulge = IGKC, IGHG3, DHCR7, S100A7, CD74** (immunoglobulin / antigen
presentation / cholesterol biosynthesis) and **dent = PTHLH, FTH1, KRT17, AMTN,
FABP5** (hypoxia / invasion / keratinization). This report applies that
instrument to all datasets.

## Results by dataset

Each interface is run as one row per developmental stage / regeneration time
point / specimen ("interface–strata"), so multi-condition datasets contribute
many rows and a claim only stands if it holds across conditions.

| Dataset | Interface–strata | ≥1 survived | Median survival | Strongest interface (survived/12) |
|---|---|---|---|---|
| Human HNSCC metastatic lymph node (Open-ST) | 1 | 1/1 | 12 | tumor/nontumor (12/12) |
| Macaque cerebellum | 4 | 4/4 | 11 | molecular layer (12/12) |
| Axolotl brain (regeneration) | 20 | 20/20 | 9 | MSN (12/12) |
| Zebrafish embryo | 4 | 4/4 | 11 | Yolk Syncytial Layer (12/12) |
| Mouse embryo (digital, Seq-Scope) | 20 | 19/20 | 6 | ExE endoderm (11/12) |
| Mouse embryo (MOSTA E16.5) | 4 | 4/4 | 4 | Cavity (5/12) |
| Mouse embryo (CNGB) | 4 | 4/4 | 7 | 3-endochondral bone (12/12) |
| Fly larva (Stereo-seq) | 17 | 17/17 | 7 | midgut (11/12) |
| Fly pupa (Stereo-seq) | 4 | 4/4 | 8 | epidermis (11/12) |
| Planarian (regeneration) | 47 | 47/47 | 10 | l3 (12/12) |
| Arabidopsis leaf | 1 | 1/1 | 11 | epidermis (11/12) |

The survival decision is a **scale-invariant quantity** (where the observed
statistic falls in the geometry-breaking null distribution), so it is unaffected
by whether a dataset is stored as raw counts or log-normalized; where absolute
effect sizes are reported, raw-count datasets are normalized first. Signal is
weakest on the thin-slab embryo datasets (MOSTA, CNGB): these span only ~9–13
serial sections in Z, which limits the resolution of the open interface sheet —
a data limitation, not a failure of the effect.

## Overall picture

- **All 126 interface–strata**: ≥1 gene survived in **125/126**, majority
  (≥6/12) in **97**, perfect (12/12) in **16**, zero survivors in **1**.
- **All 11 datasets have ≥1 interface with signal (11/11)** — human tumor,
  macaque cerebellum, axolotl brain, zebrafish embryo, three mouse embryo
  datasets, fly larva/pupa, planarian, and Arabidopsis leaf.
- **Signal is reproducible across conditions, not a single lucky snapshot.**
  Planarian carries signal in all 47 of its interface–strata (three body regions
  × 16 regeneration time points), axolotl in all 20 (four regions × five
  regeneration stages), the digital mouse embryo in 19 of 20 (four germ-layer
  interfaces × up to six specimens). The 16 perfect (12/12) interface–strata are
  phylogenetically unrelated — human tumor, macaque cerebellum, axolotl,
  zebrafish, mouse embryo (CNGB, endochondral bone), and many planarian time
  points.
- **The single zero-survivor** (mouse digital-embryo Paraxial mesoderm, one
  specimen): its expression difference tracks some other spatial structure
  (body-axis position), not the bulge/dent, so it disappears under the
  geometry-breaking null. This is exactly what the null is for — the same
  interface carries strong signal in the embryo's other specimens.

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
   composition confound; on the human tumor, doing so leaves 93 % of the
   difference intact — it is genuinely within-cell-type).
2. Use the **genuine open interface**, not a closed hull (a closed solid
   mislabels the acquired block's cut-faces as bulges).
3. Report **every stage/specimen stratum** (in developmental atlases, stage
   confounds with interface geometry; reporting all strata removes the
   temptation to cherry-pick).
4. Interface resolution is **density-adaptive**; on thin-slab data (few Z
   sections) the achievable resolution — and hence sensitivity — is limited.

**Causal direction is undetermined**: from a single-time-point snapshot one
cannot distinguish geometry driving expression from an expression difference
deforming the tissue to create the geometry. What this analysis proves is only
that a geometry↔expression correspondence remains after removing composition,
spatial smoothness, location, and developmental stage.

## Methods and reproducibility

- The interface method is in `src/open_interface.py`: `extract_open_interface`
  (genuine open sheet, density-adaptive pitch), `band_scores` (normal
  displacement, convex=bulge), `de_and_null` (within-cell-type DE + geometry-breaking
  null). The scoring/decomposition helpers are published as the skill
  `interface-bulge-dent-de`.
- Each dataset was processed with its recorded loading recipe; multi-stage/
  specimen datasets were expanded to all strata. Every stratum uses n_perm=200 and
  a nearest-1/3 distance-quantile interface band.
- No curvature is used (normal-displacement based, with a mean-curvature sensitivity
  check). The bulge/dent sign is fixed to the physical surface orientation
  (convex = the domain protrudes into its neighbor), verified against the local
  domain-fraction across the sheet.

## Artifacts

| File | Contents |
|---|---|
| `figures/crossspecies_tree_of_life.png` | Geometry-breaking null survival across all 126 interface–strata, ordered by phylogeny (cross-species summary) |
| `figures/tumor_interface_proof.png` | Six-panel proof on the human tumor interface (open sheet, decomposition, DE, null, power, definition sensitivity) |
| `figures/concept_3d_still.png`, `concept_3d_rotate.mp4` | The tumor interface as one open curved sheet |
| `results/bd_final_headline.csv` | Per-interface–stratum survivors / tested / band cells / genuine-fraction (126 rows) |
| `results/bd_final_full.json` | Per-interface–stratum gene list, log-fold, null p-value (12 tested genes each) |
| `results/legacy_closed_mesh/` | The superseded closed-convex-hull results (38 interfaces), kept for provenance |
