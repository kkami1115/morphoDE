---
name: interface-bulge-dent-de
description: Test whether gene expression differs between the bulging (protruding) and denting (recessed) parts of a 3D tissue interface — a naive, curvature-free hypothesis. Use when you have a reconstructed 3D interface mesh (e.g. tumor boundary) plus per-cell type annotations and want to prove convex vs concave regions express differently, controlling for cell-type composition, spatial autocorrelation, and mere location. Defines bulge/dent two ways (convex-hull deviation = global, smoothed-surface displacement = local), decomposes the total difference into composition vs within-cell-type parts, and tests significance with a geometry-destroying null (per-section azimuthal circular shift) instead of SAR — SAR absorbs a smooth geometric effect into its lag term and is unfair to this hypothesis. Includes null calibration + synthetic spike-in power curve and two-definition sensitivity. Causal direction (geometry drives expression vs vice versa) is not identifiable from a snapshot and is never claimed.
---

# Interface Bulge-vs-Dent Differential Expression

Prove that the **protruding (bulge)** and **recessed (dent)** parts of a 3D
tissue interface express differently, at three levels of rigor, without using
curvature. Built for serial-section spatial transcriptomics with a
reconstructed interface mesh.

## When to use
- You have a reconstructed 3D interface (mesh + per-cell signed distance / SDF)
  and per-cell **cell-type annotations**.
- You want the naive claim "bulge vs dent express differently" **proven**, not
  just observed — controlling for composition, spatial autocorrelation, and the
  fact that bulge/dent are simply different locations.

## Why NOT SAR (the key design choice)
A spatial-autoregressive (SAR / GM_Lag) model asks "does the geometry label add
predictive power *beyond the expression spatial lag*?" If geometry drives
expression, that effect is itself spatially smooth and nearly collinear with the
lag term, so SAR absorbs it and reports non-significance — an unfair test for
this hypothesis. Instead use a **geometry-destroying null**: keep the spatial
smoothness of the bulge/dent score but destroy its correspondence to actual
position (per-section azimuthal circular shift). This asks the right question:
"is the geometry-label ↔ expression correspondence stronger than a
spatially-smooth random relabeling?"

## Causal direction — always disclaimed
From a single-time-point snapshot you cannot tell whether geometry drives
expression or an expression difference shapes the geometry. Claim only that a
correspondence exists after removing composition / spatial structure / location.

## Workflow (helpers in kernel.py)

1. **Two bulge/dent definitions** — `compute_bulge_dent_scores(vertices, faces, vertex_normals)`
   returns `hull_dev` (signed distance to 3D convex hull; interior<0=dent,
   surface≈0=bulge; GLOBAL protrusion) and `smooth_disp` (outward displacement
   from a strongly Taubin-smoothed surface; +=bulge; LOCAL bumpiness). Report
   their correlation — they capture different scales.
   IMPORTANT: verify vertex normals are OUTWARD (see helper; it fixes them).

2. **Assign scores + covariates to cells** — for each interface cell, take the
   nearest-face mean score. Build covariates: SDF (depth), local density,
   vessel distance (nearest vascular-annotation cell = perfusion proxy),
   section. Tercile-label into bulge / mid / dent with `tercile(score)`.
   For a "thick tumor layer" analysis, restrict to the tumor-side SDF band
   (e.g. SDF in (-300, 0)) so tumor cells populate both bulge and dent.

3. **Region-level DE (composition-in)** — bulge-band vs dent-band mean-log-fold
   over ALL cell types; also quantify the cell-type composition difference. This
   is the TOTAL effect.

4. **Within-cell-type DE** — restrict to ONE cell type (e.g. Tumor); the DE now
   has composition confounding structurally removed = intrinsic component.

5. **Decompose total = composition + within-cell-type** with
   `decompose_de(X, labels, celltype, genes)` (symmetric Kitagawa/Oaxaca-style):
   composition term = Σ(Δ proportion)·mean-expr, within term = Σ(mean prop)·(within diff).

6. **Geometry-destroying null** — `azimuthal_shift_null(cells, X, score_col, n_perm)`
   builds the null distribution of the DE statistic by per-section azimuthal
   circular shift. Run on BOTH the composition-in and within-cell-type versions.

7. **Prove the test works** —
   `null_calibration(X, cells, score_col)`: shuffle expression across cells
   (true signal = 0) → p<0.05 must be ~nominal 5%.
   `power_curve(X, cells, score_col, effects)`: spike a known log-fold into the
   bulge side → detection rate vs effect size.

8. **Two-definition sensitivity** — rerun the null test with `label_disp`
   (smoothed-displacement) and compare to `label_hull` (convex-hull). Report
   which genes are robust across definitions and which are scale-specific.

9. **GO annotation** (optional, needs network) — g:Profiler ORA on the
   within-cell-type up/down gene lists via a delegated sub-agent (fresh kernel
   avoids the proxy-binding gotcha).

## Stratify by sample / developmental stage (REQUIRED for atlases)

The concept is NOT tumor-specific — it works on any interface between
cell-type domains (e.g. epidermis vs rest in a fly embryo). BUT: with a
single tumor specimen the "same stage, same specimen" condition holds
automatically; with a developmental atlas it does not. If the 3D stack mixes
developmental stages or specimens, the convex-hull bulge/dent score will
partly track stage/body-size, and you get a huge spurious DE that the
geometry null only partially removes. **Always run within a single
stage/specimen stratum.** Reference: flysta3d epidermis over all 90 sections
(5 stages mixed) gave max|lfc|=1.06 but only 4/10 survived the geometry null;
restricted to the E16-18h stage alone, max|lfc|=0.84 and 10/10 survived.

## Interpreting results
- A robust bulge-vs-dent signal survives the geometry-destroying null with
  observed |lfc| many-fold the null SD, is not killed by within-cell-type
  restriction, and is consistent across both definitions.
- Composition typically dominates (it did in openst: 61% composition / 39%
  within-cell-type); the within-cell-type component is the interesting part.
- Scale matters: global (convex-hull) vs local (smoothed) can disagree; hypoxia
  markers tracked the GLOBAL protrusion in openst (bulk far from vasculature),
  not local bumps.

## Reference result (openst HNSCC lymph-node, GSE251926)
Tumor-cell bulge = hypoxia/invasion (SLC2A1, EGLN3, PTHLH, KRT17, TGFBI);
dent = immunoglobulin. All target genes significant under the geometry null
(p=0.005); null calibration nominal; power 100% at lfc=0.02.
