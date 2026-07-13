#!/usr/bin/env python3
"""morphoDE — one-call demo (self-contained; for the shot-4 beat-3 recording).

Runs the whole method end to end on ONE interface and prints the result:
  labelled 3-D point cloud -> genuine open interface sheet -> bulge/dent scores
  -> within-cell-type DE -> geometry-breaking null.

Zebrafish Yolk Syncytial Layer (10 hpf) — reproduces the published 12/12,
small and fast (~2 s once the file is present).

Reproduces the exact published headline number for one interface, so it is the
ground-truth anchor for the fuller `scripts/reproduce.py` runner.

    conda activate morphode
    python demo/run_demo.py                       # looks for zf10 under data/
    python demo/run_demo.py /path/to/cache        # or point at a decompressed cache

Needs the zebrafish file `zf10_stereoseq.h5ad` (or `.h5ad.zst`); fetch it with
`python scripts/download.py zesta_zebrafish`.
"""
import os, sys, time, subprocess, glob
import anndata as ad, numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(REPO, "src"))
from open_interface import extract_open_interface, band_scores, de_and_null
from st3d_loader import assemble_3d, needs_registration, register_sections_rigid

DOMAIN = "Yolk Syncytial Layer"

# Locate the zf10 file: CLI arg > env cache > repo data/. Accept .h5ad or .h5ad.zst.
ROOTS = [a for a in sys.argv[1:] if not a.startswith("-")]
ROOTS += [os.environ.get("MORPHODE_DATA_ROOT", ""),
          os.path.join(REPO, "data")]
SRC = None
for root in ROOTS:
    if not root:
        continue
    for pat in ("zf10_stereoseq.h5ad", "zf10_stereoseq.h5ad.zst",
                "**/zf10_stereoseq.h5ad*", "**/zf10.h5ad", "**/zf10*.h5ad*"):
        hits = sorted(glob.glob(os.path.join(root, pat), recursive=True))
        # prefer a plain .h5ad over a .zst if both present
        hits = [h for h in hits if h.endswith(".h5ad")] + [h for h in hits if h.endswith(".zst")]
        if hits:
            SRC = hits[0]; break
    if SRC:
        break
if SRC is None:
    sys.exit("zf10_stereoseq.h5ad not found. Fetch it with:\n"
             "    python scripts/download.py zesta_zebrafish\n"
             "then re-run, or pass the folder holding it as an argument.")
DST = SRC
if SRC.endswith(".zst"):
    DST = os.path.join(REPO, "data", "_zf10_demo.h5ad")
    os.makedirs(os.path.dirname(DST), exist_ok=True)

t0 = time.time()
print("=" * 56)
print("  morphoDE · zebrafish · Yolk Syncytial Layer (10 hpf)")
print("=" * 56)
if not os.path.exists(DST):
    print("  decompressing dataset ...")
    subprocess.run(["zstd", "-dq", "-f", SRC, "-o", DST], check=True)

a = ad.read_h5ad(DST)
a.obs["annotation"] = a.obs["bin_annotation"].astype(str).values
stage = str(a.obs["time_point"].iloc[0]) if "time_point" in a.obs else "zf10"
a.obs["section_id"] = [f"{stage}_S{int(s):02d}" for s in a.obs["slice"].astype(int).values]
s = assemble_3d([a], "zf10", technology="stereoseq", z_spacing=15.0)
s.obs["section_id"] = a.obs["section_id"].values
s.obs["annotation"] = a.obs["annotation"].values
s.uns["registration"] = "none"
if needs_registration(s):
    s = register_sections_rigid(s)

XYZ = np.asarray(s.obsm["spatial_3d"], float)
is_domain = (s.obs["annotation"].values == DOMAIN)
section = s.obs["section_id"].values
genes = np.asarray(s.var_names)
print(f"  [1] point cloud   : {s.n_obs:,} cells, {s.obs['section_id'].nunique()} sections, "
      f"{int(is_domain.sum()):,} in domain")

V, F, gf = extract_open_interface(XYZ, is_domain)
print(f"  [2] open surface  : {len(V):,} verts, {len(F):,} faces, genuine={gf:.2f}")

band, disp, inb = band_scores(XYZ, is_domain, V, F)
Xn = s.X[band]
print(f"  [3] interface band: {int(inb.sum()):,} cells scored (bulge > 0, dent < 0)")

res = de_and_null(Xn, genes, disp, XYZ[band], section[band], n_perm=200)
print(f"  [4] DE + geometry-breaking null (n=200)")
print("=" * 56)
print(f"  RESULT: {res['n_survive']}/{res['n_tested']} genes survive the null "
      f"(p<=0.05), max|lfc|={res['max_abs_lfc']:.2f}")
print("=" * 56)
print(f"  done in {time.time()-t0:.0f}s")
