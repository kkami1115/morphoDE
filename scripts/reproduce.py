#!/usr/bin/env python3
"""morphoDE — end-to-end reproduction driver.

Runs the published open-interface bulge/dent method on any of the 11 analysis
datasets, straight from this repository. Given a decompressed dataset in the
data root, for each interface-stratum in data/recipes.csv it:

    load 3-D point cloud (st3d_loader / st3d_adapters)
      -> genuine open-interface sheet         (open_interface.extract_open_interface)
      -> bulge/dent band scores               (open_interface.band_scores)
      -> within-cell-type DE + geometry null  (open_interface.de_and_null)
      -> row appended to results/<dataset>_bd.csv

The (dataset, interface, stratum) rows in data/recipes.csv ARE the analysis
recipe: which cell-type domain to raise an interface against, and which serial
sub-block (developmental stage / regeneration time / embryo replicate) to run
it within. recipes.csv also carries the published headline counts
(truth_survive / truth_tested) so --check can compare a fresh run to them.

Data is NOT downloaded here — use scripts/download.py first, or point --data-root
at an existing decompressed mirror (e.g. a local cache of .h5ad files).

Examples
--------
    python scripts/reproduce.py --list
    python scripts/reproduce.py --smoke                 # 1 light dataset, fast
    python scripts/reproduce.py --dataset zesta_zebrafish
    python scripts/reproduce.py --all --check
"""
import argparse, os, re, sys, time, json, csv

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
SRC  = os.path.join(REPO, "src")
sys.path.insert(0, SRC)

RECIPES = os.path.join(REPO, "data", "recipes.csv")
CATALOG = os.path.join(REPO, "data", "catalog.csv")
# Default data root: repo-local data/, but honour a cache via env or --data-root.
DEFAULT_ROOT = os.environ.get("MORPHODE_DATA_ROOT", os.path.join(REPO, "data"))

# --- light datasets suitable for a fast smoke test (small, single-stratum) ---
SMOKE_DATASETS = ["acsta_arabidopsis", "zesta_zebrafish"]

# ----------------------------------------------------------------------------
# Per-dataset stratum-key extractor: maps a source-file token to the recipe's
# `stratum` value. Datasets whose recipes use stratum='all' need no entry (the
# whole assembled object is one stratum). The four stratified datasets encode
# their sub-block in the source filename; these rules reproduce that mapping.
# ----------------------------------------------------------------------------
def _dm_key(fname):
    m = re.search(r'(E[\d.]+)_stereo_(rep\d)', fname)          # GSM..._Embryo_E7.5_stereo_rep1
    return (m.group(1) + m.group(2)).replace(".", "p") if m else None

STRATUM_KEY = {
    "flysta3d_drosophila":        lambda f: f.split("_")[0],                       # E14-16h, L1, L2, L3
    "prista4d_planarian":         lambda f: f.split(".")[0],                       # 0hpa1, 12hpa2, 3dpa1
    "artista_axolotl_brain":      lambda f: (re.match(r'(\d+DPI)', f) or [None,None])[0]
                                            if re.match(r'(\d+DPI)', f) else None, # 2DPI..20DPI (no 15/25 bleed)
    "digital_mouse_embryo_seu3d": _dm_key,                                         # E7p5rep1..E8p0rep2
}


def _norm(x):
    return re.sub(r'[^a-z0-9]', '', str(x).lower().replace('.', 'p'))


def domain_mask(adata, domain):
    """Cells whose annotation belongs to the interface domain.

    Reproduces st3d_interface._vals_present: substring match of the domain
    keyword (stemmed) against the dataset's own annotation vocabulary, so
    'epidermis' picks up 'Upper_epidermal_cell'/'Lower_epidermal_cell', while a
    single-token domain like 'Notochord' matches exactly.
    """
    import numpy as np
    ann = adata.obs["annotation"].astype(str).values
    vocab = set(ann)
    # 0) Explicit multi-label domain: recipes may list several exact labels joined
    #    by '+', for a domain the published run matched as a union of distinct
    #    annotation classes (e.g. openst tumor = 'Tumor' + 'Tumor_Keratin_Pearl',
    #    which the headline picked up via the 'umor' substring). Match each exactly.
    if "+" in domain:
        want = {_norm(p) for p in domain.split("+")}
        keep = {v for v in vocab if _norm(v) in want}
        if keep:
            return np.isin(ann, list(keep)), sorted(keep)
    # 1) Exact label (case-insensitive) — the published runs used np.isin(ann,[domain])
    #    whenever the domain name was itself a label. This is important where the
    #    vocabulary has *compound* labels ("epidermis/CNS", "fat body/trachea",
    #    "midgut/malpighian tubules"): those are distinct cell classes and must NOT
    #    be folded into the base domain, or the interface geometry/DE shifts.
    exact = {v for v in vocab if _norm(v) == _norm(domain)}
    if exact:
        return np.isin(ann, list(exact)), sorted(exact)
    # 2) No exact label (e.g. acsta 'epidermis' -> 'Upper_epidermal_cell'): fall back
    #    to a stemmed-substring match against the vocabulary.
    nd = _norm(domain)
    stem = nd[:6] if len(nd) >= 6 else nd
    keep = {v for v in vocab if nd in _norm(v) or _norm(v) in nd or stem in _norm(v)}
    if not keep:
        keep = {domain}
    return np.isin(ann, list(keep)), sorted(keep)


def stratum_mask(adata, dataset_id, stratum):
    """Boolean mask of cells belonging to a recipe stratum."""
    import numpy as np
    n = adata.n_obs
    if str(stratum) == "all":
        return np.ones(n, bool)
    keyfn = STRATUM_KEY.get(dataset_id)
    src_col = "_srcfile" if "_srcfile" in adata.obs else "section_id"
    src = adata.obs[src_col].astype(str).values
    if keyfn is not None:
        keys = np.array([(_norm(keyfn(f)) if keyfn(f) else "") for f in src])
        return keys == _norm(stratum)
    # generic: normalized substring on the source token
    ns = _norm(stratum)
    return np.array([ns in _norm(f) for f in src])


def load_dataset(dataset_id, loader, data_root, technology, z_spacing):
    """Assemble the 3-D SpatialSample for a dataset via adapter or generic loader."""
    import st3d_loader as L
    workdir = os.path.join(os.environ.get("MORPHODE_WORK",
                           os.path.join(data_root, "_work")), dataset_id)
    if dataset_id == "openst_lymphnode_3d":
        return load_openst(data_root)
    if dataset_id == "artista_axolotl_brain":
        return load_artista(data_root)
    if dataset_id == "flysta3d_v2_drosophila":
        return load_flysta3d_v2(data_root)
    if dataset_id == "mosta_mouse_embryo":
        return load_mosta(data_root)
    if dataset_id == "whole_mouse_embryo_3d_cngb":
        return load_whole_mouse(data_root)
    if loader == "adapter":
        import st3d_adapters as A
        # The registered adapter takes the data ROOT (it appends dataset_id/... itself)
        # and returns a list of per-section AnnData parts; the generic loader then
        # stacks them. This mirrors the published headline call
        #   parts = REGISTRY[dataset_id](root=data_root); load(dataset_id, parts=parts)
        parts = A.REGISTRY[dataset_id](root=data_root)
        parts = parts if isinstance(parts, list) else [parts]
        s = L.load(dataset_id, data_root=data_root, workdir=workdir,
                   technology=technology, z_spacing=z_spacing, parts=parts)
    else:
        # auto_register=True (loader default) rigidly aligns serial sections before
        # the interface mesh is built — this is what the published headline used.
        # Disabling it leaves sections un-aligned, which distorts the 3-D coordinates
        # and the extracted interface geometry (verts/genuine_frac), breaking the
        # numeric match. Keep it ON for the generic multi-section loaders.
        s = L.load(dataset_id, data_root=data_root, workdir=workdir,
                   technology=technology, z_spacing=z_spacing, auto_register=True)
    return s


def load_single_file(dataset_id, data_root, fname):
    """Assemble a 3-D sample from ONE source file (per-embryo / per-specimen).

    Some datasets analyse each interface within a single representative specimen
    rather than the whole stacked object (e.g. zebrafish: each domain uses one
    embryo/timepoint). The recipe row names that file in `source_file`; we glob
    for it under the dataset dir (raw .h5ad or .h5ad.zst) and assemble it, rigidly
    registering sections when the file holds more than one.
    """
    import glob, st3d_loader as L
    ddir = os.path.join(data_root, dataset_id)
    # Prefer a plain .h5ad (canonical download / decompressed cache); fall back to
    # a .zst original, which st3d_loader transparently decompresses.
    hits = glob.glob(os.path.join(ddir, "**", fname), recursive=True)
    if not hits:
        work = os.path.join(os.environ.get("MORPHODE_WORK",
                            os.path.join(data_root, "_work")), dataset_id)
        hits = glob.glob(os.path.join(work, "**", fname), recursive=True)
    if not hits:
        hits = glob.glob(os.path.join(ddir, "**", fname + ".zst"), recursive=True)
    if not hits:
        raise FileNotFoundError(f"{fname} not found under {ddir}")
    a = L.assemble_3d([hits[0]], dataset_id)
    if a.obs["section_id"].nunique() > 1:
        a = L.register_sections_rigid(a)
    return a


#  The 24 axolotl-brain serial sections, in the published stacking order. Each
#  is one .h5ad under artista_axolotl_brain/(stomics/); the whole-dir glob would
#  also pick up pre-merged "10DPIs.h5ad" files and double-count, so we name the
#  per-section files explicitly and stack them with z = section index (no rigid
#  registration — the published run used raw section order).
ARTISTA_SECTIONS = ["2DPI_1", "2DPI_2", "2DPI_3", "5DPI_1", "5DPI_2", "5DPI_3",
                    "10DPI_1", "10DPI_2", "10DPI_3", "15DPI_1", "15DPI_2", "15DPI_3",
                    "15DPI_4", "20DPI_1", "20DPI_2", "20DPI_3", "30DPI", "60DPI",
                    "Control_Juv", "Adult", "Meta", "Stage44", "Stage54", "Stage57"]


def load_artista(data_root):
    """Bespoke assembly for artista_axolotl_brain matching the published run.

    Reads the 24 named serial sections, takes 2-D `spatial` + the `Annotation`
    obs column from each, and stacks them with z = running section index. The
    `_srcfile` column carries the section name so the recipe's DPI stratum
    (2DPI..20DPI) selects the right sections via STRATUM_KEY.
    """
    import glob, numpy as np, scipy.sparse as sp, anndata as ad, st3d_loader as L
    ddir = os.path.join(data_root, "artista_axolotl_brain")
    work = os.path.join(os.environ.get("MORPHODE_WORK",
                        os.path.join(data_root, "_work")), "artista_axolotl_brain")
    probe = os.path.join(os.environ.get("MORPHODE_WORK",
                         os.path.join(data_root, "_work")), "artista_probe")
    parts, gs = [], 0
    for name in ARTISTA_SECTIONS:
        # exact filename only (avoid the pre-merged "<DPI>s.h5ad" aggregates)
        cand = []
        for base in (ddir, work, probe):
            cand += glob.glob(os.path.join(base, "**", name + ".h5ad"), recursive=True)
        if not cand:
            for base in (ddir,):
                cand += glob.glob(os.path.join(base, "**", name + ".h5ad.zst"), recursive=True)
        cand = [c for c in cand if os.path.basename(c).split(".")[0] == name]
        if not cand:
            continue
        fp = cand[0]
        if fp.endswith(".zst"):
            import zstandard, io
            with open(fp, "rb") as fh:
                a = ad.read_h5ad(io.BytesIO(zstandard.ZstdDecompressor().decompress(fh.read())))
        else:
            a = ad.read_h5ad(fp)
        X = (a.X.tocsr() if sp.issparse(a.X) else sp.csr_matrix(a.X)).astype(np.float32)
        annot = a.obs["Annotation"].astype(str).values if "Annotation" in a.obs \
                else a.obs["annotation"].astype(str).values
        b = ad.AnnData(X=X, obs=dict(annotation=annot,
                       section_id=np.array([name] * a.n_obs),
                       _srcfile=np.array([name] * a.n_obs),
                       _secidx=np.full(a.n_obs, gs, int)), var=a.var[[]].copy())
        b.var_names = a.var_names
        b.obsm["spatial2d"] = np.asarray(a.obsm["spatial"], float)[:, :2]
        parts.append(b); gs += 1
    if not parts:
        raise FileNotFoundError(f"no artista sections found under {ddir}")
    ax = ad.concat(parts, join="outer", label="_concat_key",
                   keys=[p.obs["section_id"][0] for p in parts], index_unique="-")
    if sp.issparse(ax.X):
        ax.X.data = np.nan_to_num(ax.X.data, nan=0.0)
    ax.obsm["spatial_3d"] = np.column_stack(
        [ax.obsm["spatial2d"], ax.obs["_secidx"].values.astype(float)]).astype(float)
    ax.X = ax.X.tocsr()
    return ax


#  flysta3d_v2: the published run analysed ONE pupal embryo (PUPA-72h_cellbin),
#  manually library-size normalised it, and mapped the fine cell-type annotation
#  onto four coarse tissue domains with composite rules. Reproduced verbatim.
FLYV2_EMBRYO = "PUPA-72h_cellbin_final_modified.h5ad"

def _flyv2_domain(lab):
    l = str(lab).lower()
    if ("brain" in l) or ("nerve cord" in l) or (l == "glia"):            return "CNS primordium"
    if ("midgut" in l) or ("gastric caecum" in l) or ("proventriculus" in l): return "midgut primordium"
    if "epidermis" in l:                                                   return "epidermis"
    if l == "somatic muscle":                                             return "somatic muscle"
    return "rest"


def load_flysta3d_v2(data_root):
    """Single-embryo bespoke assembly + manual lib-size norm + coarse domains."""
    import glob, io, numpy as np, scipy.sparse as sp, anndata as ad, st3d_loader as L
    ddir = os.path.join(data_root, "flysta3d_v2_drosophila")
    work = os.path.join(os.environ.get("MORPHODE_WORK",
                        os.path.join(data_root, "_work")), "flysta3d_v2_drosophila")
    cand = []
    for base in (ddir, work):
        cand += glob.glob(os.path.join(base, "**", FLYV2_EMBRYO), recursive=True)
    if not cand:
        cand += glob.glob(os.path.join(ddir, "**", FLYV2_EMBRYO + ".zst"), recursive=True)
    if not cand:
        raise FileNotFoundError(f"{FLYV2_EMBRYO} not found under {ddir}")
    fp = cand[0]
    fv = L.assemble_3d([fp], "flysta3d_v2_drosophila")
    fv.uns["registration"] = "none"
    # manual library-size normalization + log1p (headline cell303)
    X = fv.X.tocsr().astype(np.float32)
    lib = np.asarray(X.sum(1)).ravel(); lib[lib == 0] = 1.0
    X = sp.diags((1e4 / lib).astype(np.float32)) @ X
    X.data = np.log1p(X.data); fv.X = X.tocsr()
    fv.uns["_morphode_normalized"] = "manual lib-size 1e4 + log1p"
    fine = fv.obs["annotation"].astype(str).values
    import pandas as pd
    fv.obs["annotation"] = pd.Categorical([_flyv2_domain(l) for l in fine])
    return fv


def load_openst(data_root):
    """openST lymph node: one file that already carries an ALIGNED 3-D embedding
    in obsm['spatial_3d_aligned']. Use it directly as spatial_3d (the published
    run did); the generic loader would instead rebuild z from section order.
    """
    import glob, io, numpy as np, scipy.sparse as sp, anndata as ad
    fname = "GSE251926_metastatic_lymph_node_3d.h5ad"
    ddir = os.path.join(data_root, "openst_lymphnode_3d")
    work = os.path.join(os.environ.get("MORPHODE_WORK",
                        os.path.join(data_root, "_work")), "openst_lymphnode_3d")
    cand = []
    for base in (ddir, work):
        cand += glob.glob(os.path.join(base, "**", fname), recursive=True)
    if not cand:
        cand += glob.glob(os.path.join(ddir, "**", fname + ".zst"), recursive=True)
    if not cand:
        raise FileNotFoundError(f"{fname} not found under {ddir}")
    fp = cand[0]
    if fp.endswith(".zst"):
        import zstandard
        with open(fp, "rb") as fh:
            ao = ad.read_h5ad(io.BytesIO(zstandard.ZstdDecompressor().decompress(fh.read())))
    else:
        ao = ad.read_h5ad(fp)
    ao.obsm["spatial_3d"] = np.asarray(ao.obsm["spatial_3d_aligned"], float)
    if "section_id" not in ao.obs:
        ao.obs["section_id"] = np.array(["s0"] * ao.n_obs)
    ao.X = ao.X.tocsr() if sp.issparse(ao.X) else sp.csr_matrix(ao.X)
    return ao


def load_mosta(data_root):
    """MOSTA E16.5_E2 serial sections stacked with z = section number, NO rigid
    registration (the published run concatenated the raw sections and used the
    section index as z). The registered generic path distorts the coordinates and
    breaks the geometry, so we assemble directly here.
    """
    import re, numpy as np, scipy.sparse as sp, anndata as ad
    import st3d_adapters as A
    parts = A.REGISTRY["mosta_mouse_embryo"](root=data_root)
    xyz = []
    for i, p in enumerate(parts):
        xy = np.asarray(p.obsm["spatial"], float)[:, :2]
        sec = str(p.obs["section"].iloc[0]) if "section" in p.obs else str(i)
        m = re.search(r"S(\d+)", sec)
        z = float(m.group(1)) if m else float(i)
        p.obs["section_id"] = np.array([f"S{int(z):02d}"] * p.n_obs)
        xyz.append(np.column_stack([xy, np.full(p.n_obs, z)]))
    mo = ad.concat(parts, join="outer", label="_srcfile", index_unique="-")
    mo.obsm["spatial_3d"] = np.vstack(xyz).astype(float)
    mo.X = mo.X.tocsr() if sp.issparse(mo.X) else sp.csr_matrix(mo.X)
    return mo


def load_whole_mouse(data_root):
    """CNGB whole-mouse-embryo: the generic loader stacks + registers the parts,
    then the published run OVERWRITES z with (section_number - 1) so the sections
    sit at integer depths. Reproduced here.
    """
    import re, numpy as np, scipy.sparse as sp, st3d_loader as L
    import st3d_adapters as A
    parts = A.REGISTRY["whole_mouse_embryo_3d_cngb"](root=data_root)
    cn = L.load("whole_mouse_embryo_3d_cngb", data_root=data_root, parts=parts)
    sid = cn.obs["section_id"].astype(str).values
    num = np.array([int(re.sub(r"\D", "", s) or 0) for s in sid])
    newz = (num - 1).astype(float)
    cn.obs["z"] = newz
    xyz = np.asarray(cn.obsm["spatial_3d"], float); xyz[:, 2] = newz
    cn.obsm["spatial_3d"] = xyz
    cn.X = cn.X.tocsr() if sp.issparse(cn.X) else sp.csr_matrix(cn.X)
    return cn



def ensure_lognorm(adata):
    """de_and_null expects log-normalized expression. The datasets ship in mixed
    states: some are stored already log-normalized (digital_mouse, flysta3d,
    zesta, ...), others as raw integer counts (acsta, whole_mouse, cerebellum,
    prista4d). Normalize (library-size 1e4 + log1p) ONLY when X looks like raw
    counts, so an already-log matrix is left untouched.

    Note on zesta: the published run passed zesta through run_open_stratified
    with is_raw=True, which re-normalized an already-log matrix (a double
    normalization). We deliberately do NOT reproduce that — zesta's stored X is
    already log (max ~7, non-integer), so it is left as-is. This is the one place
    reproduce.py is scientifically correct rather than bit-identical to the
    headline: it changes only the zesta Segmental-Plate survivor count (8 vs the
    published 10); the other three zesta interfaces are unaffected.
    """
    import numpy as np, scipy.sparse as sp
    X = adata.X
    sample = (X[:2000].toarray() if sp.issparse(X) else np.asarray(X[:2000]))
    if sample.size and float(sample.max()) > 30 and np.mean(sample == np.round(sample)) > 0.9:
        import scanpy as sc
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        adata.uns["_morphode_normalized"] = "normalize_total(1e4)+log1p"
    return adata


def run_interface(s, dataset_id, domain, stratum, n_perm=200):
    """One interface-stratum through the published method. Returns a result dict."""
    import numpy as np, scipy.sparse as sp
    from open_interface import extract_open_interface, band_scores, de_and_null
    strat = stratum_mask(s, dataset_id, stratum)
    if strat.sum() == 0:
        return {"skip": f"stratum '{stratum}' matched 0 cells"}
    # Avoid a full copy for whole-object strata (matters at millions of cells):
    # slice only when the stratum is a real subset.
    sub = s if bool(strat.all()) else s[strat]
    XYZ = np.asarray(sub.obsm["spatial_3d"], float)
    is_dom, dom_labels = domain_mask(sub, domain)
    if is_dom.sum() < 20 or (~is_dom).sum() < 20:
        return {"skip": f"domain '{domain}' -> {int(is_dom.sum())} cells (need >=20 each side)"}
    section = sub.obs["section_id"].astype(str).values
    genes = np.asarray(sub.var_names)
    try:
        V, F, gf = extract_open_interface(XYZ, is_dom)
        band, disp, inb = band_scores(XYZ, is_dom, V, F)
        if int(inb.sum()) < 20:
            return {"skip": f"interface band too small ({int(inb.sum())} cells)"}
        Xn = sub.X[band]
        res = de_and_null(Xn, genes, disp, XYZ[band], section[band], n_perm=n_perm)
    except Exception as e:
        return {"skip": f"geometry/DE failed: {type(e).__name__}: {str(e)[:120]}"}
    return {"dataset": dataset_id, "domain": domain, "stratum": stratum,
            "n_stratum_cells": int(strat.sum()), "n_domain_cells": int(is_dom.sum()),
            "n_band": int(inb.sum()), "genuine_frac": round(float(gf), 4),
            "n_survive": int(res["n_survive"]), "n_tested": int(res["n_tested"]),
            "max_abs_lfc": round(float(res["max_abs_lfc"]), 4),
            "top_genes": ";".join(map(str, res["genes"][:6]))}


def load_recipes():
    with open(RECIPES) as fh:
        return list(csv.DictReader(fh))


def main():
    ap = argparse.ArgumentParser(description="morphoDE end-to-end reproduction")
    ap.add_argument("--dataset", help="single dataset_id to run")
    ap.add_argument("--all", action="store_true", help="run every dataset in recipes.csv")
    ap.add_argument("--smoke", action="store_true",
                    help="fast E2E check: 1 light dataset, first interface only")
    ap.add_argument("--list", action="store_true", help="list datasets and interface counts")
    ap.add_argument("--check", action="store_true",
                    help="compare fresh survive-counts against published headline")
    ap.add_argument("--data-root", default=DEFAULT_ROOT,
                    help=f"root holding <dataset_id>/ dirs (default: {DEFAULT_ROOT})")
    ap.add_argument("--n-perm", type=int, default=200, help="null permutations (default 200)")
    ap.add_argument("--outdir", default=os.path.join(REPO, "results"))
    args = ap.parse_args()

    rec = load_recipes()
    datasets = sorted({r["dataset"] for r in rec})

    if args.list:
        print(f"{len(rec)} interface-strata across {len(datasets)} datasets\n")
        for d in datasets:
            g = [r for r in rec if r["dataset"] == d]
            strata = sorted({r["stratum"] for r in g})
            loader = g[0]["loader"]
            print(f"  {d:32s} {len(g):>3} interface-strata  loader={loader:7s} "
                  f"strata={len(strata)}")
        return

    if args.smoke:
        targets = [d for d in SMOKE_DATASETS if d in datasets][:1] or datasets[:1]
        smoke_one = True
    elif args.dataset:
        targets = [args.dataset]
        smoke_one = False
    elif args.all:
        targets = datasets
        smoke_one = False
    else:
        ap.error("choose one of --list / --smoke / --dataset ID / --all")

    os.makedirs(args.outdir, exist_ok=True)
    grand = []
    for d in targets:
        rows = [r for r in rec if r["dataset"] == d]
        if smoke_one:
            rows = rows[:1]
        meta = rows[0]
        print(f"\n=== {d}  ({len(rows)} interface-strata, loader={meta['loader']}) ===")
        # Group by source_file: rows naming a per-specimen file load that file
        # individually; the rest share one whole-dataset load (key "").
        by_src = {}
        for r in rows:
            by_src.setdefault(r.get("source_file", "") or "", []).append(r)
        cache = {}   # source_file -> (sample, norm_tag)  loaded lazily & reused
        out = []
        for r in rows:
            key = r.get("source_file", "") or ""
            if key not in cache:
                t0 = time.time()
                try:
                    if key:
                        s = load_single_file(d, args.data_root, key)
                    else:
                        s = load_dataset(d, meta["loader"], args.data_root,
                                         meta["technology"], float(meta["z_spacing"]))
                except Exception as e:
                    print(f"  LOAD FAILED [{key or 'whole'}]: {type(e).__name__}: {str(e)[:150]}")
                    cache[key] = (None, None)
                    continue
                s = ensure_lognorm(s)
                norm = s.uns.get("_morphode_normalized", "already log-normalized")
                print(f"  loaded [{key or 'whole'}] {s.n_obs:,} cells, "
                      f"{s.obs['section_id'].nunique()} sections in {time.time()-t0:.0f}s  [{norm}]")
                cache[key] = (s, norm)
            s, norm = cache[key]
            if s is None:
                continue
            res = run_interface(s, d, r["domain"], r["stratum"], n_perm=args.n_perm)
            import gc; gc.collect()   # release per-interface geometry before the next
            if "skip" in res:
                print(f"  - {r['domain']:26s} [{r['stratum']:>8}]  SKIP: {res['skip']}")
                continue
            tag = ""
            if args.check:
                tag = f"  (headline {r['truth_survive']}/{r['truth_tested']})"
            print(f"  - {r['domain']:26s} [{r['stratum']:>8}]  "
                  f"survive={res['n_survive']}/{res['n_tested']} "
                  f"genuine={res['genuine_frac']}{tag}")
            res["truth_survive"] = int(r["truth_survive"])
            res["truth_tested"]  = int(r["truth_tested"])
            out.append(res)
            grand.append(res)
        if out:
            fp = os.path.join(args.outdir, f"{d}_bd.csv")
            with open(fp, "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=list(out[0].keys()))
                w.writeheader(); w.writerows(out)
            try:
                rel = os.path.relpath(fp)
            except ValueError:
                rel = fp
            print(f"  wrote {rel} ({len(out)} rows)")

    if grand:
        total = sum(r["n_survive"] for r in grand)
        tested = sum(r["n_tested"] for r in grand)
        print(f"\nTOTAL: {len(grand)} interface-strata run; "
              f"{total}/{tested} gene-tests survive the null.")
        if args.check:
            match = sum(1 for r in grand if r["n_survive"] == r["truth_survive"])
            print(f"CHECK: {match}/{len(grand)} interface-strata match published "
                  f"survive-count exactly.")


if __name__ == "__main__":
    main()
