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
    nd = _norm(domain)
    stem = nd[:6] if len(nd) >= 6 else nd
    keep = set()
    for v in set(ann):
        nv = _norm(v)
        if nd == nv or nd in nv or nv in nd or stem in nv:
            keep.add(v)
    if not keep:                              # fall back to exact label
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
    if loader == "adapter":
        import st3d_adapters as A
        parts = A.get(dataset_id)(os.path.join(data_root, dataset_id))
        parts = parts if isinstance(parts, list) else [parts]
        s = L.load(dataset_id, data_root=data_root, workdir=workdir,
                   technology=technology, z_spacing=z_spacing,
                   auto_register=False, parts=parts)
    else:
        s = L.load(dataset_id, data_root=data_root, workdir=workdir,
                   technology=technology, z_spacing=z_spacing, auto_register=False)
    return s


def ensure_lognorm(adata):
    """de_and_null expects log-normalized expression. Datasets ship in mixed
    states (some pre-normalized, some raw counts); normalize in place only when
    X looks like raw integer counts, so a pre-normalized matrix is left alone."""
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
    import numpy as np
    from open_interface import extract_open_interface, band_scores, de_and_null
    strat = stratum_mask(s, dataset_id, stratum)
    if strat.sum() == 0:
        return {"skip": f"stratum '{stratum}' matched 0 cells"}
    sub = s[strat]
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
        t0 = time.time()
        try:
            s = load_dataset(d, meta["loader"], args.data_root,
                             meta["technology"], float(meta["z_spacing"]))
        except Exception as e:
            print(f"  LOAD FAILED: {type(e).__name__}: {str(e)[:160]}")
            continue
        s = ensure_lognorm(s)
        norm = s.uns.get("_morphode_normalized", "already log-normalized")
        print(f"  loaded {s.n_obs:,} cells, {s.obs['section_id'].nunique()} sections "
              f"in {time.time()-t0:.0f}s  [{norm}]")
        out = []
        for r in rows:
            res = run_interface(s, d, r["domain"], r["stratum"], n_perm=args.n_perm)
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
            print(f"  wrote {fp} ({len(out)} rows)")

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
