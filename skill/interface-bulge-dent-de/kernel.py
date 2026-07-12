"""Helpers for interface bulge-vs-dent differential expression.
See SKILL.md for the workflow. All heavy imports are deferred into function
bodies so loading the skill only defines names.
"""

def compute_bulge_dent_scores(vertices, faces, vertex_normals=None,
                              taubin_lamb=0.7, taubin_nu=-0.72, taubin_iters=60):
    """Two curvature-free bulge/dent scores per vertex.

    Returns dict with:
      hull_dev: signed distance to 3D convex hull (interior<0=dent, ~0=bulge). GLOBAL.
      smooth_disp: outward displacement from strongly-smoothed surface (+=bulge). LOCAL.
      vertex_normals_out: the outward normals actually used.
    Verifies/fixes that normals point outward (a common stored-mesh bug).
    """
    import numpy as np, trimesh
    from scipy.spatial import ConvexHull
    V = np.asarray(vertices, float); F = np.asarray(faces)
    mesh = trimesh.Trimesh(vertices=V, faces=F, process=False)
    mesh.fix_normals()
    VN_out = mesh.vertex_normals
    # convex-hull deviation: max over hull-face signed distances (interior => negative)
    hull = ConvexHull(V)
    eq = hull.equations
    hull_dev = ((V @ eq[:, :3].T) + eq[:, 3]).max(axis=1)
    # smoothed-surface displacement along outward normal
    sm = mesh.copy()
    trimesh.smoothing.filter_taubin(sm, lamb=taubin_lamb, nu=taubin_nu, iterations=taubin_iters)
    disp_vec = V - sm.vertices
    smooth_disp = np.einsum('ij,ij->i', disp_vec, VN_out)
    return {'hull_dev': hull_dev, 'smooth_disp': smooth_disp,
            'vertex_normals_out': VN_out}


def tercile(v):
    """Label a score vector into 'bulge' (top third), 'dent' (bottom third), 'mid'."""
    import numpy as np
    v = np.asarray(v, float)
    q1, q2 = np.quantile(v, [1/3, 2/3])
    return np.where(v <= q1, 'dent', np.where(v >= q2, 'bulge', 'mid'))


def assign_face_score(vertex_score, faces, cell_nearest_face):
    """Map a per-vertex score to cells via nearest-face mean."""
    import numpy as np
    face_score = np.asarray(vertex_score)[np.asarray(faces)].mean(1)
    return face_score[np.asarray(cell_nearest_face)]


def de_meandiff(X, mask_bulge, mask_dent, genes):
    """Mean-log-fold DE (bulge - dent). X: sparse/dense cell x gene (log-normalized)."""
    import numpy as np, scipy.sparse as sp, pandas as pd
    X = X.tocsr() if sp.issparse(X) else sp.csr_matrix(X)
    a = np.asarray(X[mask_bulge].mean(0)).ravel()
    b = np.asarray(X[mask_dent].mean(0)).ravel()
    return pd.DataFrame({'gene': np.asarray(genes), 'mean_bulge': a,
                         'mean_dent': b, 'lfc': a - b}).sort_values('lfc', ascending=False)


def decompose_de(X, labels, celltype, genes, min_cells=20):
    """Symmetric decomposition of total bulge-dent diff into composition + within-cell-type.

    total = Σ_c (p_c^b - p_c^d)·mean(m_c^b,m_c^d)   [composition]
          + Σ_c mean(p_c^b,p_c^d)·(m_c^b - m_c^d)    [within-cell-type]
    labels: array of 'bulge'/'dent'/... ; celltype: per-cell type array.
    """
    import numpy as np, scipy.sparse as sp, pandas as pd
    X = X.tocsc() if sp.issparse(X) else sp.csr_matrix(X).tocsc()
    labels = np.asarray(labels); celltype = np.asarray(celltype)
    mb = labels == 'bulge'; md = labels == 'dent'
    nb = mb.sum(); nd = md.sum()
    cts = [c for c in np.unique(celltype)
           if (mb & (celltype == c)).sum() >= min_cells and (md & (celltype == c)).sum() >= min_cells]
    ng = X.shape[1]
    comp = np.zeros(ng); within = np.zeros(ng)
    for c in cts:
        mbc = mb & (celltype == c); mdc = md & (celltype == c)
        pcb = mbc.sum() / nb; pcd = mdc.sum() / nd
        mcb = np.asarray(X[mbc].mean(0)).ravel(); mcd = np.asarray(X[mdc].mean(0)).ravel()
        comp += (pcb - pcd) * ((mcb + mcd) / 2)
        within += ((pcb + pcd) / 2) * (mcb - mcd)
    return pd.DataFrame({'gene': np.asarray(genes), 'total': comp + within,
                         'composition': comp, 'within_celltype': within})


def azimuth_order(cells):
    """Per-section vertex order by azimuth about the section centroid. Returns dict."""
    import numpy as np
    sec = cells['section'].values
    x = cells['x'].values; y = cells['y'].values
    order = {}
    for s in np.unique(sec):
        m = np.where(sec == s)[0]
        cx, cy = x[m].mean(), y[m].mean()
        az = np.arctan2(y[m] - cy, x[m] - cx)
        order[s] = (m, m[np.argsort(az)])
    return order


def azimuthal_shift_null(cells, X, score_col, target_genes, genes, n_perm=200, seed=0):
    """Geometry-destroying null: per-section azimuthal circular shift of the bulge/dent score.

    Preserves spatial smoothness, destroys score<->position correspondence.
    Returns DataFrame gene, obs_lfc, null_sd, p_geomnull.
    """
    import numpy as np, scipy.sparse as sp, pandas as pd
    rng = np.random.default_rng(seed)
    X = X.tocsc() if sp.issparse(X) else sp.csr_matrix(X).tocsc()
    genes = np.asarray(genes)
    gidx = {g: i for i, g in enumerate(genes)}
    tg = [g for g in target_genes if g in gidx]
    cols = np.array([gidx[g] for g in tg])
    Xtg = np.asarray(X[:, cols].todense())
    score = cells[score_col].values
    order = azimuth_order(cells)
    # observed
    q1, q2 = np.quantile(score, [1/3, 2/3])
    mb = score >= q2; md = score <= q1
    obs = Xtg[mb].mean(0) - Xtg[md].mean(0)
    null = np.zeros((n_perm, len(tg)))
    n = len(cells)
    for p in range(n_perm):
        sh = np.empty(n)
        for s, (m, o) in order.items():
            vals = score[o]
            off = rng.integers(1, len(vals)) if len(vals) > 2 else 0
            sh[o] = np.roll(vals, off)
        qq1, qq2 = np.quantile(sh, [1/3, 2/3])
        null[p] = Xtg[sh >= qq2].mean(0) - Xtg[sh <= qq1].mean(0)
    rows = []
    for i, g in enumerate(tg):
        pv = (np.sum(np.abs(null[:, i]) >= abs(obs[i])) + 1) / (n_perm + 1)
        rows.append({'gene': g, 'obs_lfc': obs[i], 'null_sd': null[:, i].std(), 'p_geomnull': pv})
    return pd.DataFrame(rows)


def null_calibration(X, cells, score_col, n_genes=100, n_perm=200, seed=7):
    """Shuffle expression across cells (true signal=0); geometry-null p<0.05 should be ~nominal."""
    import numpy as np, scipy.sparse as sp
    rng = np.random.default_rng(seed)
    X = X.tocsc() if sp.issparse(X) else sp.csr_matrix(X).tocsc()
    score = cells[score_col].values
    order = azimuth_order(cells)
    n = len(cells)
    q1, q2 = np.quantile(score, [1/3, 2/3])
    mb = score >= q2; md = score <= q1
    # precompute permutation masks
    masks = []
    for _ in range(n_perm):
        sh = np.empty(n)
        for s, (m, o) in order.items():
            vals = score[o]; off = rng.integers(1, len(vals)) if len(vals) > 2 else 0
            sh[o] = np.roll(vals, off)
        qq1, qq2 = np.quantile(sh, [1/3, 2/3])
        masks.append((sh >= qq2, sh <= qq1))
    samp = rng.choice(X.shape[1], min(n_genes, X.shape[1]), replace=False)
    Xs = np.asarray(X[:, samp].todense())
    pvals = []
    for gi in range(Xs.shape[1]):
        x = rng.permutation(Xs[:, gi])
        obs = x[mb].mean() - x[md].mean()
        nl = np.array([x[nb].mean() - x[nd].mean() for nb, nd in masks])
        pvals.append((np.sum(np.abs(nl) >= abs(obs)) + 1) / (n_perm + 1))
    pvals = np.array(pvals)
    return {'frac_p05': float((pvals < 0.05).mean()), 'frac_p10': float((pvals < 0.10).mean()),
            'median_p': float(np.median(pvals)), 'pvals': pvals}


def power_curve(X, cells, score_col, effects=None, n_rep=100, n_perm=200, seed=11):
    """Spike a known log-fold into the bulge side; return detection rate vs effect size."""
    import numpy as np, scipy.sparse as sp, pandas as pd
    if effects is None:
        effects = [0.0, 0.02, 0.05, 0.10, 0.15, 0.20, 0.30]
    rng = np.random.default_rng(seed)
    X = X.tocsc() if sp.issparse(X) else sp.csr_matrix(X).tocsc()
    score = cells[score_col].values
    order = azimuth_order(cells)
    n = len(cells)
    q1, q2 = np.quantile(score, [1/3, 2/3])
    mb = score >= q2; md = score <= q1
    masks = []
    for _ in range(n_perm):
        sh = np.empty(n)
        for s, (m, o) in order.items():
            vals = score[o]; off = rng.integers(1, len(vals)) if len(vals) > 2 else 0
            sh[o] = np.roll(vals, off)
        qq1, qq2 = np.quantile(sh, [1/3, 2/3])
        masks.append((sh >= qq2, sh <= qq1))
    ncol = X.shape[1]
    pool = np.asarray(X[:, rng.choice(ncol, min(200, ncol), replace=False)].todense())
    out = []
    for eff in effects:
        det = 0
        for r in range(n_rep):
            x = rng.permutation(pool[:, rng.integers(0, pool.shape[1])]).astype(float)
            x[mb] += eff
            obs = x[mb].mean() - x[md].mean()
            nl = np.array([x[nb].mean() - x[nd].mean() for nb, nd in masks])
            p = (np.sum(np.abs(nl) >= abs(obs)) + 1) / (n_perm + 1)
            if p < 0.05:
                det += 1
        out.append({'effect_lfc': eff, 'power': det / n_rep})
    return pd.DataFrame(out)
