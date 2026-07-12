"""Corrected bulge/dent on the GENUINE open tumor|domain interface.

Fixes the closed-solid convex-hull flaw: instead of closing the domain into a
watertight solid and measuring convex-hull deviation, we extract the actual
contact sheet between the domain and neighbouring SAMPLED tissue, drop faces
that border unsampled space (block cut faces / z-caps), and define bulge/dent
as local normal displacement of that open sheet.
"""
from __future__ import annotations
import numpy as np, scipy.sparse as sp
from scipy import ndimage
from scipy.spatial import cKDTree
from skimage import measure
import trimesh


def extract_open_interface(XYZ, is_domain, pitch=None, tissue_dilate=1, smooth=1.0,
                           min_axis_voxels=4):
    """Return (Vsheet, Fsheet, frac_genuine) — the open genuine-interface sheet.

    XYZ         (N,3) all-cell coords; is_domain (N,) bool for the target domain.
    Keeps only marching-cubes faces with SAMPLED tissue on both sides.

    Uses an ANISOTROPIC per-axis pitch: the in-plane pitch is span_max/90
    (or the scalar `pitch` if given), but each axis is capped so it yields at
    least `min_axis_voxels` voxels. This lets thin-slab data (z spanning only a
    few serial sections) still form a 3D iso-surface with the same open-sheet
    method used for thick tissue — one unified 3D pipeline across scales.
    """
    lo = XYZ.min(0); hi = XYZ.max(0); span = hi - lo
    if pitch is None:
        # density-adaptive base pitch: voxel edge sized so each voxel holds
        # ~target_cells_per_voxel SAMPLED cells, using the ALL-cell density in the
        # bounding box (the sheet is an iso-surface of the tumor-fraction field over
        # all tissue). Dense data -> fine sheet; sparse data -> coarse sheet, with no
        # per-dataset hand tuning. Bounded by a COST cap (total voxels <= max_voxels)
        # on the fine side and span_max/12 on the coarse side (usable 3D volume).
        target_cells_per_voxel = 4.0
        max_voxels = 6.0e6
        n_cells = max(int(XYZ.shape[0]), 1)
        vol = float(np.prod(span[span > 0]))
        dens_pitch = (vol / n_cells * target_cells_per_voxel) ** (1.0 / 3.0)
        pitch_cost_floor = (vol / max_voxels) ** (1.0 / 3.0)   # finer than this blows voxel budget
        base = float(np.clip(dens_pitch, pitch_cost_floor, np.max(span) / 12.0))
    else:
        base = float(pitch)
    # per-axis pitch: never coarser than what gives min_axis_voxels along that axis
    pitch3 = np.minimum(base, span / max(min_axis_voxels, 2))
    pitch3 = np.where(span > 0, pitch3, base)
    dims = np.ceil(span / pitch3).astype(int) + 1
    def vidx(P):
        return np.clip(((P - lo) / pitch3).astype(int), 0, dims - 1)
    vi = vidx(XYZ)
    n_dom = np.zeros(dims); n_all = np.zeros(dims)
    np.add.at(n_dom, (vi[is_domain,0], vi[is_domain,1], vi[is_domain,2]), 1)
    np.add.at(n_all, (vi[:,0], vi[:,1], vi[:,2]), 1)
    tissue = n_all > 0
    frac = np.where(tissue, n_dom / np.maximum(n_all, 1), np.nan)
    # fill non-tissue voxels with nearest tissue value so the iso-surface only
    # forms inside sampled tissue (not against empty space)
    idx = ndimage.distance_transform_edt(~tissue, return_distances=False, return_indices=True)
    frac = frac[tuple(idx)]
    if smooth > 0:
        frac = ndimage.gaussian_filter(frac, smooth)
    if not (frac.min() < 0.5 < frac.max()) or np.any(np.asarray(dims) < 2):
        return None, None, 0.0
    verts, faces, normals, _ = measure.marching_cubes(frac, level=0.5,
                                                      spacing=tuple(pitch3))
    vw = verts + lo
    tissue_d = ndimage.binary_dilation(tissue, iterations=tissue_dilate)
    fc = vw[faces].mean(1)
    fn = normals[faces].mean(1); fn /= (np.linalg.norm(fn,axis=1,keepdims=True)+1e-9)
    step = fn * pitch3[None,:] * 1.2
    def in_tissue(P):
        v = vidx(P); return tissue_d[v[:,0], v[:,1], v[:,2]]
    keep = in_tissue(fc + step) & in_tissue(fc - step)
    frac_genuine = float(keep.mean())
    if keep.sum() < 50:
        return None, None, frac_genuine
    m = trimesh.Trimesh(vertices=vw, faces=faces[keep], process=True)
    m.update_faces(m.nondegenerate_faces()); m.remove_unreferenced_vertices()
    comps = m.split(only_watertight=False)
    if len(comps) > 1:
        m = max(comps, key=lambda c: c.area)
    return np.asarray(m.vertices), np.asarray(m.faces), frac_genuine


def band_scores(XYZ, is_domain, Vsheet, Fsheet, band_um=None, pitch_ref=None,
                taubin_iters=40, band_quantile=1/3):
    """Assign each DOMAIN cell a bulge/dent score = normal displacement of the
    open sheet at its nearest vertex; keep the cells CLOSEST to the sheet.

    Band selection is scale-free: keep the `band_quantile` fraction of domain
    cells nearest the interface sheet (default nearest 1/3). If `band_um` is
    given it overrides with an absolute distance threshold.
    Returns (dom_global_idx_in_band, score, band_mask_over_domain)."""
    sheet = trimesh.Trimesh(vertices=Vsheet, faces=Fsheet, process=False)
    sm = sheet.copy(); trimesh.smoothing.filter_taubin(sm, lamb=0.7, nu=-0.72, iterations=taubin_iters)
    Vn = sheet.vertex_normals
    # Vn (marching-cubes normal of the tumor-fraction iso-surface) points toward the
    # DOMAIN side. (Vsheet - smoothed)·Vn > 0 therefore means the sheet is pushed
    # toward the domain here => the domain is locally RECESSED (dent). We negate so
    # that POSITIVE score = domain PROTRUDES into non-domain = convex bulge
    # (the physical surface orientation; positive = convex protrusion).
    disp = -np.einsum('ij,ij->i', (Vsheet - np.asarray(sm.vertices)), Vn)
    dom_idx = np.where(is_domain)[0]
    dom_xyz = XYZ[dom_idx]
    vt = cKDTree(Vsheet)
    dcell, vidx = vt.query(dom_xyz, k=1)
    if band_um is not None:
        inb = dcell < band_um
    else:
        thr = np.quantile(dcell, band_quantile)
        inb = dcell <= thr
    return dom_idx[inb], disp[vidx][inb], inb


def de_and_null(Xn_band, genes, score, xyz_band, section_band, n_perm=200,
                n_top=6, seed=0):
    """Within-stratum bulge-vs-dent DE + azimuthal geometry-breaking null.
    Xn_band: (n_band, G) log-normalized dense/sparse; returns dict."""
    q1, q2 = np.quantile(score, [1/3, 2/3])
    mb = score >= q2; md = score <= q1
    if sp.issparse(Xn_band):
        mean_b = np.asarray(Xn_band[mb].mean(0)).ravel()
        mean_d = np.asarray(Xn_band[md].mean(0)).ravel()
    else:
        mean_b = Xn_band[mb].mean(0); mean_d = Xn_band[md].mean(0)
    lfc = mean_b - mean_d
    order = np.argsort(lfc)
    test_gi = np.concatenate([order[-n_top:], order[:n_top]])   # top up + top down
    Xg = Xn_band[:, test_gi]
    Xg = np.asarray(Xg.todense()) if sp.issparse(Xg) else np.asarray(Xg)
    obs = Xg[mb].mean(0) - Xg[md].mean(0)
    # azimuthal shift within z-bands
    cx, cy = xyz_band[:,0].mean(), xyz_band[:,1].mean()
    az = np.arctan2(xyz_band[:,1]-cy, xyz_band[:,0]-cx)
    z = xyz_band[:,2]
    nb_z = min(8, max(2, len(np.unique(z))))
    zb = np.digitize(z, np.quantile(z, np.linspace(0,1,nb_z+1)[1:-1]))
    rng = np.random.default_rng(seed)
    null = np.zeros((n_perm, len(test_gi)))
    for p in range(n_perm):
        sh = score.copy()
        for b in np.unique(zb):
            mm = np.where(zb==b)[0]
            if len(mm) > 2:
                o = mm[np.argsort(az[mm])]
                sh[o] = np.roll(score[o], int(rng.integers(1, len(o))))
        pb = sh >= q2; pd_ = sh <= q1
        null[p] = Xg[pb].mean(0) - Xg[pd_].mean(0)
    pval = (np.sum(np.abs(null) >= np.abs(obs), axis=0) + 1) / (n_perm + 1)
    return dict(genes=genes[test_gi], obs_lfc=obs, pval=pval,
                n_survive=int((pval < 0.05).sum()), n_tested=len(test_gi),
                lfc_full=lfc, max_abs_lfc=float(np.abs(lfc).max()))
