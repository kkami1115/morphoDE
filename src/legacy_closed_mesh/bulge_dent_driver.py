"""bulge_dent_driver.py — per-dataset bulge-vs-dent DE with stratification + geometry null.
Reuses the interface-bulge-dent-de skill helpers (compute_bulge_dent_scores, tercile,
de_meandiff, decompose_de, azimuthal_shift_null) which must be exec'd into globals first.
"""
import numpy as np, pandas as pd, scipy.sparse as sp, re

def detect_strata(section_ids):
    """Return per-cell stratum label. If section names share a stage/specimen prefix
    (e.g. 'E16-18h_a_S01'), stratify by that prefix; else single stratum 'all'."""
    s = pd.Series([str(x) for x in section_ids])
    # prefix = everything before a trailing _S<digits> or the leading stage token
    pref = s.str.replace(r'_S\d+.*$','',regex=True)         # drop _S01 tail
    pref = pref.str.replace(r'_[a-z]$','',regex=True)       # drop replicate letter _a
    if pref.nunique() > 1 and pref.nunique() < len(s.unique()):
        return pref.values
    # try stage token like E16-18h / L1 / L2 at start
    stage = s.str.extract(r'^([EL][\d\-a-z]+h?|\d+\.\d+)')[0]
    if stage.notna().all() and stage.nunique() > 1 and stage.nunique() < 20:
        return stage.values
    return np.array(['all']*len(s))

def run_interface(adata, spec, skillns, min_band_cells=800, n_perm=200,
                  section_col='section_id', annot_col='annotation'):
    """Compute bulge-vs-dent within-cell-type DE + geometry null for one interface,
    within the largest stratum. Returns (de_df, null_df, meta)."""
    from st3d_interface import compute_interface
    ext = compute_interface(adata, spec)
    g = ext.geom
    mesh = g.mesh
    if mesh is None or len(mesh.vertices) < 50:
        return None, None, {'skip':'no_mesh'}
    cbd = skillns['compute_bulge_dent_scores'](mesh.vertices, mesh.faces)
    face_hull = cbd['hull_dev'][mesh.faces].mean(1)
    face_disp = cbd['smooth_disp'][mesh.faces].mean(1)
    from scipy.stats import spearmanr
    rho = spearmanr(cbd['hull_dev'], cbd['smooth_disp']).correlation
    annot = adata.obs[annot_col].values
    sdf = g.cell_sdf; nf = g.cell_nearest_face
    coords = adata.obsm['spatial_3d'] if 'spatial_3d' in adata.obsm else adata.obsm['spatial_3d_aligned']
    sec = adata.obs[section_col].values if section_col in adata.obs else np.zeros(adata.n_obs)
    # target domain = inside label(s) of the spec  (same cell type => intrinsic)
    inside = set(spec.a_values)
    dom = np.isin(annot, list(inside))
    idx = np.where(dom)[0]
    if len(idx) < min_band_cells:
        return None, None, {'skip':f'too_few_domain_cells={len(idx)}'}
    strata = detect_strata(sec[idx])
    # largest stratum
    vals,cnts = np.unique(strata, return_counts=True)
    big = vals[np.argmax(cnts)]
    m = strata==big
    idxs = idx[m]
    if len(idxs) < min_band_cells:
        return None, None, {'skip':f'largest_stratum_too_small={len(idxs)}'}
    cell_hull = face_hull[nf[idxs]]; cell_disp = face_disp[nf[idxs]]
    xyz = coords[idxs]
    # expression, 5% detection
    X = adata[idxs].X; X = X.tocsr() if sp.issparse(X) else sp.csr_matrix(X)
    det = np.asarray((X>0).mean(0)).ravel(); keep = det>=0.05
    genes = adata.var_names.values[keep]; Xk = X[:,keep].tocsc()
    cells = pd.DataFrame({'x':xyz[:,0],'y':xyz[:,1],'z':xyz[:,2],
                          'section':sec[idxs],'hull_dev':cell_hull,'smooth_disp':cell_disp})
    cells['label_hull'] = skillns['tercile'](cell_hull)
    mb = cells['label_hull'].values=='bulge'; md = cells['label_hull'].values=='dent'
    if mb.sum()<50 or md.sum()<50:
        return None, None, {'skip':'label_imbalance'}
    de = skillns['de_meandiff'](Xk, mb, md, genes)
    # geometry null on top up/down genes
    tg = list(de.head(6)['gene']) + list(de.tail(6)['gene'])
    nullres = skillns['azimuthal_shift_null'](cells, Xk, 'hull_dev', tg, genes, n_perm=n_perm, seed=0)
    de = de.merge(nullres[['gene','p_geomnull','null_sd']], on='gene', how='left')
    meta = {'n_domain_cells':int(len(idx)),'n_stratum_cells':int(len(idxs)),
            'stratum':str(big),'n_strata':int(len(vals)),'rho_defs':float(rho),
            'n_genes':int(keep.sum()),'max_abs_lfc':float(de.lfc.abs().max()),
            'n_null_sig':int((nullres['p_geomnull']<0.05).sum()),'n_null_tested':int(len(nullres)),
            'mesh_verts':int(len(mesh.vertices)),'watertight':bool(mesh.is_watertight)}
    return de, nullres, meta
