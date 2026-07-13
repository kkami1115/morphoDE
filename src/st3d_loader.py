"""
st3d_loader — Layer 1: 適応ローダ + 3D構成ディスパッチャ

設計契約 (3DCurv_design.md §1.1 SpatialSample3D):
  出力 AnnData に以下を保証する
    .obsm['spatial_3d'] : (N,3) float  統一3D座標
    .obs['section_id']  : カテゴリ      切片ID
    .obs['z']           : float         z座標
    .uns['dataset_id'], ['technology'], ['expr_state'], ['z_unit'], ['load_provenance']

中核思想 (§2.4): 静的クラスラベルに頼らず、実行時に obsm 次元と obs 列を
プローブして 3D 構成方式を決める。
"""
from __future__ import annotations
import os, re, glob, tempfile, subprocess, json
import numpy as np
import anndata as ad

# ---- z/section を示唆する obs 列名パターン（優先度順） ----
Z_NUMERIC_PATTS  = [r'^z$', r'new_z', r'z_coord', r'z_um', r'z_pos', r'zpos']
SECTION_PATTS    = [r'slice', r'section', r'^batch$', r'sample_?id', r'array', r'timepoint',
                    r'stage', r'z_?index', r'layer(?!_)', r'puck', r'fov_?z']
ANNOT_PATTS      = [r'annotation', r'cell[_ ]?type', r'celltype', r'region', r'domain',
                    r'layer_annotation', r'bin_annotation', r'clusters?', r'leiden']

def _match(cols, patts):
    out=[]
    for p in patts:
        for c in cols:
            if re.search(p, str(c), re.I) and c not in out:
                out.append(c)
    return out

# ---------------------------------------------------------------
# Decompressor: zst / zip / tar を作業ディレクトリへ展開して h5ad パス群を返す
# ---------------------------------------------------------------
def decompress_to_h5ad(data_dir: str, workdir: str) -> list[str]:
    """data_dir 配下の .h5ad.zst を workdir に解凍。既存 .h5ad はそのまま。
    zip/tar は展開して内部 h5ad を探す。CSVマトリクス系は別ローダ (未対応時は空)。"""
    os.makedirs(workdir, exist_ok=True)
    h5ads=[]
    zsts = glob.glob(os.path.join(data_dir,'**','*.h5ad.zst'), recursive=True)
    for z in zsts:
        out=os.path.join(workdir, os.path.basename(z)[:-4])  # strip .zst
        if not os.path.exists(out):
            subprocess.run(['zstd','-d','-q','-f','-o',out,z], check=True)
        h5ads.append(out)
    # plain h5ad already present
    for h in glob.glob(os.path.join(data_dir,'**','*.h5ad'), recursive=True):
        if not h.endswith('.zst') and h not in h5ads:
            h5ads.append(h)
    # .h5ad.gz -> gunzip into workdir
    import gzip, shutil
    for g in glob.glob(os.path.join(data_dir,'**','*.h5ad.gz'), recursive=True):
        out=os.path.join(workdir, os.path.basename(g)[:-3])  # strip .gz
        if not os.path.exists(out):
            with gzip.open(g,'rb') as fi, open(out,'wb') as fo:
                shutil.copyfileobj(fi, fo)
        if out not in h5ads: h5ads.append(out)
    # zip / tar: extract then recurse for h5ad
    if not h5ads:
        for arch in glob.glob(os.path.join(data_dir,'**','*.zip'), recursive=True):
            subprocess.run(['unzip','-n','-q',arch,'-d',workdir], check=False)
        for arch in (glob.glob(os.path.join(data_dir,'**','*.tar'), recursive=True)
                     + glob.glob(os.path.join(data_dir,'**','*.tar.gz'), recursive=True)):
            subprocess.run(['tar','xf',arch,'-C',workdir], check=False)
        h5ads = glob.glob(os.path.join(workdir,'**','*.h5ad'), recursive=True)
    return sorted(h5ads)

# ---------------------------------------------------------------
# Prober: 1つの AnnData がどの z-供給源を持つか判定
# ---------------------------------------------------------------
def probe_anndata(a: ad.AnnData) -> dict:
    """returns dict(zmode, coord_key, z_col, section_col, annot_col, spatial_dim)"""
    info=dict(zmode=None, coord_key=None, z_col=None, section_col=None,
              annot_col=None, spatial_dim=None)
    # 1) 3D coords in obsm ?
    coord3=None
    for k,v in a.obsm.items():
        arr=np.asarray(v)
        if arr.ndim==2 and arr.shape[1]>=3 and re.search(r'spat|coord|3d|umap', k, re.I) is None:
            coord3=k
        if arr.ndim==2 and arr.shape[1]>=3 and re.search(r'spat|coord|3d', k, re.I):
            coord3=k; break
    # explicit spatial 2D key
    spatial2=None
    for k in ('spatial','X_spatial','spatial_3d'):
        if k in a.obsm:
            arr=np.asarray(a.obsm[k]); 
            if arr.shape[1]>=3: coord3=coord3 or k
            elif arr.shape[1]==2: spatial2=spatial2 or k
    cols=list(a.obs.columns)
    znum=_match(cols, Z_NUMERIC_PATTS)
    sect=_match(cols, SECTION_PATTS)
    ann =_match(cols, ANNOT_PATTS)
    info['annot_col']=ann[0] if ann else None
    info['section_col']=sect[0] if sect else None
    info['z_col']=znum[0] if znum else None
    # coordinates stored as obs columns (no obsm) — synthesize obsm['spatial']
    if not coord3 and not spatial2:
        def find_col(cands):
            for c in cands:
                for col in cols:
                    if col.lower()==c: return col
            return None
        xc=find_col(['spatial_x','x','imagecol','new_x','x_coord','coord_x','array_col','pxl_col_in_fullres'])
        yc=find_col(['spatial_y','y','imagerow','new_y','y_coord','coord_y','array_row','pxl_row_in_fullres'])
        zc=find_col(['spatial_z','z','new_z','z_coord','coord_z'])
        if xc and yc:
            if zc:
                a.obsm['spatial']=a.obs[[xc,yc,zc]].astype(float).values; coord3='spatial'
                info['spatial_dim']=3
            else:
                a.obsm['spatial']=a.obs[[xc,yc]].astype(float).values; spatial2='spatial'
                info['spatial_dim']=2
            info['obs_coord_cols']=[xc,yc]+([zc] if zc else [])
    if coord3:
        info['zmode']='native_3d'; info['coord_key']=coord3; info['spatial_dim']=info['spatial_dim'] or 3
    elif spatial2 and (znum or sect):
        # 2D coords + explicit z/section column in same file
        info['zmode']='coord2d_plus_zcol'; info['coord_key']=spatial2; info['spatial_dim']=2
    elif spatial2:
        # single 2D section, z must come from filename / stacking across files
        info['zmode']='section_2d_only'; info['coord_key']=spatial2; info['spatial_dim']=2
    else:
        info['zmode']='no_spatial'
    return info

# ---------------------------------------------------------------
# Assembler3D: プローブ結果に応じて spatial_3d を構成し契約を満たす
# ---------------------------------------------------------------
def _section_z_from_name(fname: str) -> float | None:
    """ファイル名から切片インデックスを推定（S01, _z10, slice3 等）"""
    for p in [r'[_.]S(\d+)', r'[_.]z(\d+)', r'slice[_-]?(\d+)', r'sec[_-]?(\d+)']:
        m=re.search(p, os.path.basename(fname), re.I)
        if m: return float(m.group(1))
    return None

def assemble_3d(h5ads: list[str], dataset_id: str, technology: str='',
                z_spacing: float=1.0, z_unit: str='section_index',
                max_files: int|None=None) -> ad.AnnData:
    """h5ad 群を SpatialSample3D 契約の単一 AnnData へ。
    - native_3d: そのまま座標流用（複数ファイルは concat）
    - coord2d_plus_zcol: 2D + z/section 列で3D化
    - section_2d_only: 1ファイル1切片 → ファイル順で z 割当し積層
    """
    if max_files: h5ads=h5ads[:max_files]
    parts=[]; prov=[]
    for i,h in enumerate(h5ads):
        if isinstance(h, ad.AnnData):
            a=h; hname=a.uns.get('_srcname', f'part{i}')
        else:
            a=ad.read_h5ad(h); hname=os.path.basename(h)
        p=probe_anndata(a)
        prov.append(dict(file=hname, **{k:p[k] for k in ('zmode','coord_key','z_col','section_col')}))
        xy=None; z=None; sect=None
        if p['zmode']=='native_3d':
            c=np.asarray(a.obsm[p['coord_key']]).astype(float)
            xy=c[:,:2]; z=c[:,2]
            sect=a.obs[p['section_col']].astype(str).values if p['section_col'] else \
                 np.round(z/ (np.median(np.diff(np.unique(z))) or 1)).astype(int).astype(str)
        elif p['zmode']=='coord2d_plus_zcol':
            xy=np.asarray(a.obsm[p['coord_key']]).astype(float)[:,:2]
            if p['z_col'] is not None:
                z=a.obs[p['z_col']].astype(float).values
                sect=(a.obs[p['section_col']].astype(str).values if p['section_col']
                      else np.round(z).astype(int).astype(str))
            else:  # section col only -> map categories to integer z * spacing
                s=a.obs[p['section_col']].astype(str)
                nsec=s.nunique()
                if nsec<=1 and len(h5ads)>1:
                    # part は単一切片。複数part間ではファイル通し番号 i で z を振る
                    z=np.full(a.n_obs, float(i)*z_spacing); sect=s.values
                else:
                    codes=s.astype('category').cat.codes.values
                    z=codes.astype(float)*z_spacing; sect=s.values
        elif p['zmode']=='section_2d_only':
            xy=np.asarray(a.obsm[p['coord_key']]).astype(float)[:,:2]
            zi=_section_z_from_name(hname); zi = zi if zi is not None else float(i)
            z=np.full(a.n_obs, zi*z_spacing); sect=np.full(a.n_obs, f"sec{int(zi):03d}")
        else:
            a.file if False else None
            raise ValueError(f"{h}: no usable spatial coords (zmode={p['zmode']})")
        a.obsm['spatial_3d']=np.column_stack([xy, z]).astype(float)
        a.obs['section_id']=sect
        a.obs['z']=z.astype(float)
        # keep annotation under canonical name if found
        if p['annot_col'] and p['annot_col']!='annotation':
            a.obs['annotation']=a.obs[p['annot_col']].astype(str).values
        parts.append(a)
    # concat (outer join on genes to be safe)
    if len(parts)==1:
        out=parts[0]
    else:
        keys=[(h.uns.get('_srcname', f'part{i}') if isinstance(h, ad.AnnData) else os.path.basename(h))
              for i,h in enumerate(h5ads)]
        out=ad.concat(parts, join='outer', label='_srcfile', keys=keys, index_unique='-')
        # outer join で「その切片に無い遺伝子」は NaN 埋めされる → count/expr データでは
        # 未測定=0 として扱うのが正しい（NaN が下流の正規化/DEを壊す）。
        import scipy.sparse as _sp
        if _sp.issparse(out.X):
            out.X.data=np.nan_to_num(out.X.data, nan=0.0)
        else:
            out.X=np.nan_to_num(out.X, nan=0.0)
    out.uns['dataset_id']=dataset_id
    out.uns['technology']=technology
    out.uns['z_unit']=z_unit
    out.uns['expr_state']='unknown'  # caller sets after inspecting .X
    out.uns['load_provenance']=json.dumps(prov)[:20000]
    return out

def read_matrix_plus_coords(counts_path: str, coords_path: str, meta_path: str=None,
                            genes_are_rows: bool=True, dataset_id: str='') -> ad.AnnData:
    """汎用リーダ: (遺伝子×細胞 or 細胞×遺伝子) の csv/tsv(.gz) 行列 + 座標CSV(+メタ) → AnnData。
    座標CSVは index=cell、列に row/col or x/y or imagerow/imagecol、section 列を含む想定。
    counts の向きは genes_are_rows で指定（Seurat輸出は遺伝子が行）。"""
    import pandas as pd
    # sniff delimiter from the first line (R write.table often uses space + quotes;
    # GEO tsv uses tab; csv uses comma). engine='python' handles regex whitespace.
    import gzip as _gz
    opn=_gz.open if counts_path.endswith('.gz') else open
    with opn(counts_path,'rt') as f: head=f.readline()
    if '\t' in head: sep='\t'
    elif ',' in head: sep=','
    else: sep=r'\s+'
    cnt=pd.read_csv(counts_path, sep=sep, index_col=0, quotechar='"',
                    engine='python' if sep==r'\s+' else 'c')
    cnt.index=cnt.index.astype(str).str.strip('"')
    cnt.columns=[str(c).strip('"') for c in cnt.columns]
    if genes_are_rows: cnt=cnt.T          # -> cells x genes
    if coords_path is None:
        # coordinates embedded in spot names like "10x13" (col x row) — common in ST arrays
        import re as _re
        m=[_re.match(r'^(\d+\.?\d*)x(\d+\.?\d*)$', str(s)) for s in cnt.index]
        if all(m):
            a=ad.AnnData(X=cnt.values.astype('float32'))
            a.obs_names=list(cnt.index); a.var_names=list(cnt.columns)
            a.obsm['spatial']=np.array([[float(g.group(1)),float(g.group(2))] for g in m])
            a.uns['_srcname']=os.path.basename(counts_path)
            return a
        raise ValueError("coords_path=None but spot names are not 'NxM' coordinates")
    coo=pd.read_csv(coords_path, index_col=0)
    # align cells present in both
    cells=[c for c in cnt.index if c in set(coo.index)]
    cnt=cnt.loc[cells]; coo=coo.loc[cells]
    a=ad.AnnData(X=cnt.values.astype('float32'))
    a.obs_names=list(cells); a.var_names=list(cnt.columns)
    # pick xy columns
    def pick(cols, cands):
        for c in cands:
            for col in cols:
                if col.lower()==c: return col
        return None
    xcol=pick(coo.columns, ['imagecol','x','col','new_x','pxl_col_in_fullres'])
    ycol=pick(coo.columns, ['imagerow','y','row','new_y','pxl_row_in_fullres'])
    seccol=pick(coo.columns, ['section','slice','z','sample','array'])
    if xcol and ycol:
        a.obsm['spatial']=coo[[xcol,ycol]].values.astype(float)
    if seccol: a.obs['section']=coo[seccol].astype(str).values
    if meta_path:
        meta=pd.read_csv(meta_path, index_col=0).reindex(cells)
        for c in meta.columns:
            a.obs[c]=meta[c].values
    a.uns['_srcname']=os.path.basename(counts_path)
    return a


def read_h5ad_with_meta_csv(h5ad_path: str, meta_csv: str, cell_col: str=None,
                            xyz_cols=('x','y','z'), annot_col: str=None,
                            max_cells: int=None) -> ad.AnnData:
    """発現 h5ad と、座標(+annotation)を持つ別 metadata CSV を cell id で結合。
    Allen ABC MERFISH のように obsm 座標を持たず cell_metadata.csv に x,y,z がある形式向け。
    cell_col 未指定なら CSV の最初の列を cell id とみなし h5ad.obs_names に合わせる。"""
    import pandas as pd
    a=ad.read_h5ad(h5ad_path)
    if max_cells and a.n_obs>max_cells:
        a=a[:max_cells].copy()
    meta=pd.read_csv(meta_csv, index_col=(cell_col or 0), low_memory=False)
    meta=meta.reindex(a.obs_names)
    xc,yc,zc=xyz_cols
    # detect a discrete section label; if present, z should come from it (serial sections),
    # not from a per-cell continuous z. Keep xy 2D + section column.
    seccol=None
    for sc in ('brain_section_label','section','slice','z_section'):
        if sc in meta.columns: seccol=sc; break
    if seccol:
        a.obs['section']=meta[seccol].astype(str).values
        if xc in meta.columns and yc in meta.columns:
            a.obsm['spatial']=meta[[xc,yc]].astype(float).values   # 2D + section drives z
    else:
        have=[c for c in (xc,yc,zc) if c in meta.columns]
        if len(have)>=2:
            cols=[xc,yc]+([zc] if zc in meta.columns else [])
            a.obsm['spatial']=meta[cols].astype(float).values
    if annot_col and annot_col in meta.columns:
        a.obs['annotation']=meta[annot_col].astype(str).values
    a.uns['_srcname']=os.path.basename(h5ad_path)
    return a


def read_10x_h5(h5_path: str, positions_csv: str=None) -> ad.AnnData:
    """10x Visium .h5 (filtered_feature_bc_matrix.h5) を読む。
    positions_csv があれば tissue_positions を obsm['spatial'] に付与。"""
    import scanpy as sc
    import pandas as pd
    a=sc.read_10x_h5(h5_path); a.var_names_make_unique()
    if positions_csv and os.path.exists(positions_csv):
        pos=pd.read_csv(positions_csv, header=None, index_col=0)
        pos=pos.reindex(a.obs_names)
        # Visium tissue_positions: cols [in_tissue,row,col,imagerow,imagecol]
        if pos.shape[1]>=5:
            a.obsm['spatial']=pos.iloc[:,[4,3]].values.astype(float)
            a.obs['section']=os.path.basename(os.path.dirname(h5_path))
    a.uns['_srcname']=os.path.basename(h5_path)
    return a


def read_mtx_triplet(dir_or_prefix: str) -> ad.AnnData:
    """10x mtx トリプレット (matrix.mtx.gz + features/genes + barcodes) を読む。"""
    import scanpy as sc
    d=dir_or_prefix if os.path.isdir(dir_or_prefix) else os.path.dirname(dir_or_prefix)
    a=sc.read_mtx(glob.glob(os.path.join(d,'*matrix.mtx*'))[0]).T
    import pandas as pd
    feat=glob.glob(os.path.join(d,'*features*'))+glob.glob(os.path.join(d,'*genes*'))
    bc=glob.glob(os.path.join(d,'*barcodes*'))
    if feat: a.var_names=pd.read_csv(feat[0],sep='\t',header=None)[ (1 if 'features' in feat[0] else 0) ].values[:a.n_vars]
    if bc: a.obs_names=pd.read_csv(bc[0],sep='\t',header=None)[0].values[:a.n_obs]
    a.var_names_make_unique(); a.uns['_srcname']=os.path.basename(d)
    return a


def register_sections_rigid(sample: ad.AnnData) -> ad.AnnData:
    """Class B フォールバック整列: 各切片を重心中心化 + PCA主軸整列（剛体）。
    spatial_3d の xy を上書きし、uns に registration=rigid_centroid_pca を記録。
    連続2D切片の粗整列。より高精度が要るときは PASTE/STAligner に差し替える。"""
    xyz=sample.obsm['spatial_3d'].astype(float).copy()
    sec=sample.obs['section_id'].astype(str).values
    ref_y=ref_x=None
    for z in sorted(set(sec)):
        sel=sec==z; pts=xyz[sel,:2]; pts=pts-pts.mean(0)
        _,_,vt=np.linalg.svd(pts-pts.mean(0), full_matrices=False)
        ang=np.arctan2(vt[0][1],vt[0][0])
        R=np.array([[np.cos(-ang),-np.sin(-ang)],[np.sin(-ang),np.cos(-ang)]])
        pts=pts@R.T
        if ref_y is not None:
            if np.sign(np.mean(pts[:,1]**3))!=ref_y: pts[:,1]*=-1
            if np.sign(np.mean(pts[:,0]**3))!=ref_x: pts[:,0]*=-1
        else:
            ref_y=np.sign(np.mean(pts[:,1]**3)); ref_x=np.sign(np.mean(pts[:,0]**3))
        xyz[sel,:2]=pts
    sample.obsm['spatial_3d']=xyz
    sample.uns['registration']='rigid_centroid_pca'
    return sample

def needs_registration(sample: ad.AnnData, thresh: float=None) -> bool:
    """切片重心のばらつきが大きければ registration 要と判定。
    thresh 未指定なら切片内スケール(中央値)の10%を基準にする。"""
    xyz=sample.obsm['spatial_3d']; sec=sample.obs['section_id'].astype(str).values
    cents=np.array([xyz[sec==z,:2].mean(0) for z in sorted(set(sec))])
    if len(cents)<2: return False
    spread=np.hypot(cents[:,0].std(), cents[:,1].std())
    if thresh is None:
        scale=np.median([np.ptp(xyz[sec==z,0]) for z in sorted(set(sec))])
        thresh=0.1*scale
    return spread>thresh

def load(dataset_id: str, data_root: str='/mnt/f/3DCurv/data', workdir: str|None=None,
         technology: str='', z_spacing: float=1.0, max_files: int|None=None,
         auto_register: bool=True, parts: list=None) -> ad.AnnData:
    """トップレベル: dataset_id -> SpatialSample3D 契約の AnnData
    auto_register=True で、連続2D切片が整列を要する場合に剛体整列を自動適用。
    parts: h5ad以外の形式は read_matrix_plus_coords / read_10x_h5 / read_mtx_triplet で
           AnnData を作り、そのリストを parts= に渡せば h5ad 展開を飛ばして組み上げる。"""
    data_dir=os.path.join(data_root, dataset_id)
    # 展開先は /mnt/f を汚さないため既定で /mnt/d/3DCurv_tmp/_work を使う
    # (/mnt/d が無ければ従来の data_root/_st3d_work にフォールバック)
    if workdir is None:
        _tmp="/mnt/d/3DCurv_tmp/_work"
        _base=_tmp if os.path.isdir("/mnt/d/3DCurv_tmp") else os.path.join(data_root,'_st3d_work')
        workdir=os.path.join(_base, dataset_id)
    if parts is not None:
        items=parts
    else:
        items=decompress_to_h5ad(data_dir, workdir)
        if not items:
            raise FileNotFoundError(f"{dataset_id}: no h5ad found under {data_dir} "
                                    f"(non-h5ad formats: read via read_* helpers and pass parts=)")
    s=assemble_3d(items, dataset_id, technology=technology, z_spacing=z_spacing, max_files=max_files)
    s.uns['registration']='none'
    if auto_register and needs_registration(s):
        s=register_sections_rigid(s)
    # validate contract
    assert 'spatial_3d' in s.obsm and s.obsm['spatial_3d'].shape[1]==3
    assert 'section_id' in s.obs and 'z' in s.obs
    return s
