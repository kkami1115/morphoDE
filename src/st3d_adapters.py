"""
st3d_adapters — データセット専用アダプタのレジストリ

3DCurv の各データセットは寄託形式がまちまち（GEO RAW tsv / 10x Visium outs / Slide-seq puck /
DBiT / 別CSV座標 / 入れ子tar …）。汎用リーダで推測せず、**1データセット=1アダプタ**で
確実に読む。各アダプタは `-> list[AnnData]`（切片ごとの parts）を返し、
`st3d_loader.load(dataset_id, parts=adapter())` に渡すと契約 SpatialSample3D になる。

各 AnnData part の最低要件:
  .obsm['spatial'] (N,2 or 3)  +  .obs['section']（切片ID; 単一切片なら定数でも可）
  可能なら .obs['annotation']（細胞型/領域） も付与。
"""
from __future__ import annotations
import os, glob, gzip, tarfile, zipfile, tempfile, subprocess, shutil, re
import numpy as np, pandas as pd, anndata as ad, scanpy as sc

ROOT_DEFAULT="/mnt/f/3DCurv/data"
REGISTRY={}   # dataset_id -> adapter function

def adapter(dataset_id):
    def deco(fn):
        REGISTRY[dataset_id]=fn
        fn._dataset_id=dataset_id
        return fn
    return deco

def get(dataset_id):
    if dataset_id not in REGISTRY:
        raise KeyError(f"no dedicated adapter for '{dataset_id}'. registered: {sorted(REGISTRY)}")
    return REGISTRY[dataset_id]

# ---------- 共通ヘルパ ----------
def _gunzip_to(src, dstdir):
    os.makedirs(dstdir, exist_ok=True)
    out=os.path.join(dstdir, os.path.basename(src)[:-3] if src.endswith('.gz') else os.path.basename(src))
    if not os.path.exists(out):
        with gzip.open(src,'rb') as fi, open(out,'wb') as fo: shutil.copyfileobj(fi,fo)
    return out

def _read_csv_any(path, **kw):
    """csv/tsv(.gz) を区切り自動判定で読む。"""
    opn=gzip.open if path.endswith('.gz') else open
    with opn(path,'rt') as f: head=f.readline()
    sep='\t' if '\t' in head else (',' if ',' in head else r'\s+')
    eng='python' if sep==r'\s+' else 'c'
    return pd.read_csv(path, sep=sep, engine=eng, **kw)

def _extract_members(archive, dstdir, pattern):
    """tar/zip から pattern(正規表現) にマッチするメンバだけ dstdir へ展開し、パス群を返す。"""
    os.makedirs(dstdir, exist_ok=True); out=[]
    rx=re.compile(pattern)
    def _junk(name):
        base=os.path.basename(name)
        return '__MACOSX' in name or base.startswith('._') or base=='.DS_Store'
    if archive.endswith(('.tar','.tar.gz','.tgz')):
        with tarfile.open(archive) as tf:
            for m in tf.getmembers():
                if m.isfile() and rx.search(m.name) and not _junk(m.name):
                    tf.extract(m, dstdir); out.append(os.path.join(dstdir,m.name))
    else:
        with zipfile.ZipFile(archive) as zf:
            for n in zf.namelist():
                if not n.endswith('/') and rx.search(n) and not _junk(n):
                    zf.extract(n, dstdir); out.append(os.path.join(dstdir,n))
    return sorted(out)


def _spotname_xy(index):
    """spot名 'AxB'（'41x50'）→ (N,2) 座標。全一致しなければ None。"""
    m=[re.match(r'^(\d+\.?\d*)x(\d+\.?\d*)$', str(s)) for s in index]
    if all(m):
        return np.array([[float(g.group(1)), float(g.group(2))] for g in m])
    return None


# ==================================================================
# GEO tsv / Slide-seq puck ファミリ
# ==================================================================
@adapter("slideseq2_kidney")
def _slideseq2_kidney(root=ROOT_DEFAULT, work=None, max_sections=None):
    """Slide-seqV2 腎: puck 毎に BeadLocations(barcode,x,y) + MappedDGE(gene×barcode)。
    各 puck = 1 切片。"""
    d=os.path.join(root,"slideseq2_kidney")
    tar=glob.glob(d+"/*.tar")[0]
    work=work or os.path.join(root,"_adapt_work","slideseq2_kidney")
    beads=_extract_members(tar, work, r'BeadLocationsForR\.csv\.gz$')
    beads=[b for b in beads if '_qc' not in b]
    dges =_extract_members(tar, work, r'MappedDGEForR\.csv\.gz$')
    def puck_id(p): return re.search(r'(Puck_\d+_\d+)', os.path.basename(p)).group(1)
    dge_by={puck_id(p):p for p in dges}
    parts=[]
    for b in (beads[:max_sections] if max_sections else beads):
        pid=puck_id(b)
        if pid not in dge_by: continue
        loc=pd.read_csv(b)
        loc.columns=[str(c).strip().strip('"').strip("'") for c in loc.columns]
        bc=[c for c in loc.columns if c.lower().startswith('barcode')]
        xc=[c for c in loc.columns if c.lower() in ('xcoord','x')]
        yc=[c for c in loc.columns if c.lower() in ('ycoord','y')]
        if not (bc and xc and yc):   # 異形式の bead ファイルはスキップ
            continue
        loc=loc.set_index(bc[0])
        dge=pd.read_csv(dge_by[pid], index_col=0)      # gene x barcode
        dge=dge.T                                       # barcode x gene
        common=[c for c in dge.index if c in set(loc.index)]
        dge=dge.loc[common]; loc=loc.loc[common]
        a=ad.AnnData(X=dge.values.astype('float32'))
        a.obs_names=[f"{pid}_{c}" for c in common]; a.var_names=list(dge.columns)
        a.obsm['spatial']=loc[[xc[0],yc[0]]].values.astype(float)
        a.obs['section']=pid; a.uns['_srcname']=pid
        parts.append(a)
    return parts


@adapter("dbit_seq_mouse_embryo")
def _dbit(root=ROOT_DEFAULT, work=None, max_sections=None):
    """DBiT-seq 胚: GSM 毎に spot×gene tsv、spot名 'AxB' に座標埋め込み。GSM=切片。"""
    d=os.path.join(root,"dbit_seq_mouse_embryo")
    tar=glob.glob(d+"/*.tar")[0]
    work=work or os.path.join(root,"_adapt_work","dbit")
    tsvs=_extract_members(tar, work, r'GSM\d+_[^/]*\.tsv\.gz$')
    parts=[]
    for t in (tsvs[:max_sections] if max_sections else tsvs):
        sec=re.search(r'(GSM\d+)', os.path.basename(t)).group(1)
        m=pd.read_csv(t, sep='\t', index_col=0)
        xy=_spotname_xy(m.index)
        if xy is None: continue
        a=ad.AnnData(X=m.values.astype('float32'))
        a.obs_names=[f"{sec}_{s}" for s in m.index]; a.var_names=list(m.columns)
        a.obsm['spatial']=xy; a.obs['section']=sec; a.uns['_srcname']=sec
        parts.append(a)
    return parts


@adapter("scc_ji_serial")
def _scc(root=ROOT_DEFAULT, work=None, max_sections=None):
    """cSCC 連続切片: GSM 毎に stdata.tsv(spot×gene, spot名 'AxB') + spot selection(精密座標)。"""
    d=os.path.join(root,"scc_ji_serial")
    tar=glob.glob(d+"/*.tar")[0]
    work=work or os.path.join(root,"_adapt_work","scc")
    std=_extract_members(tar, work, r'_stdata\.tsv\.gz$')
    sel=_extract_members(tar, work, r'spot_data-selection-.*\.tsv\.gz$')
    def key(p):
        mm=re.search(r'(GSM\d+)', os.path.basename(p)); return mm.group(1) if mm else p
    sel_by={key(s):s for s in sel}
    parts=[]
    for t in (std[:max_sections] if max_sections else std):
        sec=key(t)
        m=pd.read_csv(t, sep='\t', index_col=0)   # spot 'AxB' x gene
        # coords: prefer precise selection file (has x,y matching 'AxB')
        xy=None
        if sec in sel_by:
            sd=pd.read_csv(sel_by[sec], sep='\t')
            sd.index=[f"{int(round(r.x))}x{int(round(r.y))}" for r in sd.itertuples()]
            common=[s for s in m.index if s in set(sd.index)]
            if common:
                m=m.loc[common]; xy=sd.loc[common][['new_x','new_y']].values.astype(float)
        if xy is None:
            xy=_spotname_xy(m.index)
            if xy is None: continue
        a=ad.AnnData(X=m.values.astype('float32'))
        a.obs_names=[f"{sec}_{s}" for s in m.index]; a.var_names=list(m.columns)
        a.obsm['spatial']=xy; a.obs['section']=sec; a.uns['_srcname']=sec
        parts.append(a)
    return parts


@adapter("asp_human_heart")
def _asp_heart(root=ROOT_DEFAULT, work=None, max_sections=None):
    """発生心臓 ST: exprMat 単一tsv、列名 'W x X x Y'（週×x×y）。週=切片(z)。"""
    d=os.path.join(root,"asp_human_heart")
    zp=os.path.join(d,"Developmental_heart_ST_raw.tsv.zip")
    work=work or os.path.join(root,"_adapt_work","asp_heart")
    tsv=_extract_members(zp, work, r'exprMat_all_weeks_raw\.tsv$')[0]
    try:
        m=pd.read_csv(tsv, sep='\t', index_col=0)
    except UnicodeDecodeError:
        m=pd.read_csv(tsv, sep='\t', index_col=0, encoding='latin-1')   # gene x spot
    m=m.T                                        # spot x gene
    parsed=[re.match(r'^(\d+)x(\d+\.?\d*)x(\d+\.?\d*)$', str(s)) for s in m.index]
    keep=[i for i,g in enumerate(parsed) if g]
    m=m.iloc[keep]; parsed=[parsed[i] for i in keep]
    week=np.array([int(g.group(1)) for g in parsed])
    xy=np.array([[float(g.group(2)), float(g.group(3))] for g in parsed])
    parts=[]
    weeks=sorted(set(week))
    for w in (weeks[:max_sections] if max_sections else weeks):
        sel=week==w
        a=ad.AnnData(X=m.values[sel].astype('float32'))
        a.obs_names=[str(s) for s in m.index[sel]]; a.var_names=list(m.columns)
        a.obsm['spatial']=xy[sel]; a.obs['section']=f"W{w}"; a.uns['_srcname']=f"W{w}"
        parts.append(a)
    return parts


def _iriseq(root=ROOT_DEFAULT, work=None, max_sections=None):
    """IRISeq 老化脳 — 要調査（未対応）。
    tar 内は GSM 毎の `*_connection.spatial.csv.gz` のみで、実測 3 列
    （barcode + 2 本のシーケンス文字列）＝IRISeq のリード接続情報であり、
    **座標も発現行列も含まない**（当初 docstring は誤りだった）。
    発現/座標は demo.mtx / demo.spatial.txt 系（全体で1セット）に別形式で入っており、
    barcode を介した連続切片25枚への割り当て規則が非自明。専用の再構築が必要なため
    レジストリには登録せず、要調査扱いとする。"""
    raise NotImplementedError("iriseq_aging_brain: non-standard IRISeq layout (connection.spatial "
                              "has no coords/expression); needs dedicated reconstruction — see docstring.")


# ==================================================================
# 10x Visium ファミリ
# ==================================================================
def _read_tissue_positions(path):
    """Visium tissue_positions(_list).csv(.gz) -> DataFrame index=barcode, cols with imagerow/imagecol.
    ヘッダ有無の両対応。返り値の 'X','Y' は fullres pixel 座標。"""
    opn=gzip.open if path.endswith('.gz') else open
    with opn(path,'rt') as f: first=f.readline()
    has_header='barcode' in first.lower()
    if has_header:
        df=pd.read_csv(path)
        df.columns=[c.strip() for c in df.columns]
        bc=[c for c in df.columns if 'barcode' in c.lower()][0]
        df=df.set_index(bc)
        xcol=[c for c in df.columns if c.lower() in ('pxl_col_in_fullres','imagecol','x')][0]
        ycol=[c for c in df.columns if c.lower() in ('pxl_row_in_fullres','imagerow','y')][0]
    else:
        df=pd.read_csv(path, header=None, index_col=0)
        # legacy: [in_tissue,row,col,imagerow,imagecol]
        df.columns=['in_tissue','row','col','imagerow','imagecol'][:df.shape[1]]
        xcol,ycol='imagecol','imagerow'
    out=pd.DataFrame({'X':df[xcol].astype(float),'Y':df[ycol].astype(float)}, index=df.index)
    return out


def _visium_part(h5=None, mtx_dir=None, positions=None, section=None):
    """1 Visium 切片を AnnData に。h5 か mtx_dir のいずれかを与える。"""
    import scanpy as sc
    if h5:
        a=sc.read_10x_h5(h5)
    else:
        a=sc.read_mtx(glob.glob(os.path.join(mtx_dir,'*matrix.mtx*'))[0]).T
        feat=(glob.glob(os.path.join(mtx_dir,'*features*'))+glob.glob(os.path.join(mtx_dir,'*genes*')))[0]
        bc=glob.glob(os.path.join(mtx_dir,'*barcodes*'))[0]
        fdf=pd.read_csv(feat,sep='\t',header=None)
        a.var_names=fdf[1].values[:a.n_vars] if fdf.shape[1]>1 else fdf[0].values[:a.n_vars]
        a.obs_names=pd.read_csv(bc,sep='\t',header=None)[0].values[:a.n_obs]
    a.var_names_make_unique()
    if positions:
        pos=_read_tissue_positions(positions).reindex(a.obs_names)
        a.obsm['spatial']=pos[['X','Y']].values.astype(float)
    a.obs['section']=section; a.uns['_srcname']=section
    return a


@adapter("colon_cancer_visiumhd")
def _colon(root=ROOT_DEFAULT, work=None, max_sections=None):
    """大腸がん Visium HD: tar 内に GSM 毎 filtered .h5 + 別 tissue_positions.csv.gz。GSM=切片。"""
    d=os.path.join(root,"colon_cancer_visiumhd")
    tar=[f for f in glob.glob(d+"/**/*",recursive=True) if f.endswith(('.tar','.tar.gz'))][0]
    work=work or os.path.join(root,"_adapt_work","colon")
    h5s=_extract_members(tar, work, r'filtered_feature_bc_matrix\.h5$')
    poss=_extract_members(tar, work, r'tissue_positions.*\.csv\.gz$')
    # match by GSM sample token
    def tok(p):
        mm=re.search(r'(P\d+CRC|P\d+NAT|BC\d+)', os.path.basename(p)); return mm.group(1) if mm else None
    pos_by={}
    for p in poss:
        t=tok(p);
        if t: pos_by.setdefault(t,p)
    parts=[]
    for h in (h5s[:max_sections] if max_sections else h5s):
        sec=re.search(r'(GSM\d+[^/]*?)_count', os.path.basename(h))
        sec=sec.group(1) if sec else os.path.basename(h)[:20]
        t=tok(h); pos=pos_by.get(t)
        try:
            a=_visium_part(h5=h, positions=pos, section=sec)
            if 'spatial' in a.obsm and not np.isnan(a.obsm['spatial']).all():
                parts.append(a)
        except Exception: pass
    return parts


@adapter("human_embryonic_limb_visium")
def _limb(root=ROOT_DEFAULT, work=None, max_sections=None):
    """胚肢 Visium: tar.gz 毎に1サンプル（filtered_feature_bc_matrix/ + spatial/positions）。"""
    d=os.path.join(root,"human_embryonic_limb_visium")
    tars=[f for f in glob.glob(d+"/**/*",recursive=True) if f.endswith(('.tar','.tar.gz'))]
    work=work or os.path.join(root,"_adapt_work","limb")
    parts=[]
    for t in (tars[:max_sections] if max_sections else tars):
        sec=re.search(r'(WSSS[^/.]+)', os.path.basename(t))
        sec=sec.group(1) if sec else os.path.basename(t).split('.')[0]
        sw=os.path.join(work, sec)
        _extract_members(t, sw, r'filtered_feature_bc_matrix/.*(matrix\.mtx|barcodes|features)')
        _extract_members(t, sw, r'tissue_positions_list\.csv$')
        mtxdir=glob.glob(sw+"/**/filtered_feature_bc_matrix",recursive=True)
        pos=glob.glob(sw+"/**/tissue_positions_list.csv",recursive=True)
        if not mtxdir or not pos: continue
        try:
            parts.append(_visium_part(mtx_dir=mtxdir[0], positions=pos[0], section=sec))
        except Exception: pass
    return parts


@adapter("sc3d_mouse_embryo")
def _sc3d(root=ROOT_DEFAULT, work=None, max_sections=None):
    """sc3D マウス胚 (Slide-seq): GSM 毎に digital_expression.txt.gz (gene×cell) +
    matched_bead_locations.txt.gz (行=cell 順で col1=x col2=y)。GSM=切片。"""
    d=os.path.join(root,"sc3d_mouse_embryo")
    tar=glob.glob(d+"/*.tar")[0]
    work=work or os.path.join(root,"_adapt_work","sc3d")
    dges=_extract_members(tar, work, r'digital_expression\.txt\.gz$')
    beads=_extract_members(tar, work, r'matched_bead_locations\.txt\.gz$')
    def gsm(p): return re.search(r'(GSM\d+_[\d_]+?)[._]', os.path.basename(p)+'.').group(1)
    bead_by={gsm(b):b for b in beads}
    parts=[]
    for dg in (dges[:max_sections] if max_sections else dges):
        g=gsm(dg)
        if g not in bead_by: continue
        m=pd.read_csv(dg, sep='\t', index_col=0)      # gene x cell
        m=m.T                                          # cell x gene
        bead=pd.read_csv(bead_by[g], sep='\t', header=None)   # positional: [_, x, y]
        if bead.shape[0]!=m.shape[0]: continue
        a=ad.AnnData(X=m.values.astype('float32'))
        a.obs_names=[f"{g}_{c}" for c in m.index]; a.var_names=list(m.columns)
        a.obsm['spatial']=bead.iloc[:, [1,2]].values.astype(float)
        a.obs['section']=g; a.uns['_srcname']=g
        parts.append(a)
    return parts


@adapter("human_intestine_dev_spatial")
def _intestine_dev(root=ROOT_DEFAULT, work=None, max_sections=None):
    """腸発生 Visium: 外 RAW.tar 内に GSM_A?.tar.gz（各=1切片; raw_feature_bc_matrix/ +
    spatial/tissue_positions_list.csv）。二段展開して mtx+positions で組む。"""
    d=os.path.join(root,"human_intestine_dev_spatial")
    outer=glob.glob(d+"/*RAW.tar")[0]
    work=work or os.path.join(root,"_adapt_work","intestine_dev")
    inners=_extract_members(outer, work, r'GSM\d+_[A-Z]\d+\.tar\.gz$')
    parts=[]
    for inner in (inners[:max_sections] if max_sections else inners):
        sec=re.search(r'(GSM\d+_[A-Z]\d+)', os.path.basename(inner)).group(1)
        sw=os.path.join(work, sec)
        _extract_members(inner, sw, r'(raw|filtered)_feature_bc_matrix/.*(matrix\.mtx|barcodes|features)')
        _extract_members(inner, sw, r'tissue_positions_list\.csv$')
        mtxdir=glob.glob(sw+"/**/*feature_bc_matrix",recursive=True)
        pos=glob.glob(sw+"/**/tissue_positions_list.csv",recursive=True)
        if not mtxdir or not pos: continue
        try:
            a=_visium_part(mtx_dir=mtxdir[0], positions=pos[0], section=sec)
            # raw matrix は組織外spotを含む → 座標NaNでない spot に限定
            keep=~np.isnan(a.obsm['spatial']).any(axis=1)
            a=a[keep].copy()
            parts.append(a)
        except Exception: pass
    return parts


@adapter("mouse_olfactory_bulb_crossplatform")
def _olfactory(root=ROOT_DEFAULT, work=None, max_sections=None):
    """マウス嗅球 Stereo-seq: 単一切片。long-format Cell_GetExp_gene.txt (geneID,x,y,UMICount,label)
    を cell(label) × gene のマトリクスへ pivot。position.tsv に cell 座標。"""
    d=os.path.join(root,"mouse_olfactory_bulb_crossplatform")
    targz=glob.glob(d+"/*.tar.gz")[0]
    work=work or os.path.join(root,"_adapt_work","olfactory")
    exp=_extract_members(targz, work, r'Cell_GetExp_gene\.txt$')
    pos=_extract_members(targz, work, r'position\.tsv$')
    if not exp: return []
    long=pd.read_csv(exp[0], sep='\t')          # geneID,x,y,UMICount,label
    # pivot to cell(label) x gene
    mat=long.pivot_table(index='label', columns='geneID', values='UMICount', aggfunc='sum', fill_value=0)
    a=ad.AnnData(X=mat.values.astype('float32'))
    a.obs_names=[f"cell{int(l)}" for l in mat.index]; a.var_names=list(mat.columns)
    # coords: prefer position.tsv (label->x,y); else centroid of gene coords per label
    if pos:
        p=pd.read_csv(pos[0], sep='\t')
        p=p.set_index('label').reindex(mat.index)
        a.obsm['spatial']=p[['x','y']].values.astype(float)
    else:
        cen=long.groupby('label')[['x','y']].mean().reindex(mat.index)
        a.obsm['spatial']=cen.values.astype(float)
    a.obs['section']='OB'; a.uns['_srcname']='mouse_olfactory_bulb'
    return [a]


@adapter("human_fetal_cortex_merfish")
def _fetal_cortex(root=ROOT_DEFAULT, work=None, max_sections=None):
    """ヒト胎児皮質 MERFISH: zip 毎（=切片, FB080_F1 等）に cell_by_gene.csv (cell×gene) +
    cell_metadata.csv (EntityID, center_x, center_y)。"""
    d=os.path.join(root,"human_fetal_cortex_merfish")
    zips=sorted(glob.glob(d+"/*.zip"))
    work=work or os.path.join(root,"_adapt_work","fetal_cortex")
    parts=[]
    for zp in (zips[:max_sections] if max_sections else zips):
        sec=os.path.basename(zp)[:-4]
        sw=os.path.join(work, sec)
        cg=_extract_members(zp, sw, r'cell_by_gene\.csv$')
        cm=_extract_members(zp, sw, r'cell_metadata\.csv$')
        if not (cg and cm): continue
        m=pd.read_csv(cg[0], index_col=0)          # cell x gene
        meta=pd.read_csv(cm[0], index_col=0).reindex(m.index)
        xc=[c for c in meta.columns if c.lower()=='center_x'][0]
        yc=[c for c in meta.columns if c.lower()=='center_y'][0]
        a=ad.AnnData(X=m.values.astype('float32'))
        a.obs_names=[f"{sec}_{i}" for i in m.index]; a.var_names=list(m.columns)
        a.obsm['spatial']=meta[[xc,yc]].values.astype(float)
        a.obs['section']=sec; a.uns['_srcname']=sec
        parts.append(a)
    return parts


@adapter("liver_zonation_st")
def _liver(root=ROOT_DEFAULT, work=None, max_sections=None):
    """肝ゾーネーション ST: RAW.tar 内 GSM_Liver_*_stdata.tsv.gz（spot×gene, spot名 'AxB'）。
    GSM=切片。座標は spot 名から。"""
    d=os.path.join(root,"liver_zonation_st")
    tar=glob.glob(d+"/*RAW.tar")[0]
    work=work or os.path.join(root,"_adapt_work","liver")
    std=_extract_members(tar, work, r'_stdata\.tsv\.gz$')
    parts=[]
    for t in (std[:max_sections] if max_sections else std):
        sec=re.search(r'(GSM\d+)', os.path.basename(t)).group(1)
        m=pd.read_csv(t, sep='\t', index_col=0)
        xy=_spotname_xy(m.index)
        if xy is None: continue
        a=ad.AnnData(X=m.values.astype('float32'))
        a.obs_names=[f"{sec}_{s}" for s in m.index]; a.var_names=list(m.columns)
        a.obsm['spatial']=xy; a.obs['section']=sec; a.uns['_srcname']=sec
        parts.append(a)
    return parts


@adapter("human_placenta_spatial_multiomic")
def _placenta(root=ROOT_DEFAULT, work=None, max_sections=None):
    """胎盤 STARmap: sample(W7/W11…) 毎に raw_expression.csv (cell×gene) +
    cell_metadata.csv (x_um,y_um, cell_id)。sample=切片。"""
    d=os.path.join(root,"human_placenta_spatial_multiomic")
    exprs=sorted(glob.glob(d+"/STARmap-ISS_sample_*_raw_expression.csv"))
    parts=[]
    for ex in (exprs[:max_sections] if max_sections else exprs):
        sec=re.search(r'sample_([A-Za-z0-9]+)_raw', os.path.basename(ex)).group(1)
        cm=ex.replace('_raw_expression.csv','_cell_metadata.csv')
        if not os.path.exists(cm): continue
        m=pd.read_csv(ex, index_col=0)          # cell x gene
        meta=pd.read_csv(cm)
        idcol=[c for c in meta.columns if c.lower()=='cell_id'][0]
        meta=meta.set_index(idcol).reindex(m.index)
        xc=[c for c in meta.columns if c.lower() in ('x_um','x')][0]
        yc=[c for c in meta.columns if c.lower() in ('y_um','y')][0]
        a=ad.AnnData(X=m.values.astype('float32'))
        a.obs_names=[f"{sec}_{i}" for i in m.index]; a.var_names=list(m.columns)
        a.obsm['spatial']=meta[[xc,yc]].values.astype(float)
        a.obs['section']=sec; a.uns['_srcname']=sec
        for anno in ('cell_type','celltype','annotation'):
            if anno in meta.columns: a.obs['annotation']=meta[anno].astype(str).values; break
        parts.append(a)
    return parts


@adapter("moffitt_hypothalamus_merfish")
def _moffitt(root=ROOT_DEFAULT, work=None, max_sections=None):
    """Moffitt 視床下部 MERFISH（この配布版は2Dのみ）: merfish_expression.csv (cell×gene) +
    merfish_spatial.xlsx (x_microns,y_microns)。z(bregma)はこの束に無く単一切片扱い。"""
    d=os.path.join(root,"moffitt_hypothalamus_merfish")
    work=work or os.path.join(root,"_adapt_work","moffitt")
    zp=glob.glob(d+"/*.zip")[0]
    _extract_members(zp, work, r'merfish_expression\.csv$')
    _extract_members(zp, work, r'merfish_spatial\.xlsx$')
    ex=glob.glob(work+"/**/merfish_expression.csv",recursive=True)
    sp=glob.glob(work+"/**/merfish_spatial.xlsx",recursive=True)
    if not (ex and sp): return []
    m=pd.read_csv(ex[0], index_col=0)
    coo=pd.read_excel(sp[0], index_col=0)
    n=min(len(m), len(coo))
    a=ad.AnnData(X=m.iloc[:n].values.astype('float32'))
    a.obs_names=[f"cell{i}" for i in range(n)]; a.var_names=list(m.columns)
    xc=[c for c in coo.columns if 'x' in c.lower()][0]; yc=[c for c in coo.columns if 'y' in c.lower()][0]
    a.obsm['spatial']=coo.iloc[:n][[xc,yc]].values.astype(float)
    a.obs['section']='A1'; a.uns['_srcname']='moffitt_2d'
    return [a]


def _visium_mtx_by_gsm(dataset_id, root, work, max_sections, sample_rx):
    """GEO Visium: GSM 毎に <prefix>_matrix.mtx.gz + barcodes + features + tissue_positions(_list).csv.gz。
    ファイル名の共通プレフィックス（GSM+サンプル）でグループ化し、各グループ=1切片。"""
    d=os.path.join(root, dataset_id)
    tar=[f for f in glob.glob(d+"/**/*",recursive=True) if f.endswith(('.tar','.tar.gz'))][0]
    work=work or os.path.join(root,"_adapt_work",dataset_id)
    mats=_extract_members(tar, work, r'matrix\.mtx\.gz$')
    def prefix(p):
        b=os.path.basename(p)
        return re.sub(r'[_-]?matrix\.mtx\.gz$','',b)
    parts=[]
    for mtx in (mats[:max_sections] if max_sections else mats):
        pref=prefix(mtx)
        # find sibling barcodes/features/positions sharing the prefix (allow _ or - joiner)
        stem=re.escape(pref)
        bc=_extract_members(tar, work, stem+r'[_-]?barcodes\.tsv\.gz$')
        ft=_extract_members(tar, work, stem+r'[_-]?features\.tsv\.gz$')
        pos=_extract_members(tar, work, stem+r'[_-]?tissue_positions(_list)?\.csv\.gz$')
        if not (bc and ft and pos): continue
        # build AnnData from triplet
        a=sc.read_mtx(mtx).T
        fdf=pd.read_csv(ft[0],sep='\t',header=None)
        names=(fdf[1] if fdf.shape[1]>1 else fdf[0]).astype(str).values[:a.n_vars]
        a.var_names=pd.Index(names).astype(str)
        a.obs_names=pd.read_csv(bc[0],sep='\t',header=None)[0].astype(str).values[:a.n_obs]
        a.var_names_make_unique()
        p=_read_tissue_positions(pos[0]).reindex(a.obs_names)
        a.obsm['spatial']=p[['X','Y']].values.astype(float)
        sec_m=re.search(sample_rx, pref) if sample_rx else None
        sec=sec_m.group(1) if sec_m else pref   # フォールバックは完全プレフィックス（切片一意）
        a.obs['section']=sec; a.uns['_srcname']=sec
        parts.append(a)
    return parts


@adapter("aging_mouse_brain_visium")
def _aging_brain(root=ROOT_DEFAULT, work=None, max_sections=None):
    """老化マウス脳 Visium: サンプル毎 mtx トリプレット + tissue_positions。ファイルprefix=切片。"""
    return _visium_mtx_by_gsm("aging_mouse_brain_visium", root, work, max_sections, None)


@adapter("human_hippocampus_visium")
def _hippocampus(root=ROOT_DEFAULT, work=None, max_sections=None):
    """ヒト海馬 Visium: capture area 毎（GSM_slide_A1/B1/…）mtx トリプレット + positions。prefix=切片。"""
    return _visium_mtx_by_gsm("human_hippocampus_visium", root, work, max_sections, None)


@adapter("human_fetal_lung_visium")
def _lung(root=ROOT_DEFAULT, work=None, max_sections=None):
    """胎児肺 Visium: tar.gz 毎に1サンプル（<id>/outs/filtered_feature_bc_matrix/ + outs/spatial）。"""
    d=os.path.join(root,"human_fetal_lung_visium")
    tars=[f for f in glob.glob(d+"/**/*",recursive=True) if f.endswith(('.tar','.tar.gz'))]
    work=work or os.path.join(root,"_adapt_work","lung")
    parts=[]
    for t in (tars[:max_sections] if max_sections else tars):
        sec=re.search(r'(\d+STDY\d+)', os.path.basename(t))
        sec=sec.group(1) if sec else os.path.basename(t).split('.')[0]
        sw=os.path.join(work, sec)
        _extract_members(t, sw, r'filtered_feature_bc_matrix/.*(matrix\.mtx|barcodes|features)')
        _extract_members(t, sw, r'tissue_positions_list\.csv$')
        mtxdir=glob.glob(sw+"/**/filtered_feature_bc_matrix",recursive=True)
        pos=glob.glob(sw+"/**/tissue_positions_list.csv",recursive=True)
        if not mtxdir or not pos: continue
        try:
            parts.append(_visium_part(mtx_dir=mtxdir[0], positions=pos[0], section=sec))
        except Exception: pass
    return parts


# ==================================================================
# 追加アダプタ（ロードエラー修正 — 2026-07 catalog全件対応）
# ==================================================================
def _zst_to(src, dstdir):
    """*.zst を dstdir に展開してパスを返す。"""
    os.makedirs(dstdir, exist_ok=True)
    out=os.path.join(dstdir, os.path.basename(src)[:-4])
    if not os.path.exists(out) or os.path.getsize(out)==0:
        subprocess.run(["zstd","-d","-q","-f","-o",out,src], check=True)
    return out


@adapter("flysta3d_v2_drosophila")
def _flysta3d_v2(root=ROOT_DEFAULT, work=None, max_sections=None):
    """Flysta3D v2 (Cell 2025): 各 *.h5ad.zst が1胚/1時点。obsm['align_spatial'] が
    pre-aligned 3D 座標、obs['align_z'] が胚内の切片。annotation あり。
    max_sections はファイル(=胚)数の上限として解釈する。"""
    d=os.path.join(root,"flysta3d_v2_drosophila")
    work=work or os.path.join(root,"_adapt_work","flysta3d_v2")
    zs=sorted(f for f in glob.glob(d+"/*.h5ad.zst"))
    if max_sections: zs=zs[:max_sections]
    parts=[]
    for z in zs:
        emb=os.path.basename(z).replace(".h5ad.zst","")
        h=_zst_to(z, work)
        a=ad.read_h5ad(h)
        # 座標キーのフォールバック: align_spatial(3D) -> spatial(2D)+align_z -> raw_spatial+align_z
        if 'align_spatial' in a.obsm and np.asarray(a.obsm['align_spatial']).shape[1]>=3:
            xyz=np.asarray(a.obsm['align_spatial'], float)
        else:
            base=None
            for k in ('align_spatial','spatial','raw_spatial'):
                if k in a.obsm and np.asarray(a.obsm[k]).shape[1]>=2:
                    base=np.asarray(a.obsm[k], float)[:,:2]; break
            if base is None or 'align_z' not in a.obs:
                try: os.remove(h)
                except: pass
                continue
            zc=np.asarray(a.obs['align_z'], float).reshape(-1,1)
            xyz=np.hstack([base, zc])
        if xyz.shape[1]<3: continue
        # keep only 'spatial' obsm — other keys (X_pca/X_umap/contour/bbox) have
        # per-embryo-inconsistent shapes and break outer ad.concat across files
        for k in list(a.obsm.keys()):
            if k!='spatial': del a.obsm[k]
        a.obsm['spatial']=xyz
        zc=xyz[:,2]
        a.obs['section']=[f"{emb}_z{int(round(v))}" for v in zc]
        if 'annotation' not in a.obs and 'Annotation_2_tissue' in a.obs:
            a.obs['annotation']=a.obs['Annotation_2_tissue'].astype(str)
        if 'annotation' not in a.obs:
            a.obs['annotation']=pd.Series(['NA']*a.n_obs,index=a.obs_names)
        a.obs['annotation']=a.obs['annotation'].astype(str)
        a.uns['_srcname']=emb
        parts.append(a)
        try: os.remove(h)
        except: pass
    return parts


@adapter("cerebellum_crossspecies_spatial")
def _cerebellum(root=ROOT_DEFAULT, work=None, max_sections=None):
    """小脳cross-species (Hao et al.): 展開済み225 h5adに複数種(Macaque/Marmoset/Mouse)×
    複数個体×多数切片が混在。3D連続切片再構成には単一個体の連番切片を使う。切片数最多の
    Macaque1 (T40..T110の35枚)を採用。annotation は小脳層(molecular/granular/purkinje/white)。
    z は切片T番号。33万細胞/切片と重いので max_sections で切片数を絞る(連番の等間隔サンプル)。"""
    import re
    cands=["/mnt/d/3DCurv_tmp/_work/cerebellum_crossspecies_spatial",
           os.path.join(root,"_work","cerebellum_crossspecies_spatial"),
           os.path.join(root,"cerebellum_crossspecies_spatial")]
    d=next((c for c in cands if os.path.isdir(c)), cands[0])
    pat=re.compile(r'^Macaque1_T(\d+)\.h5ad$')
    sel=sorted([(int(pat.match(f).group(1)),f) for f in os.listdir(d) if pat.match(f)])
    if max_sections and len(sel)>max_sections:
        idx=np.linspace(0,len(sel)-1,max_sections).round().astype(int)
        sel=[sel[i] for i in idx]
    parts=[]
    for tnum,f in sel:
        a=ad.read_h5ad(os.path.join(d,f))
        if 'annotation' not in a.obs or 'spatial' not in a.obsm: continue
        xy=np.asarray(a.obsm['spatial'],float)[:,:2]
        z=np.full((a.n_obs,1), float(tnum))
        ann=a.obs['annotation'].astype(str).values
        import anndata as _ad
        b=_ad.AnnData(X=a.X, var=a.var.copy())
        b.obs_names=a.obs_names; b.obs['annotation']=ann
        b.obs['section']=f"Macaque1_T{tnum}"
        b.obsm['spatial']=np.hstack([xy,z]); b.uns['_srcname']=f"Macaque1_T{tnum}"
        parts.append(b); del a
    return parts


@adapter("mosta_mouse_embryo")
def _mosta(root=ROOT_DEFAULT, work=None, max_sections=None):
    """MOSTA (Chen et al. Cell 2022): 展開済みキャッシュに全ステージ×全胚×全切片の
    h5ad が65個(201GB)混在する。3D連続切片再構成には『同一ステージ・同一胚の連続切片』
    だけを使う。切片数最多で3Dに適した E16.5 胚E2 (S1..S13) を採用。
    各切片は obs に Regulon/Module 列を400以上持ち重いので annotation+spatial だけ残す。
    z は切片番号(S番号)。max_sections は切片数上限。"""
    import re
    # キャッシュ場所(展開済み)。workdir 直下か _work 配下
    cands=[os.path.join(root,"_work","mosta_mouse_embryo"),
           "/mnt/d/3DCurv_tmp/_work/mosta_mouse_embryo",
           os.path.join(root,"mosta_mouse_embryo")]
    d=next((c for c in cands if os.path.isdir(c)), cands[0])
    allf=[f for f in os.listdir(d) if f.endswith('.h5ad')]
    # E16.5_E2 の連続切片 (E16.5_E2S<番号>.MOSTA.h5ad) を番号順に
    pat=re.compile(r'^E16\.5_E2S(\d+)\.MOSTA\.h5ad$')
    sel=sorted([(int(pat.match(f).group(1)),f) for f in allf if pat.match(f)])
    if max_sections: sel=sel[:max_sections]
    parts=[]
    for snum,f in sel:
        a=ad.read_h5ad(os.path.join(d,f))
        if 'annotation' not in a.obs or 'spatial' not in a.obsm: continue
        xy=np.asarray(a.obsm['spatial'],float)[:,:2]
        z=np.full((a.n_obs,1), float(snum))
        # 重い obs 列(Regulon-*/Module_*/QC)を捨て annotation だけ残す
        ann=a.obs['annotation'].astype(str).values
        import anndata as _ad
        b=_ad.AnnData(X=a.X, var=a.var.copy())
        b.obs_names=a.obs_names; b.obs['annotation']=ann
        b.obs['section']=f"E16.5_E2S{snum}"
        b.obsm['spatial']=np.hstack([xy, z])
        b.uns['_srcname']=f"E16.5_E2S{snum}"
        parts.append(b); del a
    return parts


@adapter("whole_mouse_embryo_3d_cngb")
def _cngb(root=ROOT_DEFAULT, work=None, max_sections=None):
    """CNGB 全胚 3D (GSE237308): counts.tsv(gene x cell) + coords.csv(section,imagerow/col) +
    meta.csv(spatial_domain)。E135 sec1_sec10 セット。"""
    d=os.path.join(root,"whole_mouse_embryo_3d_cngb")
    cnt=glob.glob(d+"/*sec1_sec10_count_mtx.tsv.gz")
    crd=glob.glob(d+"/*sec1_sec10_sections_coordinates_merge.csv.gz")
    meta=glob.glob(d+"/*sec1_sec10_meta.csv.gz")
    if not (cnt and crd and meta): return []
    coords=_read_csv_any(crd[0], index_col=0)
    md=_read_csv_any(meta[0], index_col=0)
    cm=_read_csv_any(cnt[0], index_col=0)
    # strip stray double-quotes from labels (this GEO deposit quotes every field)
    _dq=lambda s: str(s).strip().strip('"')
    cm.index=[_dq(x) for x in cm.index]; cm.columns=[_dq(x) for x in cm.columns]
    coords.index=[_dq(x) for x in coords.index]; md.index=[_dq(x) for x in md.index]
    X=cm.T
    cells=[c for c in X.index if c in coords.index]
    X=X.loc[cells]
    a=ad.AnnData(X.values.astype('float32'))
    a.obs_names=list(cells); a.var_names=list(X.columns)
    a.obsm['spatial']=coords.loc[cells,['imagecol','imagerow']].values.astype(float)
    a.obs['section']=coords.loc[cells,'section'].astype(str).values
    if 'spatial_domain' in md.columns:
        a.obs['annotation']=md.reindex(cells)['spatial_domain'].astype(str).values
    return [a]


def _allen_common(root, sub, work, max_sections, per_section_cap=1500):
    """Allen ABC atlas MERFISH: expression_matrices/<brain>/<date>/<brain>-log2.h5ad.zst +
    metadata/<brain>/*/cell_metadata.csv (x,y,z,brain_section_label) +
    metadata/<brain>/*/views/cell_metadata_with_cluster_annotation.csv (class/subclass)。
    座標・切片は cell_metadata から取得（reconstructed_coordinates は 638850 のみ）。
    数百万細胞なので brain_section_label ごとに per_section_cap で層別サブサンプル。
    max_sections=切片数上限（1切片あたり少数残す）。"""
    d=os.path.join(root, sub)
    parts=[]
    exps=sorted(glob.glob(d+"/expression_matrices/*/*/*-log2.h5ad.zst"))
    for exp in exps:
        br=os.path.basename(os.path.dirname(os.path.dirname(exp)))  # dir name (metadata key)
        cm =glob.glob(d+"/metadata/"+br+"/*/cell_metadata.csv")
        ann=glob.glob(d+"/metadata/"+br+"/*/views/cell_metadata_with_cluster_annotation.csv")
        if not cm: continue
        mdf=pd.read_csv(cm[0], index_col=0); mdf.index=mdf.index.astype(str)
        if not set(['x','y']).issubset(mdf.columns): continue
        seccol='brain_section_label' if 'brain_section_label' in mdf.columns else None
        # pick sections + stratified subsample by section (vectorized groupby.sample)
        if seccol:
            secs=sorted(mdf[seccol].dropna().astype(str).unique())
            if max_sections: secs=secs[:max_sections]
            mdf=mdf[mdf[seccol].astype(str).isin(secs)]
            mdf=mdf.groupby(seccol, group_keys=False).apply(
                lambda g: g.sample(min(len(g), per_section_cap), random_state=0))
        h=_zst_to(exp, os.path.join(work, br))
        a=ad.read_h5ad(h); a.obs_names=a.obs_names.astype(str)
        # subsample h5ad to the chosen cells BEFORE copy (vectorized isin)
        sel=a.obs_names.isin(mdf.index)
        a=a[sel].copy()
        mdf=mdf.reindex(a.obs_names)          # align to h5ad order, drop non-overlap
        keep=mdf[['x','y']].notna().all(axis=1).values
        a=a[keep].copy(); mdf=mdf.loc[a.obs_names]
        # true 3D coords: use z if present (finer than section-index z)
        cols=['x','y','z'] if 'z' in mdf.columns else ['x','y']
        a.obsm['spatial']=mdf[cols].values.astype(float)
        a.obs['section']=(mdf[seccol].astype(str).values if seccol else br)
        if ann:
            adf=pd.read_csv(ann[0], index_col=0); adf.index=adf.index.astype(str)
            for col in ['class','subclass','parcellation_structure','parcellation_division']:
                if col in adf.columns:
                    a.obs['annotation']=adf.reindex(a.obs_names)[col].astype(str).values; break
        a.uns['_srcname']=br
        parts.append(a)
        try: os.remove(h)
        except: pass
        break  # one brain per load (each brain is a full 3D volume)
    return parts


@adapter("allen_merfish_638850")
def _allen_638850(root=ROOT_DEFAULT, work=None, max_sections=None):
    work=work or os.path.join(root,"_adapt_work","allen_638850")
    return _allen_common(root, "allen_merfish_638850", work, max_sections)


@adapter("allen_merfish_zhuang")
def _allen_zhuang(root=ROOT_DEFAULT, work=None, max_sections=None):
    work=work or os.path.join(root,"_adapt_work","allen_zhuang")
    return _allen_common(root, "allen_merfish_zhuang", work, max_sections)


@adapter("adult_human_small_intestine_spatial")
def _small_intestine(root=ROOT_DEFAULT, work=None, max_sections=None):
    """成人小腸 (Nature 2024): cells.csv(centroid_x/y, cell_type, sample) + raw_mat.csv(gene x cell)。"""
    d=os.path.join(root,"adult_human_small_intestine_spatial")
    cells=pd.read_csv(os.path.join(d,"cells.csv"), index_col=0)
    mat=pd.read_csv(os.path.join(d,"raw_mat.csv"), index_col=0)
    X=mat.T
    common=[c for c in X.index if c in cells.index]
    X=X.loc[common]; cd=cells.loc[common]
    if max_sections:
        keep=sorted(cd['sample'].astype(str).unique())[:max_sections]
        m=cd['sample'].astype(str).isin(keep)
        X=X.loc[m.values]; cd=cd.loc[m.values]; common=list(cd.index)
    a=ad.AnnData(X.values.astype('float32'))
    a.obs_names=list(common); a.var_names=list(X.columns)
    a.obsm['spatial']=cd[['centroid_x','centroid_y']].values.astype(float)
    a.obs['section']=cd['sample'].astype(str).values
    a.obs['annotation']=cd['cell_type'].astype(str).values
    return [a]


@adapter("dlpfc_visium_spatiallibd")
def _dlpfc(root=ROOT_DEFAULT, work=None, max_sections=None):
    """DLPFC spatialLIBD (Nat Neurosci 2021): DLPFC12.zip 内に 12 切片、各 <id>/ に
    filtered_feature_bc_matrix.h5 + spatial/tissue_positions_list.csv。皮質層(Layer1-6/WM)。"""
    import scanpy as sc
    d=os.path.join(root,"dlpfc_visium_spatiallibd")
    zf=glob.glob(d+"/*.zip")[0]
    work=work or os.path.join(root,"_adapt_work","dlpfc")
    with zipfile.ZipFile(zf) as z:
        names=[n for n in z.namelist() if '__MACOSX' not in n and not os.path.basename(n).startswith('._')]
    secs=sorted({n.split('/')[1] for n in names if n.startswith('DLPFC12/') and len(n.split('/'))>2 and n.split('/')[1]})
    if max_sections: secs=secs[:max_sections]
    parts=[]
    for sec in secs:
        _extract_members(zf, work, r'DLPFC12/'+sec+r'/'+sec+r'_filtered_feature_bc_matrix\.h5$')
        _extract_members(zf, work, r'DLPFC12/'+sec+r'/spatial/tissue_positions.*\.csv$')
        h5=glob.glob(work+"/DLPFC12/"+sec+"/*filtered_feature_bc_matrix.h5")
        pos=glob.glob(work+"/DLPFC12/"+sec+"/spatial/tissue_positions*.csv")
        if not h5: continue
        try:
            p=_visium_part(h5=h5[0], positions=(pos[0] if pos else None), section=sec)
            parts.append(p)
        except Exception: pass
    return parts


@adapter("iriseq_aging_brain")
def _iriseq_v2(root=ROOT_DEFAULT, work=None, max_sections=None):
    """IRISeq 老化脳 (GSE270383): 全体で1セット。count.mtx.gz(gene x cell) + barcodes/genes +
    meta_data.csv.gz(UMAP1_spatial/UMAP2_spatial 座標, Annotation=領域_sectionN)。"""
    import scipy.io
    d=os.path.join(root,"iriseq_aging_brain")
    work=work or os.path.join(root,"_adapt_work","iriseq")
    mtx=_gunzip_to(glob.glob(d+"/*count.mtx.gz")[0], work)
    bc=_read_csv_any(glob.glob(d+"/*barcodes.tsv.gz")[0], header=None)[0].astype(str).values
    gn=_read_csv_any(glob.glob(d+"/*genes.tsv.gz")[0], header=None)
    genes=gn[gn.columns[-1]].astype(str).values
    md=_read_csv_any(glob.glob(d+"/*meta_data.csv.gz")[0], index_col=0)
    md.index=md.index.astype(str)
    M=scipy.io.mmread(mtx).tocsr()
    X=M.T.tocsr()
    a=ad.AnnData(X)
    a.obs_names=bc[:a.n_obs]; a.var_names=genes[:a.n_vars]
    a.var_names_make_unique()
    md=md.reindex(a.obs_names)
    xcol=[c for c in md.columns if c.lower()=='umap1_spatial'][0]
    ycol=[c for c in md.columns if c.lower()=='umap2_spatial'][0]
    a.obsm['spatial']=md[[xcol,ycol]].values.astype(float)
    ann=md['Annotation'].astype(str)
    sec=ann.str.extract(r'(section[\d_]+)$')[0].fillna('s0')
    reg=ann.str.replace(r'_?section[\d_]+$','',regex=True)
    a.obs['section']=sec.values
    a.obs['annotation']=reg.values
    keep=np.isfinite(a.obsm['spatial']).all(1)
    a=a[keep].copy()
    if max_sections:
        ks=sorted(pd.unique(a.obs['section']))[:max_sections]
        a=a[a.obs['section'].isin(ks)].copy()
    try: os.remove(mtx)
    except: pass
    return [a]


@adapter("her2_breast_andersson")
def _her2(root=ROOT_DEFAULT, work=None, max_sections=None):
    """HER2+ breast (Andersson Nat Commun 2021): Her2_tumor.zip 内
    ST-cnts/<sec>.tsv (spot x gene; spot id='COLxROW') + ST-pat/lbl/<sec>_labeled_coordinates.tsv
    (x,y,label=invasive cancer/immune infiltrate/... の一部切片のみ)。患者A-H×複数切片。
    label があれば annotation に、無ければ座標のみ。max_sections=切片数上限。"""
    d=os.path.join(root,"her2_breast_andersson")
    zf=glob.glob(d+"/*.zip")[0]
    work=work or os.path.join(root,"_adapt_work","her2")
    with zipfile.ZipFile(zf) as z:
        cnts=sorted(n for n in z.namelist()
                    if re.search(r'ST-cnts/[A-H]\d+\.tsv$', n) and '__MACOSX' not in n)
    secs=[os.path.basename(n)[:-4] for n in cnts]
    if max_sections: secs=secs[:max_sections]
    parts=[]
    for sec in secs:
        _extract_members(zf, work, r'ST-cnts/'+sec+r'\.tsv$')
        _extract_members(zf, work, r'ST-pat/lbl/'+sec+r'_labeled_coordinates\.tsv$')
        cp=glob.glob(work+"/Her2_tumor/ST-cnts/"+sec+".tsv")
        if not cp: continue
        cm=pd.read_csv(cp[0], sep='\t', index_col=0)
        # spot id 'COLxROW' -> coords
        xy=np.array([[float(s.split('x')[0]), float(s.split('x')[1])] for s in cm.index])
        a=ad.AnnData(cm.values.astype('float32'))
        a.obs_names=list(cm.index); a.var_names=list(cm.columns)
        a.obsm['spatial']=xy
        a.obs['section']=sec; a.uns['_srcname']=sec
        lp=glob.glob(work+"/Her2_tumor/ST-pat/lbl/"+sec+"_labeled_coordinates.tsv")
        if lp:
            lab=pd.read_csv(lp[0], sep='\t')
            # match by rounded array coords
            lab['key']=lab['x'].round().astype(int).astype(str)+"x"+lab['y'].round().astype(int).astype(str)
            spotkey=[str(int(round(float(s.split('x')[0]))))+"x"+str(int(round(float(s.split('x')[1])))) for s in cm.index]
            lmap=dict(zip(lab['key'], lab['label'].astype(str)))
            a.obs['annotation']=[lmap.get(k,'unlabeled') for k in spotkey]
        parts.append(a)
    return parts


@adapter("human_placenta_spatial_multiomic")
def _placenta_v2(root=ROOT_DEFAULT, work=None, max_sections=None):
    """胎盤 STARmap 空間マルチオミクス: sample(W7/W11) 毎に
    cell_metadata.csv(cell_id,x,y,sample) + raw_expression.csv(cell x gene)。sample=切片。"""
    d=os.path.join(root,"human_placenta_spatial_multiomic")
    cms=sorted(glob.glob(d+"/STARmap-ISS_sample_*_cell_metadata.csv"))
    if max_sections: cms=cms[:max_sections]
    parts=[]
    for cmp in cms:
        samp=re.search(r'sample_([^_]+)_cell_metadata', os.path.basename(cmp)).group(1)
        exp=glob.glob(d+f"/STARmap-ISS_sample_{samp}_raw_expression.csv")
        if not exp: continue
        cm=pd.read_csv(cmp)
        if 'cell_id' in cm.columns: cm=cm.set_index('cell_id')
        X=pd.read_csv(exp[0], index_col=0)
        common=[c for c in X.index if c in cm.index]
        X=X.loc[common]; cd=cm.loc[common]
        a=ad.AnnData(X.values.astype('float32'))
        a.obs_names=list(common); a.var_names=list(X.columns)
        xc='x' if 'x' in cd.columns else 'x_um'; yc='y' if 'y' in cd.columns else 'y_um'
        a.obsm['spatial']=cd[[xc,yc]].values.astype(float)
        a.obs['section']=f"W{samp}" if not str(samp).startswith('W') else str(samp)
        a.uns['_srcname']=str(samp)
        parts.append(a)
    return parts
