# Data sources

This study analyzes **only publicly available** spatial-transcriptomics datasets.
Raw `.h5ad` files are **not** redistributed in this repository (size and upstream
licensing). Each dataset below lists its accession, DOI, and download URL so the
analysis can be reproduced from primary sources. Only the derived result tables and
figures are included here (see `results/` and `figures/`).

The open interface sheet is extracted directly from each dataset's labelled 3D
point cloud by `src/open_interface.py` (no pre-built mesh is required). All
datasets are Stereo-seq unless noted.

| # | Dataset (this repo) | Name | Species | Technology | Accession | DOI | Source |
|---|---|---|---|---|---|---|---|
| 1 | `openst_lymphnode_3d` | Open-ST human metastatic lymph node 3D (Schott 2024) | human | Open-ST | GSE251926 / PRJNA1097649 | [10.1016/j.cell.2024.05.055](https://doi.org/10.1016/j.cell.2024.05.055) | [GEO / BioProject](https://rajewsky-lab.github.io/openst/) |
| 2 | `cerebellum_crossspecies_spatial` | Cross-species cerebellar cortex spatial atlas (mouse/marmoset/macaque) | mouse + marmoset + macaque | Stereo-seq + snRNA | CNP0003779 / CNP0005746 | [10.1126/science.ado3927](https://doi.org/10.1126/science.ado3927) | [verify (CNGB/GSA)](https://www.science.org/doi/10.1126/science.ado3927) |
| 3 | `artista_axolotl_brain` | ARTISTA: Axolotl telencephalon development/regeneration | axolotl | Stereo-seq | CNP0002068 / STDS0000056 | [10.1126/science.abp9444](https://doi.org/10.1126/science.abp9444) | [STOmicsDB/CNGB](https://db.cngb.org/stomics/artista/) |
| 4 | `zesta_zebrafish` | ZESTA: Zebrafish Embryogenesis Spatiotemporal Transcriptomic Atlas | zebrafish | Stereo-seq | CNP0002220 / STDS0000057 | [10.1016/j.devcel.2022.04.009](https://doi.org/10.1016/j.devcel.2022.04.009) | [STOmicsDB/CNGB](https://db.cngb.org/stomics/) |
| 5 | `digital_mouse_embryo_seu3d` | Digital reconstruction of full mouse embryos E7.5-E8.0 (SEU-3D, Xie 2025) | mouse | Stereo-seq | GSE278603 | [10.1016/j.cell.2025.05.035](https://doi.org/10.1016/j.cell.2025.05.035) | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE278603) |
| 6 | `mosta_mouse_embryo` | MOSTA: Mouse Organogenesis Spatiotemporal Transcriptomic Atlas | mouse | Stereo-seq | CNP0001543 / STDS0000058 | [10.1016/j.cell.2022.04.003](https://doi.org/10.1016/j.cell.2022.04.003) | [STOmicsDB/CNGB](https://db.cngb.org/stomics/mosta/) |
| 7 | `whole_mouse_embryo_3d_cngb` | 3D transcriptomics maps of whole mouse embryo during organogenesis (bioRxiv 2024) | mouse | Stereo-seq | CNP0005981 | [10.1101/2024.08.17.608366](https://doi.org/10.1101/2024.08.17.608366) | [CNGB (CNSA)](https://www.biorxiv.org/content/10.1101/2024.08.17.608366v1) |
| 8 | `flysta3d_drosophila` | Flysta3D: 3D transcriptomic maps of Drosophila embryo/larva | drosophila | Stereo-seq | CNP0002189 / STDS0000060 | [10.1016/j.devcel.2022.04.006](https://doi.org/10.1016/j.devcel.2022.04.006) | [STOmicsDB/CNGB](https://db.cngb.org/stomics/flysta3d/) |
| 9 | `flysta3d_v2_drosophila` | Flysta3D-v2: Drosophila single-cell 3D multi-omics atlas (embryo-pupa) | drosophila | Stereo-seq + scRNA + scATAC | CNGB STDS0000398 (flysta3d-v2) + Mendeley 10.17632/tvvjfr3c6j.1 | [10.1016/j.cell.2025.05.047](https://doi.org/10.1016/j.cell.2025.05.047) | [STOmicsDB/CNGB + Mendeley](https://db.cngb.org/stomics/flysta3d-v2/) |
| 10 | `prista4d_planarian` | PRISTA4D: 4D single-cell spatiotemporal atlas of planarian regeneration | planaria (S. mediterranea) | Stereo-seq | CNGB STT0000028 | [10.1093/gigascience/giag064](https://doi.org/10.1093/gigascience/giag064) | [CNGB](https://db.cngb.org/stomics/prista4d/) |
| 11 | `acsta_arabidopsis` | ACSTA: Arabidopsis Cell-type-specific Spatiotemporal Transcriptomic Atlas | arabidopsis | Stereo-seq | STDS0000104 | [10.1016/j.devcel.2022.04.011](https://doi.org/10.1016/j.devcel.2022.04.011) | [STOmicsDB/CNGB](https://db.cngb.org/stomics/) |

## Licensing notes

- **CNGB / STOmicsDB** datasets are subject to CNGB database terms — see the portal
  for each accession before redistribution.
- **GEO** datasets (`digital_mouse_embryo_seu3d`, `openst_lymphnode_3d`) follow GEO
  terms. Note: Open-ST human tumor sequence-level data (GSE251926) may be
  controlled-access; this study used processed cell-level matrices only.
- The `cerebellum_crossspecies_spatial` license is listed as *verify* in the source
  catalog — confirm terms at the publisher before redistribution.

## Reconstruction notes

- **Serial sections are the upstream authors' own selection, not a contiguous
  z-stack.** For example, the Open-ST metastatic lymph node (`openst_lymphnode_3d`)
  ships **19 non-consecutive sections** (`n_section` = 2, 3, 4, 5, 6, 7, 9, 11,
  17, 18, 19, 23, 24, 25, 26, 28, 33, 34, 36); intervening section indices are
  absent from the deposited data, so z-spacing between reconstructed sections is
  uneven (≈29 µm typical, up to ≈174 µm across the widest gap). This is a
  property of the public dataset, not of this pipeline.
- **The interface is a genuine OPEN surface, not a closed solid.** The bulge/dent
  score is computed on the true domain / non-domain boundary. We voxelize both
  the domain and its non-domain neighborhood, extract an iso-surface of the
  domain-fraction field (marching cubes), and **keep only faces that have sampled
  tissue on both sides** — this discards the physical cut-faces of the acquired
  tissue block (the flat top/bottom of the serial-section slab), which are not a
  real tissue interface. The surface voxel resolution is **density-adaptive**
  (edge length sized to local cell density; no per-dataset hand-tuning). Bulge/
  dent is the normal displacement of this open sheet, with the sign fixed to the
  physical surface orientation (positive = convex = the domain protrudes into its
  neighbor).
- **Method correction (mid-project).** An earlier version of this analysis used
  the convex hull of the domain, which closes the acquired block into a solid and
  mislabels its cut-faces as "bulges." That flaw was caught and the entire
  analysis (126 interface–strata) was re-run on the open interface described
  above. The superseded closed-hull result tables (38 interfaces) are retained in
  `results/legacy_closed_mesh/` for provenance only. On thin-slab datasets (MOSTA,
  CNGB: only ~9–13 serial sections in Z) the open-sheet resolution — and hence
  sensitivity — is limited by the sparse Z sampling; this is a data limitation,
  documented rather than hidden.

This repository's own code, skill, and derived tables/figures are MIT-licensed
(see `LICENSE`); each raw dataset retains its upstream license.