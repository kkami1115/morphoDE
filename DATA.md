# Data sources & reproduction

This study analyzes **only publicly available** spatial-transcriptomics datasets.
Raw `.h5ad` files are **not** redistributed here (size and upstream licensing).
Only derived result tables and figures are included (`results/`, `figures/`).

Two files make the data layer reproducible:

- **`data/catalog.csv`** — machine-readable ledger of the 11 analysis datasets:
  `dataset_id, name, species, technology, n_sections, accession, doi, url,
  license, source, source_arg, source_arg2, access, notes`. The `source*` columns
  drive the downloader; `access` states whether a dataset is one-command
  downloadable (`auto`) or must be obtained through a portal/controlled channel
  (`manual`).
- **`scripts/download.py`** — catalog-driven downloader. Resumable, and verified:
  the first successful fetch of each file records its `sha256` in
  `data/<dataset_id>/CHECKSUMS.sha256`; re-running verifies existing files against
  that lockfile. Where the upstream repository publishes a checksum (Zenodo,
  Figshare), the fetched file is additionally checked against it.

Datasets download as **plain `.h5ad`** and the analysis reads plain `.h5ad`.
(`download.py --compress` is an optional local disk-space convenience only and is
not part of the pipeline.)

## Quick start

```bash
# See the catalog and each dataset's access class
python3 scripts/download.py --list

# Fetch one dataset into ./data/<dataset_id>/
python3 scripts/download.py mosta_mouse_embryo

# Fetch every automatically-downloadable dataset
python3 scripts/download.py --all

# Preview source URLs without downloading
python3 scripts/download.py flysta3d_drosophila --dry-run

# Re-check local files against the recorded sha256 lockfile
python3 scripts/download.py mosta_mouse_embryo --verify
```

Requires `python3` and `curl` (optional `zstd` only for `--compress`).

## The 11 analysis datasets

| # | `dataset_id` | Name | Species | Accession | Source | Access |
|---|---|---|---|---|---|---|
| 1 | `openst_lymphnode_3d` | Open-ST human metastatic lymph node 3D (Schott 2024) | human | GSE251926 / PRJNA1097649 | GEO | processed only¹ |
| 2 | `cerebellum_crossspecies_spatial` | Cross-species cerebellar cortex atlas | mouse/marmoset/macaque | CNP0003779 / CNP0005746 | CNGB/GSA | manual² |
| 3 | `artista_axolotl_brain` | ARTISTA: axolotl telencephalon dev/regen | axolotl | CNP0002068 / STDS0000056 | CNGB | auto |
| 4 | `zesta_zebrafish` | ZESTA: zebrafish embryogenesis atlas | zebrafish | CNP0002220 / STDS0000057 | CNGB | auto |
| 5 | `digital_mouse_embryo_seu3d` | Digital mouse embryos E7.5–E8.0 (SEU-3D, Xie 2025) | mouse | GSE278603 | GEO | auto |
| 6 | `mosta_mouse_embryo` | MOSTA: mouse organogenesis atlas | mouse | CNP0001543 / STDS0000058 | CNGB | auto |
| 7 | `whole_mouse_embryo_3d_cngb` | 3D whole mouse embryo organogenesis (bioRxiv 2024) | mouse | CNP0005981 | CNGB (CNSA) | manual² |
| 8 | `flysta3d_drosophila` | Flysta3D: Drosophila embryo/larva 3D maps | drosophila | CNP0002189 / STDS0000060 | CNGB | auto |
| 9 | `flysta3d_v2_drosophila` | Flysta3D-v2: Drosophila 3D multi-omics atlas | drosophila | STDS0000398 (+ Mendeley) | CNGB | auto |
| 10 | `prista4d_planarian` | PRISTA4D: planarian regeneration 4D atlas | planaria | STT0000028 | CNGB | auto |
| 11 | `acsta_arabidopsis` | ACSTA: Arabidopsis spatiotemporal atlas | arabidopsis | STDS0000104 | CNGB | auto |

¹ Open-ST human tumor **sequence-level** data (GSE251926) may be
controlled-access; this study used the **processed cell-level matrices** only.
`download.py` fetches the processed GEO supplementary files.

² **manual** = not fetchable by a single command. `cerebellum_crossspecies_spatial`
and `whole_mouse_embryo_3d_cngb` are CNGB **raw-data projects** (CNP/CNSA), not
STOmicsDB collection slugs, so there is no stable direct file listing; obtain them
through the CNGB portal per the publisher's data-availability statement.
`download.py <id>` prints the accession, DOI and portal URL for these and stops —
it will not produce a partial mirror.

## Licensing notes

- **CNGB / STOmicsDB** datasets are subject to CNGB database terms — check the
  portal for each accession before redistribution.
- **GEO** datasets (`digital_mouse_embryo_seu3d` is CC-BY; `openst_lymphnode_3d`
  follows GEO terms) — see each series page.
- `cerebellum_crossspecies_spatial` license is listed as **verify** in the
  catalog — confirm terms at the publisher before redistribution.
- `flysta3d_v2_drosophila` multi-omics supplements are partly on Mendeley Data
  (`doi:10.17632/tvvjfr3c6j.1`); the 3D Stereo-seq `.h5ad` used here come from the
  STOmicsDB collection.

This repository's own code, skill, and derived tables/figures are MIT-licensed
(see `LICENSE`); each raw dataset retains its upstream license.

## What "reproducible" means here

`download.py` fetches from the **original repositories**, so the exact bytes are
whatever those repositories currently serve. The per-dataset
`CHECKSUMS.sha256` lockfile is written on **your** first run — it makes your own
re-downloads verifiable and lets you detect a corrupted or upstream-changed file,
but it is not shipped pre-populated (we do not host the data, so we cannot certify
a canonical byte-for-byte hash for datasets served without upstream checksums).
Where the source repository exposes a per-file checksum in its API — **Zenodo and
Figshare** — the download is additionally checked against that md5 immediately.
GEO and CNGB do not expose per-file checksums this way, so those downloads are
covered by the sha256 lockfile (self-consistency across your own re-runs) rather
than an upstream hash.
