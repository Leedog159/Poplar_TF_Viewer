# Poplar TF Activation Domain Atlas

Interactive, genome-wide access to transcriptional **activation domains (ADs)** across
*Populus trichocarpa* transcription factors, with comparative analyses against *Arabidopsis thaliana*
orthologs. The app shows per-residue ADHunter activity predictions, empirical yeast sort-seq tile
scores, disorder, and Pfam domains for single TFs. It also compares ADs across orthogroups and
exports AD calls, sequences, and full datasets.

## Live access

**[Live Web App](https://share.streamlit.io/...)** (link coming soon)

## Features

| Page | What it does |
|---|---|
| **TF Viewer** | One TF at a time: predicted activity, optional empirical tiles, AD and Pfam interval lanes, disorder, plus copyable protein/CDS sequences for every AD and domain. Filter by species and TF family. |
| **Orthogroup Viewer** | All poplar and Arabidopsis orthologs/paralogs in a TF's orthogroup, stacked on a shared scale. Tracks can use MSA columns or ungapped residue positions. Each gene links to the TF Viewer. |
| **Bulk Lookup & Export** | Paste or upload a list of gene IDs and download a summary CSV plus FASTA files: AD peptides, AD CDS, full proteins, full CDS. |
| **About** | Citation, methods, score interpretation, full-dataset downloads, and contacts. |

AD calls: contiguous runs with predicted activity **≥ 1.0 (z-score) for ≥ 10 amino acids**.
Coordinates are 1-based and inclusive throughout.

## Local installation

Requires **Python ≥ 3.10**. Tested with Python 3.10 and Streamlit 1.56 / 1.65.

**Step 1: Clone the repository**

```bash
git clone https://github.com/<org>/<repo>.git
cd <repo>
```

**Step 2: Create and activate a Conda environment**

```bash
conda create -n poplar_viewer python=3.10 -y
conda activate poplar_viewer
pip install -r requirements.txt
```

To skip Conda, use a plain virtual environment instead (`python -m venv .venv`, activate it, then
`pip install -r requirements.txt`).

**Step 3: Run the application locally**

```bash
streamlit run app.py
```

The app opens at <http://localhost:8501>. The first load reads about 140 MB of data and takes a few
seconds; later page changes are cached.

## Deploying to Streamlit Community Cloud

1. Push the repository, including `data_for_viewer/`, to GitHub.
2. On [share.streamlit.io](https://share.streamlit.io), create an app from the repo with
   **Main file path** `app.py`, and pick **Python 3.10 or newer** under *Advanced settings*.
3. Cloud installs `requirements.txt` automatically. The app uses about 0.5 GB of memory once
   everything is loaded.

The largest data file, `ptr_all_preds.pkl`, is about 67 MB. That is under GitHub's 100 MB per-file
limit, but GitHub warns on push for files over 50 MB. Use [Git LFS](https://git-lfs.com) if you prefer.

## Project structure

```
.
├── app.py                  # Entry point: page config, global CSS/header, navigation
├── data_loader.py          # All data access: loading, AD calling, Pfam merging, orthogroups,
│                           #   CDS slicing, bulk/FASTA/full-dataset exports
├── plots.py                # Shared Plotly colors, theme handling, and figure helpers
├── views/
│   ├── nav.py              # Page registry (st.navigation)
│   ├── common.py           # Cached data access, shared UI (header, cards, badges), state helpers
│   ├── single_tf.py        # TF Viewer
│   ├── orthogroup.py       # Orthogroup Viewer
│   ├── bulk.py             # Bulk Lookup & Export
│   └── about.py            # About, citation, full-dataset downloads
├── data_for_viewer/        # Source data read at runtime (see below)
├── .streamlit/config.toml  # Theme (dark by default; switchable in the app's settings menu)
├── requirements.txt
└── _archive/               # Unreferenced files kept for provenance (git-ignored, not used by the app)
```

### Source data (`data_for_viewer/`)

These seven files are the only ones the app reads:

| File | Contents |
|---|---|
| `ptr_all_preds.pkl` | Poplar TFs: sequence, per-residue ADHunter predictions, disorder, empirical tile scores |
| `ath_all_preds.pkl` | Arabidopsis TFs: same fields; tile data from Morffy *et al.* (2024) |
| `mapping.pkl` | Poplar–Arabidopsis orthology mapping (used to build orthogroups) |
| `ad_locus_conservation_raw.pkl` | Precomputed multiple sequence alignments for 641 orthogroups |
| `pfam_hits_all_proteins.csv` | hmmscan Pfam hits (merged into domain intervals by the app) |
| `Poplar_filteredTFsClean.csv` | Poplar TF family and coding sequences |
| `AthAllTFCDS.csv` | Arabidopsis TF coding sequences |

The app also reads `Poplar_filteredTFsClean.xlsx` and `AthAllTFCDS.xlsx` in place of the CSVs if
those are present. That needs `openpyxl`, which `requirements.txt` already includes.

## Citation

> Kirkpatrick, L. D.\*, Sarkiss, A. E.\*, Waldburger, L., Kesselhaut, M. A., Zintel, M., Liu, S.,
> Jay D. Keasling, Staller, M. V., & Shih, P. M. (2026). A genome-wide map of poplar transcription
> factor activation domains uncovers pervasive turnover and tuning. \[Manuscript in preparation /
> bioRxiv forthcoming\]. DOI: \[Pending submission\]
>
> \* These authors contributed equally.

Arabidopsis empirical sort-seq tile data: Morffy, N. *et al.* Identification of plant transcriptional
activation domains. *Nature* **632**, 166–173 (2024). <https://doi.org/10.1038/s41586-024-07707-3>

## Authors & contact

- **Liam D. Kirkpatrick**: [LDKirkpatrick@berkeley.edu](mailto:LDKirkpatrick@berkeley.edu)
- **Anthony E. Sarkiss**: [Asarkiss@berkeley.edu](mailto:Asarkiss@berkeley.edu)
- **Patrick M. Shih** (Corresponding Author): [pmshih@berkeley.edu](mailto:pmshih@berkeley.edu)
