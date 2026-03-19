# Transcriptional Maturation Trajectories in Brain Organoids vs. Fetal Brain

Comparing single-cell RNA-seq data from brain organoids and fetal brain tissue
to understand transcriptional maturation trajectories.

## Project Goals

- Preprocess and explore scRNA-seq datasets
- Compare brain organoids vs. fetal brain development
- UMAP visualization and cell clustering
- Pseudotime trajectory analysis
- Cell type classification (machine learning)

## Project Structure

```
brain-organoid-trajectories/
├── data/
│   ├── raw/          <- Original, unmodified data (not tracked by Git)
│   └── processed/    <- Preprocessed AnnData objects (.h5ad)
├── notebooks/
│   ├── 01_preprocessing.ipynb   <- QC, filtering, normalization, PCA
│   └── 02_umap_clustering.ipynb <- UMAP embedding and Leiden clustering
├── src/
│   ├── preprocessing.py  <- QC, filtering, normalization, HVG, PCA
│   ├── visualization.py  <- Plotting functions
│   └── utils.py          <- Data loading, saving, path helpers
├── results/
│   ├── figures/      <- Saved plots (.png/.pdf)
│   └── tables/       <- Saved data tables (.csv)
├── requirements.txt
└── README.md
```

## Quickstart

### 1. Set up a virtual environment

```bash
# Create the environment
python -m venv venv

# Activate it
# On Windows:
venv\Scripts\activate
# On macOS/Linux:
source venv/bin/activate
```

### 2. Install dependencies

```bash
pip install -r requirements.txt
```

### 3. Register the kernel for Jupyter

```bash
python -m ipykernel install --user --name=brain-organoid --display-name "Brain Organoid"
```

### 4. Launch Jupyter and run the notebooks in order

```bash
jupyter notebook
```

Open `notebooks/01_preprocessing.ipynb` and run all cells.
Then open `notebooks/02_umap_clustering.ipynb`.

## Test Dataset

The notebooks use `pbmc3k` (built into scanpy, ~6 MB) to verify the pipeline
before loading real brain data. No manual download required.

For real brain organoid data, see:
- [Velasco et al. 2019 (Nature)](https://www.nature.com/articles/s41586-019-1289-x) — GEO: GSE132672
- [Fleck et al. 2023 (Nature)](https://www.nature.com/articles/s41586-022-05279-8) — available via the paper's data portal

## Marker Genes (Brain)

| Cell Type         | Markers                    |
|-------------------|----------------------------|
| Radial glia       | SOX2, VIM, NES, HOPX       |
| Intermediate prog.| EOMES, PPP1R17             |
| Excitatory neuron | NEUROD2, NEUROD6, SLC17A7  |
| Inhibitory neuron | GAD1, GAD2, DLX2           |
| Astrocyte         | GFAP, AQP4, S100B          |
| Oligodendrocyte   | MBP, MOG, OLIG2            |
