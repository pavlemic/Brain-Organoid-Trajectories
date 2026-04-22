# Transcriptional Maturation Trajectories in Brain Organoids vs. Fetal Brain

Single-cell RNA-seq comparison of brain organoids and the developing human fetal cortex. Goal: do organoids follow the same transcriptional maturation trajectories as the fetal brain, and where do they diverge?

## Datasets

| Dataset | Role | Source | Cells |
|---|---|---|---|
| **Bhaduri et al. 2020** | Brain organoids (3 protocols, GW3–24) | GEO: [GSE132672](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132672) | 241,776 |
| **Bhaduri et al. 2021** | Fetal cortex atlas (GW14–25, 11 donors) | NeMO archive (raw FASTQs on dbGaP) | 396,186 |

Same lab, same 10x Chromium v2 chemistry — minimizes batch effects between organoid and fetal sides.

## Analysis Pipeline (8 steps)

1. QC — per-dataset filtering (min/max genes, MT%)
2. Normalization + HVG selection
3. PCA → UMAP
4. Leiden clustering + cell-type annotation
5. **Integration** — Harmony on balanced 100k+100k stratified subsample
6. **Trajectory inference** — PAGA + DPT (scanpy)
7. RNA velocity — scVelo (direction of transcriptional change)
8. **Trajectory comparison** — topology, gene dynamics, organoid-specific states

Current status: Step 4 complete. Both 100k balanced subsamples produced (colab_07). About to re-run integration → annotation → trajectory on balanced input.

## Compute Split

- **Local (laptop):** code authoring, `src/` modules, small test runs.
- **Google Colab (paid, high-RAM):** heavy compute — preprocessing, integration, trajectory. GPU available.
- **Google Drive:** all `.h5ad` files. Never pushed to GitHub.
- **GitHub:** code + notebooks + notes only. Empty notebooks (no outputs) in `notebooks/colab/`; run-output copies (`_WITH_OUTPUT`) kept locally in `outputs_local/` (gitignored).

## Project Structure

```
brain-organoid-trajectories/
├── data/                          <- Not tracked (lives on Drive)
│   └── processed/                 <- .h5ad files
├── notebooks/colab/               <- Colab-facing notebooks (run on GPU/high-RAM)
│   ├── colab_00_data_download.ipynb
│   ├── colab_01_preprocessing.ipynb
│   ├── colab_02_umap_clustering.ipynb
│   ├── colab_03_integration.ipynb
│   ├── colab_04_cell_type_annotation.ipynb
│   ├── colab_05_trajectory.ipynb
│   ├── colab_06_bhaduri2021_download.ipynb
│   └── colab_07_stratified_subsample.ipynb
├── outputs_local/                 <- Run-output notebooks + plots (not tracked)
├── src/                           <- Reusable modules
│   ├── preprocessing.py           <- QC, filtering, normalization, HVG, PCA
│   ├── visualization.py           <- Plotting functions
│   └── utils.py                   <- Data loading, saving, path helpers
├── results/
│   ├── figures/
│   └── tables/
├── NOTES.md                       <- Session-by-session lab log
├── requirements.txt
└── README.md
```

## Notebook Workflow

Notebooks follow a three-phase lifecycle:

1. **Empty/authoring** — committed to GitHub in `notebooks/colab/`. Contains code cells with cell-ID headers (`### 6a — ...`) and short before-cell explanations of what each cell does. No output, no findings.
2. **Run phase** — execute on Colab with Drive mounted. Outputs are generated.
3. **Final `_WITH_OUTPUT`** — download the run notebook into `outputs_local/`. Add interpretive markdown cells after each code cell describing the actual observed results. This is the canonical record for each session.

## Reproducing

The notebooks are designed to run top-to-bottom on Google Colab with the project Drive folder mounted at `/content/drive/MyDrive/brain-organoid-trajectories/`. Local `src/` modules are reused where relevant.

```bash
# Local dev environment
python -m venv venv
source venv/bin/activate         # Windows: venv\Scripts\activate
pip install -r requirements.txt
```

## Marker Genes (Brain)

| Cell Type          | Markers                    |
|--------------------|----------------------------|
| Radial glia (vRG)  | SOX2, VIM, NES, PAX6       |
| Outer RG (oRG)     | HOPX, TNC, PTPRZ1          |
| Intermediate prog. | EOMES, PPP1R17             |
| Excitatory neuron  | NEUROD2, NEUROD6, SLC17A7, TBR1, SATB2 |
| Inhibitory neuron  | GAD1, GAD2, DLX2           |
| Astrocyte          | GFAP, AQP4, S100B          |
| Oligodendrocyte    | MBP, MOG, OLIG1/2          |
| Microglia          | AIF1, CX3CR1, P2RY12       |
| Cycling            | MKI67, TOP2A               |

## Session Log

Detailed session-by-session notes — what was run, what was observed, what broke, what we decided — live in [`NOTES.md`](NOTES.md).
