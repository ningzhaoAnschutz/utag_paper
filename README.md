# UTag

Code repository for: **"UTag, a cysteine-free thermostable tagging system for tracking single mRNA translation live"**

**Authors:** Luis U Aguilera†, Szuhsuan Chen†, Gabriel Galindo, Rhiannon M. Sears, Jake Yarbro, Jacob DeRoo, Tristan Bunting, Hunter A. Ogg, Brian J. Geiss, Timothy J. Stasevich, Christopher D. Snow\*, and Ning Zhao\*

†Equal contribution. \*Corresponding authors.

---

## About

UTag is an epitope tagging system for tracking single mRNA translation in live cells, with a cysteine-free, thermostable intrabody variant. This repository contains the Jupyter notebooks, optimization pipelines, and reporter sequences used for the manuscript analyses. It uses [**microlive**](https://github.com/ningzhaoAnschutz/microlive) for microscopy image analysis and **tasep_models** for ribosome dynamics simulations. Microscopy files and some analysis inputs are stored separately.

---

## Repository Structure

```
utag_paper/
├── notebooks/                        # Jupyter notebooks grouped by analysis
│   ├── Figure_1/                     # UTag system characterization (1E, 1F, 1H)
│   ├── Figure_4_FRAP/                # FRAP binding affinity analysis
│   ├── Figure_5/                     # Thermostability (Western blots, CD melting)
│   ├── Figure_6/                     # Linker comparison and KDM5B localization
│   ├── Figure_7_ACF/                 # ACF analysis (7F-I) and ribosome density
│   ├── Figure_7_HT/                  # Harringtonine run-off (7J-M)
│   ├── Figure_7_and_S10_spot_distributions/  # 7C-E and S9B-F; legacy folder name
│   ├── Figure_8_Multiplexing/        # Paired-channel diffusion and visualization
│   ├── Figure_S2_ITC/                # ITC metadata; legacy folder name
│   ├── Figure_S7E_F_RNA_stability/   # ActD time courses and per-cell half-times
│   ├── Figure_S9A_puromycin/         # Puromycin spot runoff
│   ├── Figure_Sup_SSA/               # Supplementary TASEP simulations (S10)
│   └── raw_data_utag.xlsx            # Local journal source-data workbook
│
├── optimization/                     # TASEP parameter optimization (grid search)
│   ├── optimization.py               # Core χ² minimization (TASEP → ACF fitting)
│   ├── plotting_optimization.ipynb   # Saved cost surfaces and parameter ranges (S11)
│   ├── runner_mac.sh                 # Local runner for four datasets
│   ├── RUNNING_OPTIMIZATION_ON_MAC.md # Setup, inputs, and rerun instructions
│   ├── source_files_ACF/data_ACF/    # Experimental ACF inputs and count reports
│   ├── results_ACF/                  # Saved optimization grids and parameter CSVs
│   └── optimization_parameters/     # Generated cost-surface panels (local)
│
├── gene_sequences/
│   └── utag_project/                 # Reporter inputs: pNZ208 (UTag), pNZ266 (SunTag), pNZ267 (ALFA)
├── kk/                               # Local archives, including kk_old_gene_sequences/
├── export_files/                     # Local LIF-to-TIFF/AVI exporter (not tracked by Git)
│
├── src/                              # Shared project utilities
│   ├── imports.py                    # Shared imports, paths, and colormaps
│   └── utilities/
│       ├── inhibitors.py             # HT tracking-data loading; analysis for S7E-F and S9A
│       └── extract_laser_intensities.py  # LIF laser metadata extraction
├── environment.yml                   # Conda environment specification
├── LICENSE                           # BSD 3-Clause
└── README.md
```

---

## Notebooks

### Figure 1 — UTag System Characterization

| Notebook | Content |
|----------|---------|
| `Figure_1E.ipynb` | Pearson correlation analysis |
| `Figure_1F.ipynb` | Co-localization assay (nucleus/cytoplasm mEGFP ratio) |
| `Figure_1H.ipynb` | UTag characterization panel H |

### Figure 4 — FRAP Binding Affinity

| Notebook | Content |
|----------|---------|
| `FRAP_processing.ipynb` | Raw FRAP data processing and curve extraction |
| `Figure_4_FRAP_interpretation.ipynb` | Recovery curves, t½, endpoint comparison across tagging systems |
| `FRAP_representative_images.ipynb` | Representative FRAP image montages |

### Figure 5 — Thermostability

| Notebook | Content |
|----------|---------|
| `Figure_5_C.ipynb` | Western blot quantification (4–70 °C) |
| `Figure_5_D_E_TM.ipynb` | CD thermal denaturation melting curves (Tm) |

### Figure 6 — Linker Optimization & Localization

| Notebook | Content |
|----------|---------|
| `Figure_6C.ipynb` | mEGFP/mCherry intensity ratios (5/9/13 aa linkers); some output names still use 6B |
| `Figure_6D.ipynb` | Representative UTag and cysteine-free UTag KDM5B localization images |
| `Figure_Translocation.ipynb` | KDM5B localization image processing and plotting |

### Figure 7 — Translation Kinetics

| Notebook | Content |
|----------|---------|
| `Figure_7_ACF/ACF_calculation_data_and_controls.ipynb` | Experimental ACF calculation from tracking data (4 datasets + controls) |
| `Figure_7_and_S10_spot_distributions/Fig_S10_spot_distributions_removing_aggreates.ipynb` | Aggregate-filtered spot properties for 7C-E and S9B-F; file and folder names retain older figure numbers |
| `Figure_7_ACF/Figure_7_F_G_H_I__ACF.ipynb` | ACF model-vs-data comparison plots (TASEP simulation overlay) |
| `Figure_7_HT/Figure_7_J_K_L_M__HT_linear_fit.ipynb` | Published harringtonine run-off linear fits (7J-M) |
| `Figure_7_HT/Figure_7_J_K_L_M__HT.ipynb` and `Figure_7_J_K_L_M__HT_variable_ki.ipynb` | Alternative simulation-based HT fits |
| `Figure_7_ACF/calculate_ribosomal_density.ipynb` | Ribosomal density and inter-ribosome distance from saved initiation and elongation rates |

### Figure 8 — Multiplexing

| Notebook in `Figure_8_Multiplexing/` | Content |
|----------|---------|
| `Figure_8_multiplexing_diffusion_with_without_aggregates.ipynb` | Paired-channel diffusion and spot-property analysis with and without aggregate filtering |
| `Figure_8_multiplexing_diffusion_filtered.ipynb` | Single-configuration diffusion analysis |
| `Figure_8_multiplexing_visualization_v5.ipynb` | Representative cell images, trajectories, and spot montages |

The diffusion notebooks use `multiplexing_msd.py` in the same folder.

### Supplementary Figures

| Notebook | Content |
|----------|---------|
| `Figure_S2_ITC/extracting_metadata_itc.ipynb` | Isothermal titration calorimetry metadata; S2 is a legacy folder name |
| `Figure_S7E_F_RNA_stability/Figure_S7E_F_RNA_stability.ipynb` | ActD time courses (S7E) and per-cell mRNA half-times (S7F) |
| `Figure_S9A_puromycin/Figure_S9A_puromycin.ipynb` | Translation spot counts after puromycin treatment (S9A) |
| `Figure_Sup_SSA/Modeling_TASEP.ipynb` | Stochastic TASEP simulations, simulated ACFs, and parameter sweeps |

---

## Optimization Pipeline

The `optimization/` directory estimates translation kinetics by fitting a TASEP model to experimental autocorrelation functions:

1. **`optimization.py`** — 80×80 grid search over initiation rate (ki) and elongation rate (ke), minimizing a χ² cost function between simulated and experimental ACFs. Uses MD5-hashed pickle caching for simulation results.
2. **`runner_mac.sh`** — Local runner for four datasets: UTag, UTag-C-free, SunTag, ALFAtag.
3. **`plotting_optimization.ipynb`** — Redraws Figure S11 from saved grids. Reports medians of the 320 lowest-cost grid points (top 5%) and their 5th-95th percentile ranges, not calibrated confidence intervals.
4. **`source_files_ACF/data_ACF/`** — Experimental ACF inputs and retained cell/trajectory counts.
5. **`results_ACF/`** — Saved grid-search CSVs and full-precision parameter summaries. Simulation caches are under `source_files_ACF/simulation_cache/`; caches, logs, and generated plots are excluded from Git.

The HT simulation notebooks compare model predictions with experimental intensity decay. The published 7J-M panels use the separate linear-fit notebook. See [RUNNING_OPTIMIZATION_ON_MAC.md](optimization/RUNNING_OPTIMIZATION_ON_MAC.md) before a new fit: stochastic results can differ, and reruns overwrite that dataset's result tables. Use a separate checkout to preserve the published run.

---

## Environment Setup

```bash
# 1. Create the conda environment
conda env create -f environment.yml
conda activate microlive

# 2. Optional: use local development checkouts instead of the installed packages
pip install -e /path/to/microlive
pip install -e /path/to/tasep_models
```

**Python version:** ≥ 3.10 &ensp;|&ensp; **Conda environment:** `microlive`

---

## Microscopy Data

The earlier Figure 7 live-cell microscopy release is publicly available on Zenodo:

> **Download:** [https://doi.org/10.5281/zenodo.19925202](https://doi.org/10.5281/zenodo.19925202)

That release (~13 GB uncompressed) contains 187 maximum-intensity z-projected OME-TIFF time-lapse files from U-2 OS cells, organized into two directories:

| Directory | Figures | Description |
|-----------|---------|-------------|
| `Fig_7_ACF/` | 7 F–I | Extended time-lapse recordings for temporal autocorrelation analysis |
| `Fig_7_Harringtonine/` | 7 J–M | Harringtonine-treated cells for ribosome run-off kinetics |

Each directory is subdivided by tagging system (UTag, UTag_CF, SunTag, AlfaTag) and experimental date. All images were acquired on a Leica Stellaris 5 confocal microscope (63× oil objective, 512 × 512 px, 16-bit, pixel size 129.89 nm).

The expanded local package is in `export_files/outputs/final/`: 274 OME-TIFF time series and 274 display AVIs across Figure 7 ACF, Figure 7 HT, Figure 8 multiplexing, S7E-F RNA stability, and S9A puromycin. Figure 8 TIFFs and AVIs contain all 60 frames used for MSD analysis. The other collections retain their selected analysis windows. This package is not yet described by the older Zenodo release linked above.

TIFFs retain uncorrected source intensities, maximum-projected over Z where applicable. AVIs use a fixed 0.4-pixel spatial Gaussian and fixed movie/channel contrast limits for display only; do not use them for intensity quantification. `files.csv` records source scenes, frame ranges, and checksums; `analysis_links.csv` links scenes to analysis inputs. See `export_files/README.md` for export details. Export tools, microscopy files, tracking CSVs, the source-data workbook, and local archives are excluded from Git.

---

## Reproducing Figures

The Figure 7 HT notebooks and S9A use the same tracking CSV loader in `src/utilities/inhibitors.py`. The published 7J-M linear fits are computed in `Figure_7_J_K_L_M__HT_linear_fit.ipynb`, not by the utility's generic fit routine. S7E-F and S9A use the utility for analysis and plotting; each still reads its own local inputs and writes to its own `figures/` and `raw_data/`. Run these notebooks within this repository using the `microlive` kernel.

1. Activate the `microlive` conda environment.
2. Obtain the required inputs and check the notebook's local paths. The linked [Zenodo release](https://doi.org/10.5281/zenodo.19925202) contains Figure 7 images, not every input for every notebook.
3. Open the appropriate notebook from the folders listed above.
4. Check output paths before running cells. Many notebooks use `src/imports.py` for shared imports and repository paths; this does not download missing data or dependencies.
5. The published HT outputs are in `figures_linear/` and `raw_data_linear/`. S7E-F and S9A use their own `figures/` and `raw_data/` folders. The local journal source-data workbook is `notebooks/raw_data_utag.xlsx`.

---

## Key Dependencies

| Category | Libraries |
|----------|-----------|
| **Microscopy Engine** | [microlive](https://github.com/ningzhaoAnschutz/microlive) |
| **TASEP Simulations** | `tasep_models` (installed as `tasep-models`) |
| **Image Analysis** | scikit-image, big-fish, trackpy, tifffile, readlif |
| **Deep Learning** | PyTorch, Cellpose |
| **Scientific Computing** | NumPy, SciPy, statsmodels, Numba |
| **Visualization** | matplotlib, seaborn, joypy, brokenaxes, matplotlib-scalebar |
| **Biology** | dna-features-viewer, snapgene-reader, Biopython |

See [`environment.yml`](environment.yml) for the complete specification.

---

## License

Code: [BSD 3-Clause License](LICENSE), copyright (c) 2025 Ning Zhao. See the Zenodo record for the microscopy data license.
