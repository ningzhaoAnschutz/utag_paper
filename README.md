# UTag

Code repository for: **"UTag, a de novo designed cysteine-free thermostable tagging system for tracking single mRNA translation live"**

**Authors:** Luis U. Aguilera†, Szu-Hsuan (Ashlyn) Chen†, Rhiannon M. Sears, Jake Yarbro, Jacob DeRoo, Hunter A. Ogg, Brian J. Geiss, Timothy J. Stasevich, Christopher D. Snow\*, and Ning Zhao\*

†Equal contribution. \*Corresponding authors.

---

## About

UTag is a de novo designed, cysteine-free, thermostable epitope tagging system for tracking single mRNA translation in live cells. This repository contains the Jupyter notebooks, optimization pipelines, and gene sequences needed to reproduce every figure in the manuscript, using two companion libraries: [**microlive**](https://github.com/NingZhao-Lab/microlive) for microscopy image analysis and [**tasep_models**](https://github.com/NingZhao-Lab/tasep_models) for ribosome dynamics simulations.

---

## Repository Structure

```
utag_paper/
├── notebooks/                        # Jupyter notebooks (one folder per figure)
│   ├── Figure_1/                     # UTag system characterization (1E, 1F, 1H)
│   ├── Figure_4_FRAP/                # FRAP binding affinity analysis
│   ├── Figure_5/                     # Thermostability (Western blots, CD melting)
│   ├── Figure_6/                     # Linker optimization, localization, gene maps
│   ├── Figure_7/                     # ACF, spot properties, Harringtonine chase
│   ├── Figure_S2_ITC/                # ITC binding affinity (in vitro)
│   └── Figure_Sup_SSA/               # Supplementary TASEP simulations
│
├── optimization/                     # TASEP parameter optimization (grid search)
│   ├── optimization.py               # Core χ² minimization (TASEP → ACF fitting)
│   ├── plotting_optimization.ipynb   # Cost surface visualization and CI extraction
│   ├── runner*.sh                    # Cluster submission scripts (4 datasets)
│   └── results_ACF/                  # Experimental ACF data and output CSVs
│
├── gene_sequences/                   # Plasmid DNA sequences (.dna, ~20 constructs)
│   └── utag_project/                 # UTag, SunTag, ALFAtag reporter constructs
│
├── src/                              # Shared project utilities
│   ├── imports.py                    # Shared imports, paths, and colormaps
│   └── utilities/
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
| `Figure_1E.ipynb` | Co-localization assay (nucleus/cytoplasm mEGFP ratio) |
| `Figure_1F.ipynb` | Pearson correlation analysis |
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
| `Figure_6B.ipynb` | mEGFP/mCherry intensity ratios (5/9/13 aa linkers) |
| `Figure_6D.ipynb` | Gene maps and construct visualization |
| `Figure_Translocation_Final.ipynb` | KDM5B nuclear translocation analysis and publication figures |

### Figure 7 — Translation Kinetics

| Notebook | Content |
|----------|---------|
| `ACF_calculation_data_and_controls.ipynb` | Experimental ACF calculation from tracking data (4 datasets + controls) |
| `Figure_7_AC_spot_properties.ipynb` | Spot properties: intensity, FWHM, SNR distributions |
| `Figure_7_DG_ACF.ipynb` | ACF model-vs-data comparison plots (TASEP simulation overlay) |
| `Figure_7_HK_HT_NZ.ipynb` | Harringtonine run-off assay and ke optimization |
| `calculate_ribosomal_density.ipynb` | Ribosomal density, inter-ribosome distance from optimized ki/ke |

### Supplementary Figures

| Notebook | Content |
|----------|---------|
| `Figure_S2_ITC/extracting_metadata_itc.ipynb` | Isothermal Titration Calorimetry — KD, ΔH, stoichiometry |
| `Figure_Sup_SSA/Modeling_TASEP.ipynb` | Stochastic TASEP simulations, simulated ACF, parameter recovery |

---

## Optimization Pipeline

The `optimization/` directory estimates translation kinetics by fitting a TASEP model to experimental autocorrelation functions:

1. **`optimization.py`** — 80×80 grid search over initiation rate (ki) and elongation rate (ke), minimizing a χ² cost function between simulated and experimental ACFs. Uses MD5-hashed pickle caching for simulation results.
2. **`runner*.sh`** — Cluster submission scripts for 4 datasets: UTag, UTag-C-free, SunTag, ALFAtag.
3. **`plotting_optimization.ipynb`** — Cost-surface heatmaps with Profile Likelihood confidence intervals. Reports the median of the top 5% lowest-cost grid points.
4. **`results_ACF/`** — Experimental ACF data, optimization output CSVs, and simulation cache.

Optimized parameters are independently validated by simulating Harringtonine run-off assays and comparing against experimental intensity decay curves.

---

## Environment Setup

```bash
# 1. Create the conda environment
conda env create -f environment.yml
conda activate microlive

# 2. Install companion libraries (development mode)
pip install -e /path/to/microlive
pip install -e /path/to/tasep_models
```

**Python version:** ≥ 3.10 &ensp;|&ensp; **Conda environment:** `microlive`

---

## Microscopy Data

Live-cell confocal microscopy images used in this study are publicly available on Zenodo:

> **Download:** [https://doi.org/10.5281/zenodo.19925202](https://doi.org/10.5281/zenodo.19925202)

The dataset (~13 GB uncompressed) contains 187 maximum-intensity z-projected OME-TIFF time-lapse files from U2-OS cells, organized into two directories:

| Directory | Figures | Description |
|-----------|---------|-------------|
| `Fig_7_ACF/` | 7 F–I | Extended time-lapse recordings for temporal autocorrelation analysis |
| `Fig_7_Harringtonine/` | 7 J–M | Harringtonine-treated cells for ribosome run-off kinetics |

Each directory is subdivided by tagging system (UTag, UTag_CF, SunTag, AlfaTag) and experimental date. All images were acquired on a Leica Stellaris 5 confocal microscope (63× oil objective, 512 × 512 px, 16-bit, pixel size 129.89 nm).

---

## Reproducing Figures

1. Activate the `microlive` conda environment.
2. Download microscopy data from [Zenodo](https://doi.org/10.5281/zenodo.19925202) and place it in the expected data directory.
3. Open the appropriate notebook from `notebooks/Figure_N/`.
4. Run all cells — `src/imports.py` auto-configures paths and dependencies.
5. FRAP nuclei segmentation models are auto-downloaded on first use via `microlive`.

---

## Key Dependencies

| Category | Libraries |
|----------|-----------|
| **Microscopy Engine** | [microlive](https://github.com/NingZhao-Lab/microlive) |
| **TASEP Simulations** | [tasep_models](https://github.com/NingZhao-Lab/tasep_models) |
| **Image Analysis** | scikit-image, big-fish, trackpy, tifffile, readlif |
| **Deep Learning** | PyTorch, Cellpose |
| **Scientific Computing** | NumPy, SciPy, statsmodels, Numba |
| **Visualization** | matplotlib, seaborn, joypy, brokenaxes, matplotlib-scalebar |
| **Biology** | dna-features-viewer, snapgene-reader, Biopython |

See [`environment.yml`](environment.yml) for the complete specification.

---

## License

BSD 3-Clause License. Copyright (c) 2026 Luis U. Aguilera.