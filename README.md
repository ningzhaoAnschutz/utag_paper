# UTag

Analysis code for **"UTag, a cysteine-free thermostable tagging system for tracking single mRNA translation live"**.

**Authors:** Luis U Aguilera†, Szuhsuan Chen†, Gabriel Galindo, Rhiannon M. Sears, Jake Yarbro, Jacob DeRoo, Tristan Bunting, Hunter A. Ogg, Brian J. Geiss, Timothy J. Stasevich, Christopher D. Snow\*, and Ning Zhao\*

†Equal contribution. \*Corresponding authors.

## Overview

UTag is an epitope tagging system for imaging the translation of single mRNAs in live cells. It is paired with a cysteine-free, thermostable intrabody. This repository contains the notebooks used to produce the manuscript figures, the TASEP parameter-fitting code, and the reporter plasmid maps.

| Resource             | Link                                                       |
| -------------------- | ---------------------------------------------------------- |
| Microscopy data      | [Zenodo](https://doi.org/10.5281/zenodo.24091136)          |
| Image analysis       | [microlive](https://github.com/ningzhaoAnschutz/microlive) |
| Ribosome simulations | [tasep_models](https://pypi.org/project/tasep-models/)     |

## Installation

```bash
conda env create -f environment.yml
conda activate microlive
```

This creates an environment named `microlive` (Python ≥ 3.10) that includes `microlive` and `tasep_models`. To develop against local copies of those packages instead, run:

```bash
pip install -e /path/to/microlive
pip install -e /path/to/tasep_models
```

## Repository structure

```
utag_paper/
├── notebooks/                 # One folder per figure (see Notebooks below)
├── optimization/              # TASEP parameter fitting (Figure 7, Figure S11, Table S7)
├── gene_sequences/
│   └── utag_project/          # Plasmid maps: pNZ208 (UTag), pNZ266 (SunTag), pNZ267 (ALFA-tag)
├── src/
│   ├── imports.py             # Shared imports, paths, and colormaps
│   └── utilities/
│       ├── inhibitors.py                  # Loading and analysis for HT, ActD, and puromycin data
│       └── extract_laser_intensities.py   # Reads laser settings from LIF metadata
├── environment.yml
├── LICENSE
└── README.md
```

## Notebooks

All notebooks are in `notebooks/`. Each table lists the figure panels a notebook produces and what it does. A dash (—) marks preprocessing notebooks that feed other notebooks.

### Figure 1: Intrabody localization

Folder: `notebooks/Figure_1/`

| Figure | Notebook          | Description                                                                          |
| ------ | ----------------- | ------------------------------------------------------------------------------------ |
| 1E     | `Figure_1E.ipynb` | Line profiles and Pearson's r for wtscFv, Anti-UTag-IB(15F11), and Anti-UTag-IB(2E2) |
| 1F     | `Figure_1F.ipynb` | Nucleus-to-cytoplasm intensity ratio for the same three intrabodies                  |
| 1H     | `Figure_1H.ipynb` | Nucleus-to-cytoplasm ratio of Anti-UTag-IB-mEGFP with mCh-H2B or 4xUTag-mCh-H2B      |

### Figure 4: FRAP

Folder: `notebooks/Figure_4_FRAP/`

| Figure | Notebook                             | Description                                              |
| ------ | ------------------------------------ | -------------------------------------------------------- |
| —      | `FRAP_processing.ipynb`              | Extracts FRAP recovery curves from the raw recordings    |
| 4A–H   | `Figure_4_FRAP_interpretation.ipynb` | Recovery curves, t½, and endpoint comparison across tags |
| 4      | `FRAP_representative_images.ipynb`   | Representative image montages and kymographs             |

Shared plotting functions are in `frap_plotting_module.py`. The exported source-data tables are described in `raw_data/README_source_data.md`.

### Figure 5: Thermostability

Folder: `notebooks/Figure_5/`

| Figure | Notebook                | Description                                 |
| ------ | ----------------------- | ------------------------------------------- |
| 5C     | `Figure_5_C.ipynb`      | Western blot quantification from 4 to 70 °C |
| 5D–E   | `Figure_5_D_E_TM.ipynb` | Circular dichroism melting curves and Tm    |

The exported source-data tables are described in `raw_data/README_source_data.md`.

### Figure 6: Linkers and localization

Folder: `notebooks/Figure_6/`

| Figure | Notebook                     | Description                                                       |
| ------ | ---------------------------- | ----------------------------------------------------------------- |
| 6C     | `Figure_6C.ipynb`            | mEGFP/mCherry intensity ratios for 5, 9, and 13 aa linkers        |
| —      | `Figure_Translocation.ipynb` | Processes KDM5B localization images and crops regions of interest |
| 6D     | `Figure_6D.ipynb`            | Representative KDM5B localization with UTag and UTag(ΔCys)        |

The exported source-data tables are described in `raw_data/README_source_data.md`.

### Figure 7: Translation kinetics

Figure 7 notebooks are split across three folders.

**Spot properties:** `notebooks/Figure_7_and_S9_spot_distributions/`

| Figure      | Notebook                                             | Description                                                                          |
| ----------- | ---------------------------------------------------- | ------------------------------------------------------------------------------------ |
| 7C–E, S9B–D | `Fig_S9_spot_distributions.ipynb`                    | Spot properties before and after removing spots larger than 0.5 µm                   |
| S9E–F       | `Fig_S9_spot_distributions_removing_aggreates.ipynb` | Trajectories per cell before and after removing trajectories with mean size > 0.5 µm |

Shared functions are in `processing_spot_size.py`. `simulated_spots_size_snr.py` is a standalone script that illustrates how spot size and SNR are measured, using two synthetic spots.

**Autocorrelation (ACF):** `notebooks/Figure_7_ACF/`

| Figure   | Notebook                                  | Description                                                     |
| -------- | ----------------------------------------- | --------------------------------------------------------------- |
| —        | `ACF_calculation_data_and_controls.ipynb` | Computes experimental ACFs for the four tags and their controls |
| 7F–I     | `Figure_7_F_G_H_I__ACF.ipynb`             | Experimental ACFs overlaid with fitted TASEP simulations        |
| Table S7 | `calculate_ribosomal_density.ipynb`       | Ribosome density and spacing from the fitted rates              |

**Harringtonine run-off:** `notebooks/Figure_7_HT/`

| Figure | Notebook                                | Description                                                        |
| ------ | --------------------------------------- | ------------------------------------------------------------------ |
| 7J–M   | `Figure_7_J_K_L_M__HT_linear_fit.ipynb` | Linear fits to run-off curves, run-off times, and elongation rates |

### Figure 8: Multiplexing

Folder: `notebooks/Figure_8_Multiplexing/`

| Figure | Notebook                                                        | Description                                                                         |
| ------ | --------------------------------------------------------------- | ----------------------------------------------------------------------------------- |
| 8      | `Figure_8_multiplexing_diffusion_filtered.ipynb`                | Diffusion of paired SunTag and UTag(ΔCys) spots, with a selectable spot-size filter |
| 8      | `Figure_8_multiplexing_diffusion_with_without_aggregates.ipynb` | Same analysis with and without removing large (aggregate) spots                     |
| 8      | `Figure_8_multiplexing_visualization_v5.ipynb`                  | Representative cell images, trajectories, and spot montages                         |

Shared diffusion functions are in `multiplexing_msd.py`.

### Supplementary figures

| Figure | Notebook                                                      | Description                                                    |
| ------ | ------------------------------------------------------------- | -------------------------------------------------------------- |
| S2     | `Figure_S2_ITC/extracting_metadata_itc.ipynb`                 | Extracts metadata from isothermal titration calorimetry files  |
| S7E–F  | `Figure_S7E_F_RNA_stability/Figure_S7E_F_RNA_stability.ipynb` | mRNA decay after actinomycin D and per-cell half-lives         |
| S9A    | `Figure_S9A_puromycin/Figure_S9A_puromycin.ipynb`             | Translation spot counts after puromycin treatment              |
| S10    | `Figure_Sup_SSA/Modeling_TASEP.ipynb`                         | Stochastic TASEP simulations, simulated ACFs, parameter sweeps |

Paths are relative to `notebooks/`. Figure S11 is produced by `optimization/plotting_optimization.ipynb` (see [TASEP parameter fitting](#tasep-parameter-fitting)).

## TASEP parameter fitting

The `optimization/` folder estimates the translation initiation rate (ki) and elongation rate (ke) for each tag by fitting a TASEP model to the experimental ACFs.

| File or folder                   | Purpose                                                                       |
| -------------------------------- | ----------------------------------------------------------------------------- |
| `optimization.py`                | 80 × 80 grid search over ki and ke that minimizes a χ² cost                   |
| `runner_mac.sh`                  | Runs the fit for one tag: UTag, UTag(ΔCys), SunTag, or ALFA-tag               |
| `plotting_optimization.ipynb`    | Plots the cost surfaces and parameter ranges (Figure S11)                     |
| `source_files_ACF/data_ACF/`     | Experimental ACFs and the number of cells and trajectories per tag            |
| `results_ACF/`                   | Published grids (`results_df_*.csv`) and estimates (`optimized_params_*.csv`) |
| `RUNNING_OPTIMIZATION_ON_MAC.md` | Step-by-step instructions for running a new fit                               |

Reported rates are the median of the 320 lowest-cost grid points (top 5%), with the 5th–95th percentile range of those points. These ranges describe the spread of good fits and are not confidence intervals.

> **Note:** Fits are stochastic. Running a new fit overwrites that tag's files in `results_ACF/`, so work in a separate copy of the repository if you want to keep the published results.

## Microscopy data

The live-cell imaging data are available on [Zenodo](https://doi.org/10.5281/zenodo.24091136) as `utag_data_v2.zip` (about 6.3 GB). The archive contains 273 time-lapse recordings of U-2 OS cells, each provided as an OME-TIFF for analysis and an AVI movie for viewing.

| Folder                        | Figure | Recordings | Frames     | Interval |
| ----------------------------- | ------ | ---------: | ---------- | -------- |
| `Figure_7_ACF/`               | 7F–I   |         60 | 360 or 365 | 5 s      |
| `Figure_7_HT/`                | 7J–M   |        126 | 30         | 60 s     |
| `Figure_8_Multiplexing/`      | 8      |         23 | 60         | 1 s      |
| `Figure_S7E_F_RNA_Stability/` | S7E–F  |         50 | 13         | 20 min   |
| `Figure_S9A_Puromycin/`       | S9A    |         14 | 35         | 20 s     |

**Acquisition:** Leica Stellaris 5 confocal microscope, 63× oil objective, 512 × 512 pixels, 16-bit, about 130 nm per pixel (exact values are in each file's OME metadata). Z-stacks are shown as maximum-intensity projections; Figure 8 recordings are single planes. Frame numbers are zero-based.

**Processing:** The OME-TIFFs are unfiltered and have no photobleaching correction. The AVIs are lightly smoothed and contrast-adjusted for viewing only, so use the OME-TIFFs for any quantification.

**Dataset notes:**

- **Harringtonine (7J–M):** Frame 5 is the treatment reference point. Of the 126 cells, 117 responded and were used for the fits; the 9 non-responders are included to show how responders were selected.
- **RNA stability (S7E–F):** 52 cells analyzed. DMSO controls appear in S7E only.
- **Puromycin (S9A):** Frame 4 is the treatment reference point. 15 cells analyzed.
- **Multiplexing (8):** The SunTag and UTag(ΔCys) channels are paired measurements of the same cells.

**Index files in the archive:**

- `files.csv`: source recording, frame range, and SHA-256 checksum for every file.
- `analysis_links.csv`: maps each analyzed cell to its recording. There are 276 analyzed cells from 273 recordings because three recordings contain two analyzed cells.
- `README.md`: instructions for opening the files in Fiji.

## Reproducing the figures

1. Install and activate the `microlive` environment, then select it as the Jupyter kernel.
2. Obtain the input data for the figure you want to reproduce (microscopy recordings are on [Zenodo](https://doi.org/10.5281/zenodo.24091136)). The ACF inputs for parameter fitting are already included in `optimization/source_files_ACF/data_ACF/`.
3. Open the notebook and set the input and output paths in its first cells.
4. Run the notebook from its own folder, top to bottom. Run interactive cropping cells one at a time.

Most notebooks save figures and source-data tables to subfolders inside their own folder, such as `figures/`, `plots/`, `raw_data/`, or `results_*/`.

## Dependencies

| Purpose                 | Packages                                                    |
| ----------------------- | ----------------------------------------------------------- |
| Microscopy analysis     | microlive                                                   |
| TASEP simulations       | tasep_models                                                |
| Image analysis          | scikit-image, big-fish, trackpy, tifffile, readlif          |
| Segmentation            | PyTorch, Cellpose                                           |
| Numerics and statistics | NumPy, SciPy, statsmodels, Numba                            |
| Plotting                | matplotlib, seaborn, joypy, brokenaxes, matplotlib-scalebar |
| Sequences               | dna-features-viewer, snapgene-reader, Biopython             |

The full list is in [environment.yml](environment.yml).

## License

The code is released under the [BSD 3-Clause License](LICENSE), copyright (c) 2025 Ning Zhao. The license for the microscopy data is listed on the Zenodo record.
