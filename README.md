# UTag

Code repository for **"UTag, a cysteine-free thermostable tagging system for tracking single mRNA translation live"**.

**Authors:** Luis U Aguilera†, Szuhsuan Chen†, Gabriel Galindo, Rhiannon M. Sears, Jake Yarbro, Jacob DeRoo, Tristan Bunting, Hunter A. Ogg, Brian J. Geiss, Timothy J. Stasevich, Christopher D. Snow\*, and Ning Zhao\*

†Equal contribution. \*Corresponding authors.

## About

UTag is an epitope tagging system for tracking single mRNA translation in live cells. It includes a cysteine-free, thermostable intrabody variant. This repository contains the notebooks, TASEP parameter optimization code, and reporter plasmid maps used for the manuscript analyses. Image analysis uses [microlive](https://github.com/ningzhaoAnschutz/microlive), and ribosome simulations use [tasep_models](https://pypi.org/project/tasep-models/). Image data are on [Zenodo](https://doi.org/10.5281/zenodo.24091136). Tracking CSVs and most other notebook inputs are not in this repository (see [Data not in this repository](#data-not-in-this-repository)).

## Repository structure

```
utag_paper/
├── notebooks/
│   ├── Figure_1/                     # 1E, 1F, 1H
│   ├── Figure_4_FRAP/                # FRAP analysis and frap_plotting_module.py
│   ├── Figure_5/                     # Western blot and CD melting curves
│   ├── Figure_6/                     # Linker comparison and KDM5B localization
│   ├── Figure_7_ACF/                 # ACF analysis (7F-I) and ribosome density
│   ├── Figure_7_HT/                  # Harringtonine run-off linear fits (7J-M)
│   ├── Figure_7_and_S9_spot_distributions/   # 7C-E and S9B-F
│   ├── Figure_8_Multiplexing/        # Paired-channel diffusion and images
│   ├── Figure_S2_ITC/                # ITC metadata (older folder name)
│   ├── Figure_S7E_F_RNA_stability/   # ActD time courses and per-cell half-times
│   ├── Figure_S9A_puromycin/         # Puromycin spot run-off; includes figures/FigS9A.tiff
│   └── Figure_Sup_SSA/               # TASEP simulations (S10)
├── optimization/
│   ├── optimization.py               # Grid search fitting TASEP ACFs to data
│   ├── runner_mac.sh                 # macOS runner for the four datasets
│   ├── plotting_optimization.ipynb   # Cost surfaces and parameter ranges (S11)
│   ├── source_files_ACF/data_ACF/    # Experimental ACFs and cell/trajectory counts
│   └── results_ACF/                  # Saved grids and parameter estimates
├── gene_sequences/utag_project/      # pNZ208 (UTag), pNZ266 (SunTag), pNZ267 (ALFA-tag)
├── src/
│   ├── imports.py                    # Shared imports, paths, and colormaps
│   └── utilities/
│       ├── inhibitors.py             # Tracking-data loading; HT, ActD, and puromycin analysis
│       └── extract_laser_intensities.py  # Laser settings from LIF metadata
├── environment.yml
├── LICENSE
└── README.md
```

The Figure 4, 5, and 6 folders also contain `raw_data/README_source_data.md`, which describes the exported source-data tables.

## Notebooks

### Figure 1

| Notebook            | Content                                                                                                |
| ------------------- | ------------------------------------------------------------------------------------------------------ |
| `Figure_1E.ipynb` | Line profiles and Pearson r for wtscFv, Anti-UTag-IB(15F11), and Anti-UTag-IB(2E2) with 4xUTag-mCh-H2B |
| `Figure_1F.ipynb` | Nucleus/cytoplasm intensity ratio for the same three intrabodies                                       |
| `Figure_1H.ipynb` | Nucleus/cytoplasm intensity ratio of Anti-UTag-IB-mEGFP with mCh-H2B or 4xUTag-mCh-H2B                 |

### Figure 4: FRAP

| Notebook                               | Content                                                              |
| -------------------------------------- | -------------------------------------------------------------------- |
| `FRAP_processing.ipynb`              | Processes raw FRAP data and extracts recovery curves                 |
| `Figure_4_FRAP_interpretation.ipynb` | Recovery curves, t½, and endpoint comparison across tagging systems |
| `FRAP_representative_images.ipynb`   | Representative FRAP image montages and kymographs                    |

Both plotting notebooks use `frap_plotting_module.py` in the same folder.

The raw LIF input folder is `/Volumes/Luis_DRIVE/UTag_paper_data/FRAP/`.

### Figure 5: Thermostability

| Notebook                  | Content                                 |
| ------------------------- | --------------------------------------- |
| `Figure_5_C.ipynb`      | Western blot quantification (4–70 °C) |
| `Figure_5_D_E_TM.ipynb` | CD thermal denaturation curves and Tm   |

### Figure 6: Linkers and localization

| Notebook                       | Content                                                                                    |
| ------------------------------ | ------------------------------------------------------------------------------------------ |
| `Figure_6C.ipynb`            | mEGFP/mCherry intensity ratios for 5, 9, and 13 aa linkers; some output names still use 6B |
| `Figure_6D.ipynb`            | Representative UTag and cysteine-free UTag KDM5B localization images                       |
| `Figure_Translocation.ipynb` | KDM5B localization image processing and interactive cropping                               |

The localization notebooks read LIF files from `/Volumes/Luis_DRIVE/UTag_paper_data/Translocation/`.

### Figure 7: Translation kinetics

| Notebook                                                                                  | Content                                                                                                      |
| ----------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------ |
| `Figure_7_and_S9_spot_distributions/Fig_S9_spot_distributions.ipynb`                    | Spot properties before and after excluding detections larger than 0.5 µm (7C-E and S9B-D)                   |
| `Figure_7_and_S9_spot_distributions/Fig_S9_spot_distributions_removing_aggreates.ipynb` | Total and retained trajectories per cell after excluding whole trajectories with mean size > 0.5 µm (S9E-F) |
| `Figure_7_ACF/ACF_calculation_data_and_controls.ipynb`                                  | Experimental ACFs from tracking data (four datasets and controls)                                            |
| `Figure_7_ACF/Figure_7_F_G_H_I__ACF.ipynb`                                              | Experimental ACFs overlaid with TASEP simulations                                                            |
| `Figure_7_ACF/calculate_ribosomal_density.ipynb`                                        | Ribosome density and inter-ribosome distance from the fitted initiation and elongation rates                 |
| `Figure_7_HT/Figure_7_J_K_L_M__HT_linear_fit.ipynb`                                     | Harringtonine run-off: linear fits, run-off times, and elongation rates                                      |

The spot-distribution notebooks use `processing_spot_size.py` in the same folder. `simulated_spots_size_snr.py` is a standalone script that measures two synthetic Gaussian spots with microlive to illustrate how spot size and SNR are defined.

### Figure 8: Multiplexing

| Notebook in `Figure_8_Multiplexing/`                            | Content                                                                            |
| ----------------------------------------------------------------- | ---------------------------------------------------------------------------------- |
| `Figure_8_multiplexing_diffusion_with_without_aggregates.ipynb` | Paired-channel diffusion and spot properties, with and without aggregate filtering |
| `Figure_8_multiplexing_diffusion_filtered.ipynb`                | Diffusion analysis for a single filter setting                                     |
| `Figure_8_multiplexing_visualization_v5.ipynb`                  | Representative cell images, trajectories, and spot montages                        |

The diffusion notebooks use `multiplexing_msd.py` in the same folder.

### Supplementary figures

| Notebook                                                        | Content                                                                  |
| --------------------------------------------------------------- | ------------------------------------------------------------------------ |
| `Figure_S2_ITC/extracting_metadata_itc.ipynb`                 | Isothermal titration calorimetry metadata                                |
| `Figure_S7E_F_RNA_stability/Figure_S7E_F_RNA_stability.ipynb` | ActD time courses (S7E) and per-cell mRNA half-times (S7F)               |
| `Figure_S9A_puromycin/Figure_S9A_puromycin.ipynb`             | Translation spot counts after puromycin treatment (S9A)                  |
| `Figure_Sup_SSA/Modeling_TASEP.ipynb`                         | Stochastic TASEP simulations, simulated ACFs, and parameter sweeps (S10) |

## TASEP parameter optimization

`optimization/` estimates initiation (ki) and elongation (ke) rates by fitting a TASEP model to the experimental ACFs.

- `optimization.py` runs an 80 × 80 grid search over ki and ke and minimizes a χ² cost between simulated and experimental ACFs. Simulation results are cached as pickles keyed by an MD5 hash of the inputs.
- `runner_mac.sh` runs one dataset at a time: UTag, UTag C-free, SunTag, or ALFA-tag.
- `plotting_optimization.ipynb` redraws Figure S11 from the saved grids. It reports the median of the 320 lowest-cost grid points (top 5%) and their 5th–95th percentile range. These ranges are not calibrated confidence intervals.
- `source_files_ACF/data_ACF/` holds the experimental ACF inputs and the number of cells and trajectories retained for each dataset.
- `results_ACF/` holds the saved grids (`results_df_*.csv`) and parameter estimates (`optimized_params_*.csv`) used for Figure 7 and Table S7.

Simulation caches, logs, and generated plots (`optimization_parameters/`) are not tracked. Read [RUNNING_OPTIMIZATION_ON_MAC.md](optimization/RUNNING_OPTIMIZATION_ON_MAC.md) before starting a new fit. Results are stochastic, and a rerun overwrites that dataset's tables in `results_ACF/`, so use a separate checkout if you need to keep the published run.

## Environment

```bash
conda env create -f environment.yml
conda activate microlive

# Optional: use local checkouts instead of the installed packages
pip install -e /path/to/microlive
pip install -e /path/to/tasep_models
```

The environment is named `microlive` and requires Python ≥ 3.10.

## Microscopy data

Live-cell image data are on [Zenodo](https://doi.org/10.5281/zenodo.24091136).

The deposit (`utag_data_v2.zip`, approximately 6.3 GB) contains 273 OME-TIFF time series and 273 matching AVI movies from U-2 OS cells, exported from selected recordings in 37 Leica LIF files. The LIF files themselves are not included.

| Folder                          | Time-lapse recordings | Figures | Frames and interval                                                   |
| ------------------------------- | --------------------: | ------- | --------------------------------------------------------------------- |
| `Figure_7_ACF/`               |                    60 | 7F–I   | 360 frames for 48 recordings; 365 frames for 12 recordings; 5 s apart |
| `Figure_7_HT/`                |                   126 | 7J–M   | 30 frames, 60 s; treatment reference at frame 5                       |
| `Figure_8_Multiplexing/`      |                    23 | 8       | 60 frames, 1 s; SunTag and UTag(ΔCys) channels                       |
| `Figure_S7E_F_RNA_Stability/` |                    50 | S7E–F  | 13 time points, 20 min apart, spanning 0-4 h                          |
| `Figure_S9A_Puromycin/`       |                    14 | S9A     | 35 frames, 20 s; treatment reference at frame 4                       |

The HT collection includes 117 responding cells used for the curves and fits and nine nonresponders retained to document responder selection. The RNA-stability recordings contain 52 analyzed cells; DMSO controls contribute to S7E only. The puromycin recordings contain 15 analyzed cells.

Frame numbers are zero-based. Images were acquired on a Leica Stellaris 5 confocal microscope with a 63× oil objective (512 × 512 px, 16-bit). Pixel size is about 130 nm; the exact value for each file is stored in its OME metadata. Z-stacks were maximum-intensity projected; Figure 8 retains its native single plane. No smoothing or photobleaching correction was applied to the TIFFs. The AVIs use a 0.4-pixel spatial Gaussian and fixed contrast limits per recording and channel for viewing only, so use the TIFFs for quantification.

`files.csv` gives the source LIF, original scene identifier, frame range, and SHA-256 checksum for each file. `analysis_links.csv` maps each analysis input to its recording; 276 analysis inputs come from 273 recordings because three recordings contain two analyzed cells. Figure 8's two channels are paired measurements of the same cells. The deposit's `README.md` explains how to open the files in Fiji.

An [earlier release](https://doi.org/10.5281/zenodo.19925202) containing only Figure 7 recordings is superseded by this deposit.

## Data not in this repository

The following are kept locally and excluded by `.gitignore`:

- Source LIF files, tracking CSVs, and other notebook inputs.
- `notebooks/raw_data_utag.xlsx`, the journal source-data workbook.
- `export_files/`, the code that produced the Zenodo TIFFs and AVIs.
- Generated figures and numerical source-data exports. The source-data README files and export plans are tracked. Other tracked exceptions are `Figure_S9A_puromycin/figures/FigS9A.tiff` and the optimization files listed above.

## Reproducing figures

1. Activate the `microlive` environment and select it as the notebook kernel.
2. Obtain the inputs for the notebook. The image deposit does not include the tracking CSVs, quantification workbooks, or FRAP analysis tables read by the plotting notebooks. The experimental ACF inputs needed for parameter optimization are included under `optimization/source_files_ACF/data_ACF/`.
3. Check the input and output paths in the first cells. Many notebooks find the repository root through `src/imports.py`, which does not download data.
4. Run the notebook from its folder. Run interactive cropping cells individually and check the configured output paths before saving.

FRAP and localization LIF inputs use the Luis_DRIVE folders listed above. ACF and spot-property notebooks read `/Volumes/Luis_DRIVE/UTag_paper_data/ACF/`; HT reads `/Volumes/Luis_DRIVE/UTag_paper_data/Harringtonine/`. Figure 5, Figure 6C, and ITC metadata notebooks retain OneDrive paths. RNA-stability, puromycin, and multiplexing inputs are in notebook-local data folders.

Output folders differ by analysis:

Below, `Figure_*` folders are under `notebooks/`; `optimization/` is at the repository root.

| Analysis                               | Outputs                                                                                                                       |
| -------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------- |
| Figure 1                               | `Figure_1/plots/`, with a subfolder for each panel                                                                          |
| FRAP processing                        | `Figure_4_FRAP/results/results_FRAP_<LIF name>/`                                                                            |
| FRAP quantitative plots                | `Figure_4_FRAP/Figure_4_plots/` and `Figure_4_FRAP/raw_data/`                                                             |
| FRAP representative images             | `publication_images/` under the raw FRAP input folder                                                                       |
| Localization images                    | The translocation input folder; interactive crops in its `tif_images/` subfolders                                           |
| Figures 5, 6C, 7C-I, S7E-F, S9A        | `figures/` and `raw_data/` within each notebook folder                                                                    |
| Experimental ACF preprocessing         | `optimization/source_files_ACF/data_ACF/`                                                                                   |
| HT linear fits                         | Under `Figure_7_HT/`: `figures_linear/`, `raw_data_linear/`, and `linear_fit_tables/`                                 |
| Figure 8 single-configuration analysis | `Figure_8_Multiplexing/results_normal_diffusion/` and `raw_data/`; filtered configurations use threshold-specific folders |
| Figure 8 aggregate comparison          | `Figure_8_Multiplexing/results_aggregate_comparison_0p5um/` with separate results, figures, and data for each variant       |
| Figure 8 representative images         | `Figure_8_Multiplexing/figures/` by default                                                                                 |
| S10 simulations                        | `Figure_Sup_SSA/results_simulations/`                                                                                       |
| S11 cost surfaces                      | `optimization/optimization_parameters/`                                                                                     |

The HT, S7E-F, and S9A notebooks load tracking data with `src/utilities/inhibitors.py`. The 7J-M linear fits are calculated in `Figure_7_J_K_L_M__HT_linear_fit.ipynb` itself; S7E-F and S9A also use `inhibitors.py` for analysis and plotting.

## Dependencies

| Purpose                 | Packages                                                    |
| ----------------------- | ----------------------------------------------------------- |
| Microscopy analysis     | [microlive](https://github.com/ningzhaoAnschutz/microlive)   |
| TASEP simulations       | [tasep_models](https://pypi.org/project/tasep-models/)       |
| Image analysis          | scikit-image, big-fish, trackpy, tifffile, readlif          |
| Segmentation            | PyTorch, Cellpose                                           |
| Numerics and statistics | NumPy, SciPy, statsmodels, Numba                            |
| Plotting                | matplotlib, seaborn, joypy, brokenaxes, matplotlib-scalebar |
| Sequences               | dna-features-viewer, snapgene-reader, Biopython             |

See [environment.yml](environment.yml) for the full list.

## License

Code is released under the [BSD 3-Clause License](LICENSE), copyright (c) 2025 Ning Zhao. The microscopy data license is given in the Zenodo record.
