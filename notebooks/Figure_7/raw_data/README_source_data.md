# Figure 7C–M source data

## Scope

This directory contains processed figure source data for the eleven quantitative panels in Figure 7 (7C–M). These files are not unprocessed microscopy data. Each table was exported from the same prepared arrays used to draw its corresponding native notebook panel.

Panels 7A–B are outside this package. The notebooks, panel conditions, and source-data files are:

| Panel | Notebook | Condition | Source-data file |
|---|---|---|---|
| 7C | `Figure_7_C_D_E_spot_properties.ipynb` | Four conditions | `Figure_7C_spot_intensity_individual_cell_means.csv` |
| 7D | `Figure_7_C_D_E_spot_properties.ipynb` | Four conditions | `Figure_7D_spot_SNR_individual_cell_means.csv` |
| 7E | `Figure_7_C_D_E_spot_properties.ipynb` | Four conditions | `Figure_7E_spot_size_individual_cell_means.csv` |
| 7F | `Figure_7_F_G_H_I__ACF.ipynb` | UTag | `Figure_7F_UTag_ACF_curves.csv` |
| 7G | `Figure_7_F_G_H_I__ACF.ipynb` | UTag(DeltaCys) | `Figure_7G_UTag_DeltaCys_ACF_curves.csv` |
| 7H | `Figure_7_F_G_H_I__ACF.ipynb` | SunTag | `Figure_7H_SunTag_ACF_curves.csv` |
| 7I | `Figure_7_F_G_H_I__ACF.ipynb` | ALFA-tag | `Figure_7I_ALFA_tag_ACF_curves.csv` |
| 7J | `Figure_7_J_K_L_M__HT.ipynb` | UTag | `Figure_7J_UTag_HT_curves.csv` |
| 7K | `Figure_7_J_K_L_M__HT.ipynb` | UTag(DeltaCys) | `Figure_7K_UTag_DeltaCys_HT_curves.csv` |
| 7L | `Figure_7_J_K_L_M__HT.ipynb` | SunTag | `Figure_7L_SunTag_HT_curves.csv` |
| 7M | `Figure_7_J_K_L_M__HT.ipynb` | ALFA-tag | `Figure_7M_ALFA_tag_HT_curves.csv` |

Canonical ASCII publication labels are `UTag`, `UTag_DeltaCys`, `SunTag`, and `ALFA_tag`. They correspond to the displayed labels UTag, UTag(DeltaCys) or DeltaCys, SunTag, and ALFA-tag, respectively.

## Figures 7C–E: spot properties

The primary input root is `/Volumes/LaCie/UTag_paper_data/ACF`. If that volume is unavailable, the notebook may use only the byte-identical submitted-paper OneDrive copy at `General - Zhao (NZ) Lab/Zhao lab shared folder/Our papers/Submitted/UTag paper/Raw_Data_ACF_HT/ACF_data`. One source root is selected for the complete run and is treated as read-only.

Each CSV row represents one analyzed cell, corresponding to one included `tracking_*.csv` file and one black dot in the panel. Missing source values are removed, then the arithmetic mean within that cell file is calculated. The included counts are 19 UTag, 14 UTag_DeltaCys, 16 SunTag, and 11 ALFA_tag cells (60 rows per panel).

| File | Columns | Calculation |
|---|---|---|
| `Figure_7C_spot_intensity_individual_cell_means.csv` | `condition`, `mean_spot_intensity_au` | Mean of nonmissing `spot_int_ch_0`, in arbitrary units |
| `Figure_7D_spot_SNR_individual_cell_means.csv` | `condition`, `mean_spot_snr_au` | Mean of nonmissing `snr_ch_0`, in arbitrary units |
| `Figure_7E_spot_size_individual_cell_means.csv` | `condition`, `mean_spot_size_um` | Mean of nonmissing `spot_size_ch_0` after conversion to micrometers |

Spot size is converted using `spot_size_ch_0 × 0.12989318982387477 µm/pixel`. The native plots show every per-cell value, a red median, an interquartile-range box, and whiskers at the 5th and 95th percentiles. Figure 7E remains generated with `show_stats=False`; its final-composite significance bars and stars were added during figure assembly and are not notebook or source-data outputs. Significance tests and p-values are excluded from these files.

## Figures 7F–I: autocorrelation curves

Each ACF file contains 295 displayed coordinates from 0 through 24.5 minutes inclusive at 5-second (1/12-minute) intervals. One row represents one lag coordinate for the experimental and model curves in one panel.

The shared schema is:

| Column | Definition |
|---|---|
| `lag_time_min` | Autocorrelation lag in minutes |
| `mean_experimental_autocorrelation` | Displayed experimental mean autocorrelation |
| `bootstrap_se_experimental_autocorrelation` | Symmetric bootstrap standard error used for the experimental band |
| `n_experimental_trajectories` | Post-filtering experimental trajectory count |
| `mean_model_autocorrelation` | Displayed mean simulated autocorrelation |
| `bootstrap_se_model_autocorrelation` | Symmetric bootstrap standard error used for the model band |
| `n_model_repetitions` | Simulated trajectory count (200) |

The experimental coordinates are read from the four local `optimization/results_ACF/df_ACF_*.csv` tables used by the final plotting notebook. Channel-0 spot-intensity trajectories are filtered at mean SNR ≥ 0.5, require at least 30% coverage, allow at most five missing frames during shifting, and are padded or truncated to 360 frames. Missing values are forward-filled within valid trajectory limits. Correlation outliers beyond 4 median absolute deviations are removed. Baseline correction and linear projection at lag zero are enabled. The stored error values are `np.nanstd` across 1,000 resampled mean-correlation curves and therefore estimate the standard error of the mean autocorrelation, despite the upstream column name `std_correlation`.

Experimental post-filtering counts are:

| Panel | Condition | Trajectories |
|---|---|---:|
| 7F | UTag | 302 |
| 7G | UTag_DeltaCys | 166 |
| 7H | SunTag | 136 |
| 7I | ALFA_tag | 157 |

The TASEP model uses 200 repetitions, a 2,000-second burn-in, a 2,000-second post-burn-in duration, a 5-second time grid, the same 4-MAD correlation filter, baseline correction, lag-zero projection, and 1,000 bootstrap resamples. `REPRODUCIBLE_RANDOM_SEED = 42` is the single configured seed; deterministic child RNG states for both parallel SSA repetitions and the fixed-order bootstrap are derived from it. The kinetic inputs are:

| Condition | `ki` (1/s) | `ke` (aa/s) |
|---|---:|---:|
| UTag | 0.03126582278481013 | 3.113924050632911 |
| UTag_DeltaCys | 0.06316455696202532 | 5.291139240506329 |
| SunTag | 0.033924050632911394 | 4.329113924050633 |
| ALFA_tag | 0.042784810126582286 | 4.8354430379746836 |

Autocorrelation describes the similarity of a translation-site intensity trajectory to a time-shifted version of itself. Hidden lags beyond the panel's 24.5-minute x-axis are retained for plotting behavior but are not exported.

## Figures 7J–M: harringtonine response

The required input root is `/Volumes/LaCie/UTag_paper_data/Harringtonine`; no OneDrive fallback is used. Conditions map to `UTag/HT_Analysis_GUI`, `UTag_CF/HT_Analysis_GUI`, `SunTag/HT_Analysis_GUI`, and `AlfaTag/HT_Analysis_GUI`. Only regular `.csv` files beginning with `tracking_` are accepted. AppleDouble files beginning with `._` are rejected, and each included `results_*` directory must contain exactly one real tracking CSV.

Each HT CSV contains 30 time coordinates from −5 through 24 minutes relative to harringtonine application. One row represents one displayed time coordinate. The shared schema is:

| Column | Definition |
|---|---|
| `time_min` | Time relative to harringtonine application, in minutes |
| `mean_normalized_experimental_intensity` | Mean normalized intensity across responding cells |
| `population_sd_normalized_experimental_intensity` | Population SD (`ddof=0`) used for the experimental band |
| `n_responding_cells` | Responding-cell count contributing to the experimental summary |
| `mean_model_fit_normalized_intensity` | Displayed best-candidate model mean |
| `n_model_repetitions` | Simulations contributing to the model mean (100) |

For each cell and frame, the notebook counts unique particles and sums `spot_int_ch_0`. During the five pre-treatment frames, each sum is divided by that frame's particle count. At and after treatment, each sum is divided by the mean pre-treatment particle count. The complete trajectory is then divided by its mean pre-treatment intensity. Frame 5 is the application frame and is shifted to 0 minutes.

A responding cell has a mean over its final five normalized frames below 30% of its pre-treatment baseline. Only responders contribute to the orange mean and band:

| Panel | Condition | Total cells | Responding cells |
|---|---|---:|---:|
| 7J | UTag | 38 | 32 |
| 7K | UTag_DeltaCys | 24 | 22 |
| 7L | SunTag | 33 | 32 |
| 7M | ALFA_tag | 32 | 32 |

The experimental center is the arithmetic mean and the band is population SD, not SEM. The cyan region from −1 to 0 minutes marks HT application and is an annotation rather than a measured variable.

For the model, `ki` is fixed at 0.03 1/s. Fifteen candidate `ke` values from 1.0 to 8.0 aa/s in 0.5-aa/s increments are evaluated with 95% inhibition effectiveness, 100 repetitions per candidate, a 2,000-second burn-in, and 1-second simulation steps downsampled to 60 seconds. Seed 42 is deterministically derived per condition, candidate rate, and repetition. The selected candidate minimizes the 0-to-24-minute weighted residual sum, `Σ[(experimental mean − model mean)² / (experimental population SD² + 10⁻⁶)]`. This is a finite candidate search, so there are no iterative-fit convergence states or failed-fit exclusions.

| Panel | Condition | Winning `ke` (aa/s) |
|---|---|---:|
| 7J | UTag | 3.5 |
| 7K | UTag_DeltaCys | 3.0 |
| 7L | SunTag | 4.0 |
| 7M | ALFA_tag | 3.0 |

The current deterministic seed-42 execution is the regeneration reference for the gray HT model lines. Relative to the archived April panel images, the gray curves have a small model-only displacement (median 0.0049–0.0083 and 95th percentile 0.0146–0.0250 normalized-intensity units across panels; maximum 0.0194–0.0278). The experimental means, population-SD bands, responder counts, and selected `ke` values are unchanged; the regenerated panels should not be described as pixel-identical to the April images.

The complete 15-rate diagnostic tables are written separately to `../optimization_tables/`; they are not publication source-data CSVs.

## Explicit exclusions

This package does not export panels 7A–B; image provenance; raw microscopy images or LIF metadata; source tracking filenames, folder names, experiment IDs, or acquisition dates; cell identifiers for the unpaired 7C–E plots; raw per-spot or per-frame rows; box medians, quartiles, whiskers, swarm jitter positions, p-values, Mann–Whitney statistics, significance bars, or stars; individual ACF correlations or TASEP trajectories; hidden ACF coordinates outside 0–24.5 minutes; individual or excluded nonresponding HT trajectories; HT diagnostic linear fits; the 15-rate cost sweep within `raw_data/` or `figures/`; ACF cost surfaces or parameter-search grids; plotting colors, marker sizes, axis coordinates, legends, or a combined Figure 7 canonical dataset.
