# Figure 7C-M source-data export plan

## Status and scope

This is a review plan only. It does not modify any notebook, source CSV,
processed ACF table, model parameter file, plot, or other figure folder.

The requested scope is the eleven quantitative panels in final manuscript
Figure 7:

| Panel | Notebook | Final plotted content |
|---|---|---|
| Figure 7C | `Figure_7_C_D_E_spot_properties.ipynb` | Per-cell mean translation-spot intensity |
| Figure 7D | `Figure_7_C_D_E_spot_properties.ipynb` | Per-cell mean translation-spot SNR |
| Figure 7E | `Figure_7_C_D_E_spot_properties.ipynb` | Per-cell mean translation-spot size in micrometers |
| Figure 7F | `Figure_7_F_G_H_I__ACF.ipynb` | UTag experimental and model autocorrelation curves |
| Figure 7G | `Figure_7_F_G_H_I__ACF.ipynb` | UTag(DeltaCys) experimental and model autocorrelation curves |
| Figure 7H | `Figure_7_F_G_H_I__ACF.ipynb` | SunTag experimental and model autocorrelation curves |
| Figure 7I | `Figure_7_F_G_H_I__ACF.ipynb` | ALFA-tag experimental and model autocorrelation curves |
| Figure 7J | `Figure_7_J_K_L_M__HT.ipynb` | UTag harringtonine-response mean, SD band, and model fit |
| Figure 7K | `Figure_7_J_K_L_M__HT.ipynb` | UTag(DeltaCys) harringtonine-response mean, SD band, and model fit |
| Figure 7L | `Figure_7_J_K_L_M__HT.ipynb` | SunTag harringtonine-response mean, SD band, and model fit |
| Figure 7M | `Figure_7_J_K_L_M__HT.ipynb` | ALFA-tag harringtonine-response mean, SD band, and model fit |

Panels 7A-B are outside this request. No image-provenance table is proposed.

The files described below are processed **figure source data**, not
acquisition-level raw microscopy data. The README must make that distinction.

## Fixed scope

- Export one minimal source-data CSV for each quantitative panel (7C-M), plus
  one README, inside `notebooks/Figure_7/raw_data/`.
- Use LaCie as the required 7J-M input. For 7C-E, use LaCie as the primary
  input and the byte-identical submitted-paper OneDrive copy only as an
  explicit fallback.
- Export exactly the 295 displayed ACF coordinates from 0 through 24.5 minutes
  inclusive, while retaining the complete arrays for plotting.
- Keep `Figure 7E` notebook `show_stats=False`. Its bars/stars were added in
  final-figure assembly and are not notebook or source-data outputs.
- Preserve the existing source data, calculations, fitted parameter inputs,
  panel layouts, and plot settings. The approved F-I deterministic model seed
  correction is the only intended model change.
- For 7J-M, expect small, documented differences in the gray model curve when
  regenerating from the current seeded reference; experimental data, SD bands,
  responder counts, and selected best `ke` values must remain unchanged.

The remaining verification gates are stated in the validation section rather
than repeated as resolved decision history.

---

# Recommended folder structure

All outputs should be generated inside
`notebooks/Figure_7/`:

```text
notebooks/Figure_7/
├── Figure_7_C_D_E_spot_properties.ipynb
├── Figure_7_F_G_H_I__ACF.ipynb
├── Figure_7_J_K_L_M__HT.ipynb
├── Figure_7_C_M_source_data_export_plan.md
├── figures/
│   ├── Figure_7C_spot_intensity.png
│   ├── Figure_7C_spot_intensity.svg
│   ├── Figure_7D_spot_SNR.png
│   ├── Figure_7D_spot_SNR.svg
│   ├── Figure_7E_spot_size_um.png
│   ├── Figure_7E_spot_size_um.svg
│   ├── Figure_7F_UTag_ACF.png
│   ├── Figure_7F_UTag_ACF.svg
│   ├── Figure_7G_UTag_DeltaCys_ACF.png
│   ├── Figure_7G_UTag_DeltaCys_ACF.svg
│   ├── Figure_7H_SunTag_ACF.png
│   ├── Figure_7H_SunTag_ACF.svg
│   ├── Figure_7I_ALFA_tag_ACF.png
│   ├── Figure_7I_ALFA_tag_ACF.svg
│   ├── Figure_7J_UTag_HT.png
│   ├── Figure_7J_UTag_HT.svg
│   ├── Figure_7K_UTag_DeltaCys_HT.png
│   ├── Figure_7K_UTag_DeltaCys_HT.svg
│   ├── Figure_7L_SunTag_HT.png
│   ├── Figure_7L_SunTag_HT.svg
│   ├── Figure_7M_ALFA_tag_HT.png
│   └── Figure_7M_ALFA_tag_HT.svg
├── raw_data/
│   ├── Figure_7C_spot_intensity_individual_cell_means.csv
│   ├── Figure_7D_spot_SNR_individual_cell_means.csv
│   ├── Figure_7E_spot_size_individual_cell_means.csv
│   ├── Figure_7F_UTag_ACF_curves.csv
│   ├── Figure_7G_UTag_DeltaCys_ACF_curves.csv
│   ├── Figure_7H_SunTag_ACF_curves.csv
│   ├── Figure_7I_ALFA_tag_ACF_curves.csv
│   ├── Figure_7J_UTag_HT_curves.csv
│   ├── Figure_7K_UTag_DeltaCys_HT_curves.csv
│   ├── Figure_7L_SunTag_HT_curves.csv
│   ├── Figure_7M_ALFA_tag_HT_curves.csv
│   └── README_source_data.md
└── optimization_tables/
    ├── ke_optimization_analysis_UTag.csv
    ├── ke_optimization_analysis_UTag_C_Free.csv
    ├── ke_optimization_analysis_SunTag.csv
    └── ke_optimization_analysis_AlfaTag.csv
```

This gives one CSV per final quantitative image/panel, eleven CSVs total. The
panel and single-condition identity are encoded in each filename, avoiding
redundant `figure_panel` and `condition` columns in the single-condition ACF
and HT tables.

The existing `ACF_plots/` and `HT_plots/` folders can remain as archived or
analysis outputs. The publication notebooks should write their reviewed final
PNG/SVG files to `figures/` and their source-data tables only to `raw_data/`.
The J-M diagnostic optimization cost tables should be written only to
`optimization_tables/`.

The notebooks generate native quantitative panels. The supplied assembled
Figure 7 additionally contains a shared/manual legend and manual Figure 7E
significance annotations. Those assembly elements are not generated or
exported by the three notebooks and are not numerical source data.

---

# Inputs, outputs, and reproducibility boundaries

## Figure 7C-E notebook

The active notebook uses:

```python
path_main_folder = Path('/Volumes/LaCie/UTag_paper_data/ACF')
results_folder = path_main_folder.joinpath('results_Spots_Properties')
```

These two variables have different roles:

- `path_main_folder` is the input root. `dataset_selection` reads the four
  condition folders `UTag_ACF`, `UTag_CF_ACF`, `SunTag_ACF`, and
  `AlfaTag_ACF` beneath it.
- `results_folder` is an output folder only. It receives the current C-E PNG
  and SVG files; it does not generate or supply the plotted measurements.

The current notebook does not read OneDrive. The 60 included input files
(19 UTag, 14 UTag_DeltaCys, 16 SunTag, and 11 ALFA_tag) are byte-identical to
the submitted-paper OneDrive fallback: matching relative paths, sizes, and
SHA-256 hashes. Select exactly one source root per run and keep all source
folders read-only. New plots and CSVs must be written only under
`notebooks/Figure_7/`.

For strict native-panel pixel preservation, retain the current traversal order
when constructing the plotted C-E arrays. Sorting the same files preserves all
values and box statistics but changed five antialiased pixels in the 7E swarm
test. Sorting may be used for a checksum manifest or validation only.

## Figure 7F-I notebook

The final ACF notebook reads the exact plotted experimental summaries from:

```text
optimization/results_ACF/df_ACF_utag.csv
optimization/results_ACF/df_ACF_utag_c_free.csv
optimization/results_ACF/df_ACF_suntag.csv
optimization/results_ACF/df_ACF_alfatag.csv
```

These four processed tables are the immediate experimental source of truth for
Figures 7F-I. Each contains:

```text
lags,mean_correlation,std_correlation
```

Each file has 360 lag rows at 5-second intervals. The final notebook converts
lags from seconds to minutes and displays the range from 0 to 24.5 minutes.

The local `mean_correlation` values are identical to the current submitted
OneDrive ACF result files. The `std_correlation` values differ slightly because
the upstream bootstrap was stochastic. Therefore, implementation must use the
four local tables actually read by the final plotting notebook and must not
replace their error values with a fresh upstream ACF calculation.

The model curves are newly simulated each time the notebook runs. The notebook
currently has no explicit random seed for either the ACF TASEP simulations or
their bootstrap uncertainty. This was confirmed through all relevant layers:

- `Figure_7_F_G_H_I__ACF.ipynb` does not call `np.random.seed`,
  `random.seed`, or a seeded simulation wrapper.
- Its `run_simulation` function calls `simulate_TASEP_SSA` directly. That
  parallel wrapper has no seed argument, and its Numba SSA implementation uses
  `np.random.rand()` to generate events.
- The installed `microlive.Correlation` bootstrap creates
  `np.random.default_rng()` without a seed for each bootstrap iteration.
  Consequently, adding only a global `np.random.seed(...)` line would not seed
  this bootstrap.

The deterministic seed code the user recalled is present in
`Figure_7_J_K_L_M__HT.ipynb`, where `REPRODUCIBLE_RANDOM_SEED = 42` and a
seeded TASEP wrapper are used. It is not currently present in the ACF notebook.

The four orange experimental ACF curves and their orange uncertainty bands are
fixed local CSV inputs and will repeat exactly. The newly simulated gray model
mean and its gray bootstrap band are not currently guaranteed to repeat. This
can produce small model-only differences between executions even when code and
parameters are unchanged.

The approved deterministic correction will expose one configured seed only:

```python
REPRODUCIBLE_RANDOM_SEED = 42
```

Both stochastic model stages must derive their RNG state from that single base
seed. The parallel SSA repetitions will use the same notebook-local stable child
seed strategy already used by the HT notebook so that worker scheduling does not
change the output. The model bootstrap will use a notebook-local deterministic
bootstrap helper driven from the same base seed, with 1,000 resamples generated
in a fixed order. It must preserve the current resampling, baseline-correction,
`np.nanstd`, and lag-zero policies and must not call the installed unseeded
`microlive.Correlation` bootstrap. No second user-configured seed is introduced.
A global `np.random.seed(...)` line by itself is insufficient.

The single seeded execution defines the publication model arrays. No comparison
with a previous stochastic gray model curve is required. The orange experimental
values must remain unchanged.

The ACF output directory is derived from `Path().resolve()`. Its location can
therefore depend on the directory from which Jupyter was launched. It should be
replaced with a repository-derived Figure 7 path during implementation.

## Figure 7J-M notebook

The active notebook points to an obsolete, now-missing OneDrive hierarchy. The
final optimized HT plotting function must instead use the required LaCie root:

```text
/Volumes/LaCie/UTag_paper_data/Harringtonine
```

The exact condition directories under that root are:

| Notebook condition | Required LaCie directory |
|---|---|
| `utag` | `UTag/HT_Analysis_GUI` |
| `utag_c_free` | `UTag_CF/HT_Analysis_GUI` |
| `suntag` | `SunTag/HT_Analysis_GUI` |
| `alfatag` | `AlfaTag/HT_Analysis_GUI` |

Use regular CSV files whose names begin with `tracking_`; explicitly reject
filenames beginning with `._` and assert that every included `results_*`
directory contains exactly one real tracking CSV. This is required because the
LaCie UTag tree currently contains 38 real tracking CSVs plus 38 AppleDouble
`._tracking_*.csv` sidecars. After excluding sidecars, the selected source has
127 real tracking CSVs: 38 UTag, 24 UTag_DeltaCys, 33 SunTag, and 32 ALFA_tag.

The final optimized HT plotting function must write panel PNG/SVG files only to
`notebooks/Figure_7/figures/`. The four optimization cost tables must be written
only to `notebooks/Figure_7/optimization_tables/`. The notebook currently
targets the wrong Figure 6 folder. The preceding diagnostic pass writes
separate figures to `optimization/results_ACF` and is not a final-panel source.

The J-M simulation wrapper is already seeded. A controlled current-code run
against the LaCie source reproduced the stored May notebook results to the
recorded precision and selected the same best `ke` values as the April panels
(3.5, 3.0, 4.0, 3.0 aa/s). The gray model curves are not pixel-identical to the
April reference: their median vertical movement is 0.0049-0.0083 normalized
intensity units and their 95th-percentile movement is 0.0146-0.0250 units
(maximum 0.0194-0.0278). This is a small model-only change on the 0-1.4
y-axis; experimental points and population-SD bands remain exact. Retain the
April images and cost tables as provenance, but use the current seeded run as
the regeneration reference rather than claiming April pixel identity.

All three notebooks currently have a Jupyter kernel name of `python3` with a
display name of `microlive`. Implementation should set the actual kernel name
to the installed `microlive` kernel, consistent with Figures 4-6.

---

# Canonical condition labels

Use these ASCII labels consistently in all Figure 7 CSVs and documentation:

| CSV label | Displayed label in Figure 7 |
|---|---|
| `UTag` | UTag |
| `UTag_DeltaCys` | DeltaCys or UTag(DeltaCys) |
| `SunTag` | SunTag |
| `ALFA_tag` | ALFA or ALFA-tag |

Do not use Greek delta characters, TeX strings, `utag_c_free`, `UTag_C_Free`,
or `AlfaTag` as publication CSV labels.

---

# Figures 7C-E: spot-property box-and-swarm plots

## Unit of observation

For each condition, the notebook loads every included `tracking_*.csv`. It
extracts the selected measurement from all rows, removes missing values, and
computes one arithmetic mean per tracking file. The resulting mean is one
black dot in the final panel.

The folder organization and current reports indicate that one included
tracking file represents one analyzed cell. Thus, each CSV row should represent
one analyzed cell and one displayed black dot—not one RNA spot, frame, or raw
tracking-table row.

No cell identifier is required to reproduce the unpaired box-and-swarm plots.
The filename identifies the panel, while `condition` identifies the x-axis
group.

## Source counts

| Condition | Included tracking files/cells | Dots in each of 7C, 7D, and 7E |
|---|---:|---:|
| `UTag` | 19 | 19 |
| `UTag_DeltaCys` | 14 | 14 |
| `SunTag` | 16 | 16 |
| `ALFA_tag` | 11 | 11 |
| **Total** | **60** | **60** |

The three panels use the same 60 source files but calculate different
measurements from them.

## Figure 7C calculation

Source field:

```text
spot_int_ch_0
```

One plotted value is:

\[
\text{mean spot intensity for one cell}
= \operatorname{mean}(\texttt{spot\_int\_ch\_0})
\]

Missing source values are removed before the mean is calculated. The displayed
unit is arbitrary units.

### Proposed CSV

`Figure_7C_spot_intensity_individual_cell_means.csv`

Expected rows: 60.

| Column | Definition |
|---|---|
| `condition` | Canonical ASCII condition label |
| `mean_spot_intensity_au` | Arithmetic mean of all nonmissing `spot_int_ch_0` values in one included cell tracking file, in arbitrary units |

### Mock Figure 7C table

Illustrative values only; do not copy them into the real CSV.

| condition | mean_spot_intensity_au |
|---|---:|
| UTag | 2142.61 |
| UTag_DeltaCys | 1987.34 |
| SunTag | 2310.18 |
| ALFA_tag | 2256.47 |

## Figure 7D calculation

Source field:

```text
snr_ch_0
```

One plotted value is the arithmetic mean of all nonmissing channel-0 spot SNR
values in one included cell tracking file. The final y-axis reports SNR in
arbitrary units.

### Proposed CSV

`Figure_7D_spot_SNR_individual_cell_means.csv`

Expected rows: 60.

| Column | Definition |
|---|---|
| `condition` | Canonical ASCII condition label |
| `mean_spot_snr_au` | Arithmetic mean of all nonmissing `snr_ch_0` values in one included cell tracking file |

### Mock Figure 7D table

Illustrative values only.

| condition | mean_spot_snr_au |
|---|---:|
| UTag | 2.84 |
| UTag_DeltaCys | 2.63 |
| SunTag | 3.17 |
| ALFA_tag | 3.42 |

## Figure 7E calculation

Source field:

```text
spot_size_ch_0
```

The notebook converts every nonmissing spot-size value from pixels to
micrometers using:

\[
\text{spot size in micrometers}
= \texttt{spot\_size\_ch\_0} \times 0.12989318982387477
\]

It then calculates one arithmetic mean per included cell tracking file. The
conversion factor is the microscope pixel size in micrometers per pixel.

### Proposed CSV

`Figure_7E_spot_size_individual_cell_means.csv`

Expected rows: 60.

| Column | Definition |
|---|---|
| `condition` | Canonical ASCII condition label |
| `mean_spot_size_um` | Arithmetic mean channel-0 spot size for one included cell after conversion to micrometers |

### Mock Figure 7E table

Illustrative values only.

| condition | mean_spot_size_um |
|---|---:|
| UTag | 0.3271 |
| UTag_DeltaCys | 0.3198 |
| SunTag | 0.3612 |
| ALFA_tag | 0.3235 |

## Box and whisker definitions for 7C-E

All three panels display:

- every per-cell mean as a black swarm point;
- a red median line;
- a box spanning the 25th to 75th percentiles;
- whiskers at the 5th and 95th percentiles;
- all dots, including values outside the whiskers.

Figures 7C-D have no significant comparison bars in the final composite.
The Figure 7E bars and stars were manually added during final-figure assembly;
the notebook intentionally uses `show_stats=False` and does not generate those
annotations. They are outside the source-data export and native-panel
regression scope.

The CSVs should contain only the 60 individual values represented by the dots.
Do not export medians, quartiles, whiskers, swarm x-coordinates, p-values,
Mann-Whitney statistics, or significance labels.

The notebook prints group mean plus a population-SD-divided-by-square-root-n
quantity, but those summaries are not displayed by the boxplots and should not
become publication source-data columns.

---

# Figures 7F-I: autocorrelation data and model curves

## Panel-to-condition map

The mapping follows the final notebook loop order and the panel layout in the
provided composite:

| Panel | Canonical condition | Current plot name |
|---|---|---|
| Figure 7F | `UTag` | `utag` |
| Figure 7G | `UTag_DeltaCys` | `utag_c_free` |
| Figure 7H | `SunTag` | `suntag` |
| Figure 7I | `ALFA_tag` | `alfatag` |

## Experimental ACF values

The plotted experimental series is read directly from the corresponding local
`df_ACF_*.csv` file:

- `lags` is in seconds and is divided by 60 for display in minutes;
- `mean_correlation` is the displayed orange experimental mean;
- `std_correlation` is the magnitude used above and below the mean for the
  orange uncertainty band.

The upstream `ACF_calculation_data_and_controls.ipynb` establishes the relevant
processing:

- channel-0 spot intensity is used;
- trajectories with mean SNR below 0.5 are removed;
- trajectories require at least 30% coverage and allow at most five missing
  frames during shifting;
- trajectories are padded or truncated to 360 frames;
- the acquisition interval is 5 seconds;
- missing values are forward-filled within valid trajectory limits;
- correlation outliers beyond 4 median absolute deviations are removed;
- baseline correction is enabled;
- lag zero is replaced by linear projection;
- 1,000 bootstrap resamples are used for uncertainty.

Although the source column is named `std_correlation`, the ACF-generation
notebook invokes `microlive.Correlation(..., use_bootstrap=True)`. The installed
implementation creates 1,000 resampled mean-correlation curves and calculates
`np.nanstd` across those curves at each lag. It therefore represents a
**bootstrap estimate of the standard error of the mean autocorrelation**, not
the standard deviation of individual trajectory correlations. This definition
must be documented explicitly in the README.

Current submitted-source reports give these post-filtering trajectory counts:

| Panel | Condition | Experimental trajectories after filtering |
|---|---|---:|
| Figure 7F | `UTag` | 302 |
| Figure 7G | `UTag_DeltaCys` | 166 |
| Figure 7H | `SunTag` | 136 |
| Figure 7I | `ALFA_tag` | 157 |

These counts should be validated against the exact ACF inputs associated with
the four local summary tables before they are written into the publication
tables or README.

## Model ACF values

The gray model line and gray uncertainty band are generated by a TASEP SSA
simulation and the same correlation workflow. Current model settings are:

| Setting | Value |
|---|---:|
| Simulation repetitions | 200 |
| Burn-in time | 2,000 s |
| Simulated duration after burn-in | 2,000 s |
| Simulation/correlation interval | 5 s |
| Correlation outlier threshold | 4 MAD |
| Bootstrap iterations | 1,000 |
| Baseline correction | enabled |
| Lag-zero linear projection | enabled |

The model mean is `mean_correlation_ssa`. The model band uses
`std_correlation_ssa`, which is likewise a bootstrap estimate of the standard
error of the model mean autocorrelation.

The fitted kinetic parameters read by the notebook are:

| Condition | `ki` (1/s) | `ke` (aa/s) |
|---|---:|---:|
| `UTag` | 0.03126582278481013 | 3.113924050632911 |
| `UTag_DeltaCys` | 0.06316455696202532 | 5.291139240506329 |
| `SunTag` | 0.033924050632911394 | 4.329113924050633 |
| `ALFA_tag` | 0.042784810126582286 | 4.8354430379746836 |

These parameters and model settings belong in the README, not in separate
publication CSVs, because the final panels display the curve coordinates rather
than a parameter table.

## Approved ACF export range

The experimental input has 360 rows from 0 to 29.9167 minutes, and the model
has 400 simulated lag points. The final axes display only 0 to 24.5 minutes.

Consistent with the approved requirement to export only values displayed in
the figure, every Figure 7F-I table will contain the 295 common time points
from 0 to 24.5 minutes inclusive at 5-second intervals. Hidden coordinates
beyond the final x-axis limit will not be included.

## Proposed ACF files and schemas

The same schema applies to:

- `Figure_7F_UTag_ACF_curves.csv`
- `Figure_7G_UTag_DeltaCys_ACF_curves.csv`
- `Figure_7H_SunTag_ACF_curves.csv`
- `Figure_7I_ALFA_tag_ACF_curves.csv`

Expected rows: exactly 295 per file.

| Column | Definition |
|---|---|
| `lag_time_min` | Displayed autocorrelation lag in minutes |
| `mean_experimental_autocorrelation` | Experimental mean autocorrelation displayed as the orange point/line series |
| `bootstrap_se_experimental_autocorrelation` | Symmetric bootstrap standard error used for the orange band |
| `n_experimental_trajectories` | Number of post-filtering experimental trajectories contributing to the curve |
| `mean_model_autocorrelation` | Mean simulated autocorrelation displayed as the gray model curve |
| `bootstrap_se_model_autocorrelation` | Symmetric bootstrap standard error used for the gray model band |
| `n_model_repetitions` | Number of simulated trajectories contributing to the model curve; currently 200 |

Separate lower and upper band columns are redundant because each band is
plotted as mean minus and plus the same error magnitude.

### Mock Figure 7F-I ACF table

Illustrative values only. Each panel will have its own CSV and its own actual
values.

| lag_time_min | mean_experimental_autocorrelation | bootstrap_se_experimental_autocorrelation | n_experimental_trajectories | mean_model_autocorrelation | bootstrap_se_model_autocorrelation | n_model_repetitions |
|---:|---:|---:|---:|---:|---:|---:|
| 0.000000 | 0.0762 | 0.0000 | 302 | 0.0481 | 0.0000 | 200 |
| 0.083333 | 0.0680 | 0.0025 | 302 | 0.0469 | 0.0012 | 200 |
| 0.166667 | 0.0613 | 0.0024 | 302 | 0.0458 | 0.0012 | 200 |

The mock uses an illustrative Figure 7F-like count only to demonstrate the
schema. It must not be copied into any final CSV.

---

# Figures 7J-M: harringtonine-response data and optimized model curves

## Panel-to-condition map

| Panel | Canonical condition | Current optimization plot name |
|---|---|---|
| Figure 7J | `UTag` | `UTag` |
| Figure 7K | `UTag_DeltaCys` | `UTag_C_Free` |
| Figure 7L | `SunTag` | `SunTag` |
| Figure 7M | `ALFA_tag` | `AlfaTag` |

The final composite uses the four plots generated by
`evaluate_harringtonine_ke_optimization`, not the earlier diagnostic plots from
`process_HT_data`.

## Experimental trajectory calculation

Each included cell tracking file is reduced to a 30-frame normalized
trajectory:

1. At each frame, count unique particles and sum `spot_int_ch_0` across those
   particles.
2. For the five pre-treatment frames, divide the summed intensity at each frame
   by that frame's particle count.
3. Calculate the mean particle count across the five pre-treatment frames.
4. For treatment and post-treatment frames, divide summed intensity by that
   fixed mean pre-treatment particle count.
5. Concatenate the pre- and post-treatment values.
6. Divide the complete trajectory by the mean of the five pre-treatment
   intensity values.
7. Shift time so frame 5, the harringtonine application frame, is 0 minutes.

This produces time points from -5 through 24 minutes at one-minute intervals.

## Responder filter

Only responding cells contribute to the displayed orange mean and band.

For each cell:

- baseline is the mean of normalized frames 0-4 before time shifting;
- the response threshold is 30% of that baseline;
- the final response is the mean of the last five frames;
- the cell is retained if the final response is below the threshold.

Current counts are:

| Panel | Condition | Total cells | Responding cells used in plot | Excluded nonresponding cells |
|---|---|---:|---:|---:|
| Figure 7J | `UTag` | 38 | 32 | 6 |
| Figure 7K | `UTag_DeltaCys` | 24 | 22 | 2 |
| Figure 7L | `SunTag` | 33 | 32 | 1 |
| Figure 7M | `ALFA_tag` | 32 | 32 | 0 |

The orange center line is the arithmetic mean across responding cells at each
time point. The orange shaded band is mean plus/minus **population standard
deviation**, calculated by `np.std(..., axis=0)` with `ddof=0`, because the
final optimization call uses `use_sem=False`.

The individual cell trajectories are not displayed and therefore are not
proposed as separate publication CSVs.

## Optimized model line

The final gray line is the mean of the best TASEP SSA simulation selected from
15 candidate elongation rates between 1.0 and 8.0 aa/s in 0.5-aa/s steps.

For each candidate, the notebook compares the model and experimental mean
trajectories at the 25 time points from 0 through 24 minutes after
harringtonine application, with `ki` fixed at 0.03 1/s. It selects the
candidate with the lowest weighted residual cost:

\[
C(k_e)=\sum_{t=0}^{24}
\frac{[\bar I_{\mathrm{experimental}}(t)-\bar I_{\mathrm{model}}(t;k_e)]^2}
{[\mathrm{population\ SD}_{\mathrm{experimental}}(t)]^2+10^{-6}}.
\]

This is a finite candidate search, not a nonlinear fit with convergence
states. There is no failed-fit exclusion rule: a panel would return no model
only if no valid candidate cost could be calculated, which does not occur for
the four current Figure 7 panels.

Current settings are:

| Setting | Value |
|---|---:|
| Fixed initiation rate `ki` | 0.03 1/s |
| Candidate elongation rates | 1.0 to 8.0 aa/s, 15 values |
| Inhibition effectiveness | 95% |
| Simulation repetitions per candidate | 100 |
| Random seed | 42, deterministically derived per condition/rate/repetition |
| Burn-in time | 2,000 s |
| Simulation time step | 1 s |
| Displayed/downsampled interval | 60 s |
| HT application time before shifting | 5 min |
| Cost comparison window | 0 through 24 min after HT application |

The best current elongation rates from the notebook's stored execution are:

| Panel | Condition | Best `ke` (aa/s) |
|---|---|---:|
| Figure 7J | `UTag` | 3.5 |
| Figure 7K | `UTag_DeltaCys` | 3.0 |
| Figure 7L | `SunTag` | 4.0 |
| Figure 7M | `ALFA_tag` | 3.0 |

The final plot displays only the best model mean. Although model SD values are
calculated during optimization, no model uncertainty band is drawn in Figures
7J-M, so model SD should not be exported.

The best-fit rate, fixed initiation rate, seed, and model settings should be
documented in the README. The complete 15-rate optimization cost sweep is a
diagnostic analysis and is not proposed for the publication source-data
package.

## Proposed HT files and schemas

The same schema applies to:

- `Figure_7J_UTag_HT_curves.csv`
- `Figure_7K_UTag_DeltaCys_HT_curves.csv`
- `Figure_7L_SunTag_HT_curves.csv`
- `Figure_7M_ALFA_tag_HT_curves.csv`

Expected rows: 30 per file.

| Column | Definition |
|---|---|
| `time_min` | Time relative to harringtonine application in minutes; -5 through 24 |
| `mean_normalized_experimental_intensity` | Arithmetic mean normalized intensity across responding cells, displayed as the orange point/line series |
| `population_sd_normalized_experimental_intensity` | Symmetric population SD used for the orange shaded band |
| `n_responding_cells` | Number of responding cells contributing to the experimental mean and SD |
| `mean_model_fit_normalized_intensity` | Mean of the best-fit 100-repetition simulation displayed as the gray line |
| `n_model_repetitions` | Number of simulations contributing to the displayed model mean; currently 100 |

The cyan region from -1 to 0 minutes is an HT-application annotation and does
not require a numerical column. Its definition should be stated in the README.

### Mock Figure 7J-M HT table

Illustrative values only. Each panel will have its own CSV.

| time_min | mean_normalized_experimental_intensity | population_sd_normalized_experimental_intensity | n_responding_cells | mean_model_fit_normalized_intensity | n_model_repetitions |
|---:|---:|---:|---:|---:|---:|
| -5 | 0.982 | 0.148 | 32 | 1.000 | 100 |
| -4 | 1.011 | 0.132 | 32 | 0.995 | 100 |
| 0 | 0.947 | 0.401 | 32 | 0.962 | 100 |
| 1 | 0.891 | 0.388 | 32 | 0.901 | 100 |

---

# Proposed file manifest

| File | Expected rows | Plotted values represented |
|---|---:|---|
| `Figure_7C_spot_intensity_individual_cell_means.csv` | 60 | Every black dot in Figure 7C |
| `Figure_7D_spot_SNR_individual_cell_means.csv` | 60 | Every black dot in Figure 7D |
| `Figure_7E_spot_size_individual_cell_means.csv` | 60 | Every black dot in Figure 7E |
| `Figure_7F_UTag_ACF_curves.csv` | 295 | Visible experimental mean/band and model mean/band in Figure 7F |
| `Figure_7G_UTag_DeltaCys_ACF_curves.csv` | 295 | Visible experimental mean/band and model mean/band in Figure 7G |
| `Figure_7H_SunTag_ACF_curves.csv` | 295 | Visible experimental mean/band and model mean/band in Figure 7H |
| `Figure_7I_ALFA_tag_ACF_curves.csv` | 295 | Visible experimental mean/band and model mean/band in Figure 7I |
| `Figure_7J_UTag_HT_curves.csv` | 30 | Experimental mean/SD band and model line in Figure 7J |
| `Figure_7K_UTag_DeltaCys_HT_curves.csv` | 30 | Experimental mean/SD band and model line in Figure 7K |
| `Figure_7L_SunTag_HT_curves.csv` | 30 | Experimental mean/SD band and model line in Figure 7L |
| `Figure_7M_ALFA_tag_HT_curves.csv` | 30 | Experimental mean/SD band and model line in Figure 7M |
| `README_source_data.md` | n/a | Data dictionary, transformations, filters, fits, labels, and exclusions |

---

# Explicit exclusions

Do not export:

- panels 7A-B or image provenance;
- raw microscopy images or LIF metadata;
- source tracking filenames or folder names;
- experiment IDs or acquisition dates;
- cell IDs for the unpaired 7C-E boxplots;
- raw per-spot or per-frame rows underlying the per-cell means in 7C-E;
- box medians, quartiles, whiskers, or swarm jitter positions;
- p-values, Mann-Whitney statistics, significance bars, or stars;
- individual ACF trajectory correlations;
- individual TASEP simulation trajectories;
- hidden ACF lag coordinates outside the final 0-24.5-minute x-axis range;
- individual HT cell trajectories, because they are not displayed;
- nonresponding HT trajectories as a separate publication table;
- HT linear-fit diagnostics produced by the earlier non-final plotting pass;
- the 15-rate HT cost-sweep tables from `raw_data/` or `figures/`; retain them
  only as diagnostics in `optimization_tables/`;
- ACF cost surfaces or optimized-parameter search grids;
- plotting colors, marker sizes, axis coordinates, or legend annotations;
- a combined Figure 7 canonical dataset.

---

# Recommended implementation strategy

## Shared setup

1. Record hashes of all three notebooks, source/processed files, and archived
   final PNGs before changes.
2. Resolve `main_dir` and `notebooks/Figure_7` from the repository rather than
   from the Jupyter launch directory.
3. Create `figures/`, `raw_data/`, and `optimization_tables/` automatically
   with `pathlib`.
4. Set the notebook kernel metadata to the installed `microlive` kernel.
5. Preserve all source data as read-only inputs.
6. Preserve current labels, plot dimensions, styles, statistics, and final
   panel order unless the user separately approves a scientific change.

## Figures 7C-E

1. Resolve the mounted LaCie source root first; use the checksum-verified
   submitted-paper OneDrive copy only if the primary root is unavailable.
2. Select one source root for the entire run and print which root was selected.
3. Keep the current source traversal order for arrays passed to the plot; use
   a separate sorted inventory only for validation or checksum reporting.
4. Verify exactly one tracking CSV per included `results_*` folder.
5. Keep the existing missing-value removal, arithmetic mean, `show_stats`
   settings, and spot-size conversion unchanged.
6. Construct each panel's two-column export dataframe from the exact per-cell
   mean array passed to the box/swarm function.
7. Save each CSV immediately before or after its corresponding plot from that
   same dataframe with `index=False`.
8. Save only the three reviewed panel images to `figures/`; do not write new
   publication outputs to either source-data root.

## Figures 7F-I

1. Continue reading the four local `optimization/results_ACF/df_ACF_*.csv`
   tables used by the current final notebook.
2. Do not rerun or substitute the upstream experimental ACF calculation during
   source-data export.
3. Set the single configured model seed to
   `REPRODUCIBLE_RANDOM_SEED = 42`. Derive deterministic worker seeds for the
   parallel SSA repetitions from that base seed using the HT notebook's stable
   seed strategy.
4. Do not use the installed unseeded `microlive.Correlation` bootstrap for the
   model band. Use a notebook-local deterministic 1,000-resample bootstrap
   driven by the same base seed in a fixed order, while preserving the current
   resampling, baseline-correction, `np.nanstd`, and lag-zero behavior.
5. Build each model curve once. Use the exact returned mean, bootstrap error,
   and lag arrays for both plotting and export.
6. Keep the current full 360 experimental and 400 model coordinates when
   plotting, including the existing axis clipping behavior.
7. Convert lags to minutes and create a separate 295-row, 0-24.5-minute view
   only for the CSV export.
8. Build one panel-specific wide dataframe containing the exact experimental
   and model mean/error values used by the corresponding plot.
9. Save each CSV with `index=False` and save each PNG/SVG to `figures/`.
10. Treat the single seeded execution as the publication source. Do not compare
    the gray model curve with a prior stochastic execution.

## Figures 7J-M

1. Replace the obsolete missing input path with the required LaCie root
   `/Volumes/LaCie/UTag_paper_data/Harringtonine`; do not use OneDrive for the
   final J-M run.
2. Map the four notebook condition names explicitly to the LaCie directories
   `UTag`, `UTag_CF`, `SunTag`, and `AlfaTag`, each beneath its
   `HT_Analysis_GUI` directory.
3. Select only regular `.csv` files whose names begin with `tracking_`, reject
   every name beginning with `._`, and require exactly one real tracking CSV in
   each included `results_*` directory.
4. Preserve the current intensity normalization, 30% responder threshold,
   final-five-frame classification, and responder-only summary.
5. Preserve `use_sem=False`, so the experimental band remains population SD.
6. Preserve the seeded TASEP wrapper, fixed `ki`, candidate `ke` grid, cost
   definition, and best-fit selection. Treat the current stored notebook
   output as the regeneration reference and retain the April output as
   provenance for the documented small model-only difference.
7. Construct each panel's 30-row dataframe from the exact `full_frames`,
   experimental mean, experimental SD, and selected best-model mean arrays used
   by the final optimization plot.
8. Do not rerun the winning simulation solely for export. Reuse the already
   generated `list_ssa_arrays[best_sim_idx]` curve.
9. Save each CSV with `index=False` and redirect the final panel PNG/SVG files
   from the erroneous Figure 6 path to `notebooks/Figure_7/figures/`.
10. Create `notebooks/Figure_7/optimization_tables/` automatically and write
    the four 15-rate `ke_optimization_analysis_*.csv` tables only there. Keep
    these diagnostic tables out of both `raw_data/` and `figures/`.

No new standalone plotting package or generalized framework is needed. The
smallest safe change is to create export dataframes beside the arrays already
used by each notebook's plotting code.

---

# Separate code-simplification backlog

Do not combine this cleanup with the source-data export implementation. First
complete the export and visual regression work; then make one small cleanup
change at a time and rerun the affected native-panel regression. None of the
items below changes a scientific calculation when removed as described.

| Notebook | Candidate | Risk and condition |
|---|---|---|
| 7C-E, cell 2 | Remove `read_lif_file_print_intensity`; it is never called and is unrelated to the plotted source data. | Low |
| 7C-E, cells 0-1 | Remove the inactive OneDrive assignment/comment after an explicit LaCie-first fallback resolver is in place. | Low |
| 7C-E, cells 4-6 | Remove unused raw-row counts, empty dataframe allocation, and unused `list_merged_data` concatenation. | Low |
| 7C-E, cells 4-6 | Remove dormant control-spot and per-frame paths only if this notebook will not be reused for those analyses. | Medium |
| 7C-E, cell 7 | Stop returning `p_values`, which callers discard. Retain the Mann-Whitney calculation if its console output remains the source for manually assembled 7E annotations. | Low |
| 7F-I, cells 3, 7, 9, 13-16 | Remove duplicate rcParam assignment, unused globals/arguments, empty cells, commented-out saves, discarded parameter dataframe, unused loop index, and debug-only prints. | Low |
| 7F-I | Do not remove `calculate_codon_elongation_rates`, `create_probe_vector`, or the inactive multi-tau settings during export work. | Preserve |
| 7J-M, cell 10 and cells 8-9 | Move the non-final diagnostic pass and its helpers to an archival diagnostic notebook after the final J-M reference issue is resolved. | Medium; not before regression |
| 7J-M, cells 1, 6, 8-9, 12 | Remove current-directory variable, unused sequence bindings, unused diagnostics variables, commented cost/interpolation code, unused model-error list, unused `sim_time`, and unused `results_folder` function argument. | Low after targeted regression |
| 7J-M | Preserve the seeded SSA wrapper, responder classification, and final optimization code during source-data export. The small expected April-to-current gray-model difference is already documented. | Preserve |

---

# README requirements

`raw_data/README_source_data.md` should document:

1. The final Figure 7C-M panel-to-notebook and panel-to-condition maps.
2. The complete eleven-file CSV manifest.
3. The canonical ASCII label mapping.
4. What one row represents in every file.
5. That 7C-E dots are arithmetic means within one analyzed cell tracking file.
6. The source fields for intensity, SNR, and spot size.
7. The spot-size conversion factor and unit.
8. The median, IQR, and 5th-to-95th-percentile boxplot definitions.
9. That significance results are excluded.
10. The ACF acquisition interval, displayed lag range, and autocorrelation
    meaning.
11. The experimental ACF trajectory filtering, SNR threshold, missing-frame
    handling, baseline correction, lag-zero correction, outlier threshold, and
    bootstrap definition.
12. Experimental ACF trajectory counts for each panel.
13. The ACF TASEP parameters, repetitions, burn-in, time grid, bootstrap, and
    the single configured model seed `REPRODUCIBLE_RANDOM_SEED = 42`, including
    that both SSA and bootstrap randomness are derived from that seed.
14. The HT per-frame particle/intensity aggregation and two-stage
    normalization.
15. The HT application frame, shifted time grid, response threshold, and
    final-five-frame classification.
16. Total and responding HT cell counts for each condition.
17. That the HT experimental center is an arithmetic mean and its band is
    population SD, not SEM.
18. The HT model seed, repetitions, fixed `ki`, candidate `ke` grid,
    inhibition effectiveness, winning `ke` for each panel, 0-to-24-minute
    weighted-residual cost definition, and the absence of failed-fit
    exclusions.
19. That the cyan -1-to-0-minute region is an annotation, not a measured
    variable.
20. That the exports use the same prepared arrays as the plotting functions.
21. That the files are processed figure source data rather than unprocessed
    microscopy data.
22. That LaCie is the required 7J-M source root. For 7C-E, document LaCie as the
    primary source and the submitted-paper OneDrive tree as the byte-identical
    fallback.
23. Every explicit exclusion listed in this plan.

---

# Validation checklist

## Source integrity

1. Confirm the LaCie primary source selected for 7C-E contains exactly 19, 14,
   16, and 11 included tracking files. If fallback is exercised, confirm the
   same counts under OneDrive and record that the checksum-verified copy was
   selected.
2. Confirm that the required LaCie HT root is selected and that its mapped
   `UTag`, `UTag_CF`, `SunTag`, and `AlfaTag` directories contain exactly 38,
   24, 33, and 32 real tracking CSVs after excluding all `._*` AppleDouble
   files. Confirm exactly one real tracking CSV per included `results_*`
   directory. Do not use the OneDrive tree for the final J-M run.
3. Hash all source tracking CSVs, local ACF summary CSVs, optimized parameter
   files, gene sequence files, and the active `tasep_models` and `microlive`
   source revisions before and after execution.
4. Confirm no source file contents or modification times change.
5. Confirm no other manuscript figure folder changes.

## File manifest and schemas

6. Confirm exactly eleven CSVs and one README are generated in `raw_data/`.
7. Confirm exactly eleven PNG/SVG pairs are generated in `figures/`.
   Confirm that exactly four `ke_optimization_analysis_*.csv` diagnostic tables
   are generated in `optimization_tables/` and none are written to `figures/`
   or `raw_data/`.
8. Confirm no publication output is written to the LaCie drive, OneDrive,
   `optimization/results_ACF`, `Figure_6`, or a launch-directory-dependent
   folder.
9. Confirm every CSV has the exact approved columns and column order.
10. Confirm every numerical column parses as numeric and retains full
    computational precision.
11. Confirm there is no index column such as `Unnamed: 0`.
12. Confirm all filenames and condition labels are ASCII-only.

## Figures 7C-E reconciliation

13. Confirm 60 rows per file and condition counts of 19, 14, 16, and 11.
14. Confirm one exported row for every displayed black dot.
15. Recalculate every per-cell mean directly from its source tracking CSV and
    reconcile at full precision.
16. For 7E, separately verify the 0.12989318982387477 micrometer-per-pixel
   conversion.
17. Recalculate plotted medians and box percentiles from the exported values
   and compare with the rendered panels.
18. Verify the native 7E notebook panel without bars; manual final-composite
    significance annotations are outside this notebook's regression target.

## Figures 7F-I reconciliation

19. Confirm exactly 295 unique `lag_time_min` rows from 0 through 24.5 minutes
    inclusive at 5-second intervals independently for each panel.
20. Confirm every experimental mean and error value exactly matches the local
    `df_ACF_*.csv` value used by the plot.
21. Confirm the error column is described as a bootstrap estimate of the
    standard error: the standard deviation across 1,000 resampled
    mean-correlation curves.
22. Confirm experimental trajectory counts against the exact source reports.
23. Confirm every model mean and model error value is the same in the CSV and
    plotted arrays from that execution.
24. Confirm the model SSA and 1,000-resample bootstrap both derive from the
    single configured seed 42 and that the unseeded installed bootstrap is not
    called. No historical gray-model comparison is required.
25. Confirm full arrays—not the 295-row export subset—are still passed to the
    plotting function.

## Figures 7J-M reconciliation

26. Confirm exactly 30 unique `time_min` rows from -5 through 24 in each file.
27. Confirm total/responding counts of 38/32, 24/22, 33/32, and 32/32.
28. Recalculate the responder classification from the 30% threshold and last
    five frames.
29. Recalculate every experimental mean and population SD from responding
    trajectories and reconcile at full precision.
30. Confirm the four selected best `ke` values and that each exported model
    curve is the already selected best curve.
31. Confirm that the current seed-42 run reproduces the stored May notebook
    best-fit selections and recorded costs. Document the measured small
    gray-model displacement from the April reference; do not describe the
    panels as pixel-identical to April.
32. Run the selected J-M configuration twice and confirm identical HT CSVs,
    model curves, and best-fit selections.

## Notebook and plot regression

33. Validate notebook JSON and parse every code cell.
34. Execute each notebook from top to bottom with the installed `microlive`
    kernel from a clean output location.
35. Compare the regenerated 7C-E and 7J-M PNGs with their appropriate reference
    panels. For 7F-I, validate the generated CSV-to-plot identity and fixed
    experimental inputs; do not compare the seeded gray model with a historical
    stochastic model curve.
36. Require 7C-D and the native no-bar 7E panel to remain scientifically and
    visually unchanged, apart from output locations. Preserve manual composite
    annotations outside the notebook.
37. For 7F-I, require the experimental curves and bands to remain unchanged.
    The single seed-42 model execution defines the gray model output and needs
    no historical model-curve comparison.
38. For J-M, require exact experimental curves and SD bands. The gray model
    curve is expected to show the documented small April-to-current
    displacement; do not claim April pixel identity.
39. Treat SVG timestamps and internal Matplotlib IDs as metadata differences,
    and use rendered PNG comparison as the visual source of truth.
40. Run `git diff --check`, inspect the targeted notebook diff, and review
    `git status --short` without touching unrelated user files.
