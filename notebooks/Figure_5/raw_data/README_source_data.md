# Figure 5C-E source data

These files contain the processed numerical values displayed in Figure 5C,
Figure 5D, and Figure 5E. They are figure source data, not acquisition-level raw
western-blot images or raw circular-dichroism instrument files.

## Source notebooks

| Panel | Notebook | Measurement |
|---|---|---|
| Figure 5C | `Figure_5_C.ipynb` | Soluble-protein ratio after incubation at 4, 50, 60, or 70 °C |
| Figure 5D | `Figure_5_D_E_TM.ipynb` | Wild-type UTag intrabody ellipticity at 216 nm |
| Figure 5E | `Figure_5_D_E_TM.ipynb` | UTag(ΔCys) intrabody ellipticity at 216 nm |

The notebooks read the approved source workbooks from the current manuscript
folder in OneDrive and write publication outputs only under
`notebooks/Figure_5/figures/` and `notebooks/Figure_5/raw_data/`.

## File manifest

| File | Rows | Displayed component |
|---|---:|---|
| `Figure_5C_thermal_solubility_summary.csv` | 12 | Figure 5C bar heights and symmetric error magnitudes |
| `Figure_5D_Anti_UTag_IB_CD_observed_values.csv` | 273 | Figure 5D experimental dots |
| `Figure_5D_Anti_UTag_IB_CD_fitted_curves.csv` | 600 | Figure 5D fitted lines, 200 coordinates per repetition |
| `Figure_5E_Anti_UTag_IB_DeltaCys_CD_observed_values.csv` | 273 | Figure 5E experimental dots |
| `Figure_5E_Anti_UTag_IB_DeltaCys_CD_fitted_curves.csv` | 600 | Figure 5E fitted lines, 200 coordinates per repetition |

## Canonical condition labels

CSV contents and filenames use ASCII-only condition labels:

| CSV label | Displayed figure label | Meaning |
|---|---|---|
| `wt_anti_UTag_scFv` | wt anti-UTag-scFv | Wild-type anti-UTag scFv control |
| `Anti_UTag_IB` | Anti-UTag-IB | Wild-type UTag intrabody |
| `Anti_UTag_IB_DeltaCys` | Anti-UTag-IB(ΔCys) | Cysteine-free UTag intrabody |

The source workbooks and older plot code sometimes use `FB`. Publication CSV
files use the approved term `IB`.

## Figure 5C summary

One row in `Figure_5C_thermal_solubility_summary.csv` represents one displayed
bar and its symmetric error magnitude.

| Column | Definition |
|---|---|
| `condition` | Canonical protein condition |
| `temperature_c` | Incubation temperature in degrees Celsius |
| `mean_normalized_soluble_protein_ratio` | Arithmetic mean plotted as the bar height |
| `population_sd_div_sqrt_n_normalized_soluble_protein_ratio` | Existing plotted error magnitude |
| `n_replicates` | Number of replicate ratios represented by the bar |

The soluble-protein value at each temperature \(T_i\) is normalized to the 4 °C
reference \(T_0\). The existing workbook calculates the error magnitude as:

\[
\frac{\operatorname{STDEV.P}(\text{replicate ratios})}{\sqrt{n}}
\]

This is population standard deviation divided by the square root of the sample
count. It is preserved to reproduce Figure 5C exactly and must not be described
as conventional sample-based SEM. Each displayed bar represents two
replicates.

Only the 12 displayed summary rows are exported. The 24 upstream individual
replicate values and the available but undisplayed 80 °C measurement are
intentionally excluded.

## Figures 5D and 5E observed values

One row in an `observed_values.csv` file represents one displayed experimental
dot.

| Column | Definition |
|---|---|
| `condition` | Canonical intrabody condition |
| `replicate_id` | Anonymous repetition identifier: `replicate_1`, `replicate_2`, or `replicate_3` |
| `temperature_c` | Experimental temperature in degrees Celsius |
| `ellipticity_at_216_nm` | Experimental circular-dichroism ellipticity at 216 nm |

Each panel contains three repetitions with 91 measurements per repetition.
The source workbook and notebook do not specify a more detailed physical unit
for ellipticity, so no unit is invented in the CSV.

## Figures 5D and 5E fitted curves

One row in a `fitted_curves.csv` file represents one coordinate of a displayed
fitted line.

| Column | Definition |
|---|---|
| `condition` | Canonical intrabody condition |
| `replicate_id` | Repetition fitted independently |
| `temperature_c` | Temperature coordinate of the displayed fitted line |
| `fitted_ellipticity_at_216_nm` | Fitted ellipticity displayed at that temperature |

Each repetition is fitted independently using:

\[
y(T) = d + \frac{a}{1 + \exp[-b(T-c)]}
\]

Here, \(a\) is the amplitude, \(b\) is the slope coefficient, \(c\) is the
midpoint or melting temperature, and \(d\) is the lower asymptote. Initial
values use the first percentile for \(d\), the difference between the 99th and
first percentiles for \(a\), 1.0 for \(b\), and 60 °C for \(c\). Fits are
unbounded and allow up to 20,000 function evaluations.

Each displayed line contains 200 temperatures spanning that repetition's
observed minimum and maximum. The repetitions have slightly different
temperature ranges, so the fitted values are stored in long format instead of
using one potentially inaccurate shared temperature column. A failed fit is
reported by the notebook and contributes no fitted-line coordinates.

## Export integrity and exclusions

The CSV files are written from the exact dataframe or arrays supplied to the
plotting calls. Numerical values are saved at available computational precision
without publication rounding.

The package intentionally excludes Figure 5C individual replicate values,
separate fit-parameter tables, manually entered rounded melting-temperature
summaries, p-values, statistical-test tables, image provenance, raw western
blot intensities, and raw instrument files.
