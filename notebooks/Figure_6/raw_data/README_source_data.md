# Figure 6B source data

This folder contains the processed numerical values displayed in Figure 6B.
These are figure source data, not acquisition-level raw microscopy images or
unprocessed fluorescence measurements.

## Source notebook

| Panel | Notebook | Measurement |
|---|---|---|
| Figure 6B | `Figure_6B.ipynb` | Normalized cell-level GFP/mCherry fluorescence-intensity ratio |

The notebook reads the existing source workbooks from the submitted-paper
OneDrive folder and writes publication outputs only under
`notebooks/Figure_6/figures/` and `notebooks/Figure_6/raw_data/`.

## File manifest

| File | Rows | Displayed component |
|---|---:|---|
| `Figure_6B_linker_normalized_GFP_mCh_ratio_individual_values.csv` | 306 | All black swarm points in Figure 6B |

## Columns

One row represents one analyzed cell and one displayed black dot.

| Column | Definition |
|---|---|
| `condition` | Canonical ASCII construct label |
| `normalized_gfp_to_mch_fluorescence_intensity_ratio` | Cell-level normalized GFP/mCherry fluorescence-intensity ratio |

## Conditions and sample counts

| CSV label | Displayed figure label | Cells |
|---|---|---:|
| `GFP_mCh` | GFP-mCh | 128 |
| `5aa_linker` | 5aa linker | 62 |
| `9aa_linker` | 9aa linker | 62 |
| `13aa_linker` | 13aa linker | 54 |

Cells are pooled across the approved source workbooks for each construct.

## Measurement and normalization

Each source workbook calculates:

\[
\text{Ratio} = \frac{\text{GFP intensity}}{\text{mCherry intensity}}
\]

\[
\text{Normalized ratio} =
\frac{\text{Ratio}}{\text{acquisition-specific normalization factor}}
\]

The source-workbook factors are:

| Acquisition date | Normalization factor |
|---|---:|
| 2024-06-14 | 1.10 |
| 2024-07-03 | 2.90 |
| 2024-07-17 | 0.82 |

The notebook does not recalculate the normalization. It reads the exact
precomputed `Normalized ratio` values and supplies them to both the CSV export
and the plot.

## Plot definition

Figure 6B displays:

- one black swarm point per cell;
- a red median line;
- a box spanning the 25th to 75th percentiles;
- whiskers defined using the 5th and 95th percentiles;
- all individual points, including values beyond the whiskers;
- significance bars for significant two-sided Mann–Whitney U comparisons.

Only the individual black-dot values are included in the CSV.

## Source workbook mapping

| Condition | Source workbook |
|---|---|
| `GFP_mCh` | `20240703 pNZ257.xlsx` |
| `GFP_mCh` | `20240614_pNZ257.xlsx` |
| `GFP_mCh` | `20240717 pNZ257.xlsx` |
| `5aa_linker` | `20240717 pNZ112_263.xlsx` |
| `5aa_linker` | `20240703 pNZ112_263.xlsx` |
| `9aa_linker` | `20240717 pNZ112_360.xlsx` |
| `9aa_linker` | `20240703 pNZ112_360.xlsx` |
| `13aa_linker` | `20240717 pNZ112_361.xlsx` |
| `13aa_linker` | `20240614_pNZ112_361.xlsx` |

## Exclusions

The source-data package intentionally excludes box geometry, printed means and
standard deviations, Mann–Whitney statistics, p-values, significance stars,
swarm jitter coordinates, experiment IDs, source filenames as CSV columns,
source cell labels, raw GFP and mCherry intensities, unnormalized ratios, and
image provenance.
