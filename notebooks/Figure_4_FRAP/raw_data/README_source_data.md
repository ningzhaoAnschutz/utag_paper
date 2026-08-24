# Figure 4 FRAP source data

These files contain the processed numerical data used to generate the graphs in
Figure 4A-H. They are figure source data, not unprocessed acquisition-level raw
microscopy data. Files are regenerated from the in-memory data assembled from
the experiment-level CSV files in `FRAP_quantification/`;
`combined_FRAP_data.csv` is not read or required.

## Condition labels

| CSV label | Figure label |
|---|---|
| `UTag` | UTag |
| `UTag_DeltaCys` | UTag(ΔCys), ΔCys |
| `SunTag` | SunTag |
| `ALFA_tag` | ALFA-tag, ALFA |
| `HA` | HA |

ASCII-only condition labels are used in every CSV to avoid encoding and naming
inconsistencies.

## Figure 4A-E: condition-specific FRAP trajectories

Each panel has an individual-trajectory file and a time-point summary file.
The individual files contain `condition`, anonymous `cell_id`, `time_s`, and
`normalized_mean_roi_intensity`. Each value is the mean fluorescence intensity
within one cell ROI at one time point after per-cell min-max normalization:
`(value - cell minimum) / (cell maximum - cell minimum)`.

The summary files contain `condition`, `time_s`,
`mean_normalized_intensity`, and `n_cells`. The green line is the arithmetic
mean across the available cells at each time point. Panels A-E display no error
bars or error bands, so no variability column is included.

| Panel | Condition | Individual rows | Cells |
|---|---|---:|---:|
| 4A | `UTag` | 7,000 | 50 |
| 4B | `UTag_DeltaCys` | 6,160 | 44 |
| 4C | `SunTag` | 7,980 | 57 |
| 4D | `ALFA_tag` | 7,980 | 57 |
| 4E | `HA` | 5,880 | 42 |

## Figure 4F: combined FRAP trajectories

`Figure_4F_FRAP_all_conditions_timepoint_summary.csv` contains `condition`,
`time_s`, `mean_normalized_intensity`,
`sample_sd_normalized_intensity`, and `n_cells`. The colored line is the
arithmetic mean across cells and the shaded region is mean ± sample standard
deviation (SD; pandas default `ddof=1`). It is not SEM.

The existing Figure 4F plotting call repeats the recovery-drop quality check
with a threshold of 0.8 on the prepared normalized trajectories. It retains
49 UTag, 44 UTag_DeltaCys, 57 SunTag, 57 ALFA_tag, and 42 HA cells. This is why
the UTag `n_cells` in panel F is 49 rather than the 50 used in panel A.

## Figure 4G: recovered intensity

`Figure_4G_FRAP_recovered_intensity_individual_values.csv` contains
`condition` and `normalized_recovered_intensity`. Each row is one black dot:
the final recorded normalized mean ROI intensity for one cell. Box medians,
quartiles, whiskers, statistical tests, p-values, and significance labels are
not included.

## Figure 4H: recovery half-time

`Figure_4H_FRAP_t_half_individual_values.csv` contains `condition` and
`t_half_s`. Each row is one successfully fitted black dot. Fits use only the
post-bleach interval, with the 10-second bleach point reset to t=0. The model is
`I(t) = 1 - exp(-k*t)`, with the recovery plateau fixed at 1, and
`t_half_s = ln(2) / k`. Box summaries and statistical-test results are not
included.

## Data preparation

The initial notebook quality check retains trajectories whose fluorescence drop
during the first 20 seconds is greater than 0.5. The files contain 250 cells
from 15 experiment folders: 50 UTag from 3 folders, 44 UTag_DeltaCys from 3,
57 SunTag from 4, 57 ALFA_tag from 2, and 42 HA from 3. Figure 4F has the
additional filter described above.

There are 140 time points per complete trajectory: 0-39 seconds in one-second
increments, followed by 44-539 seconds in five-second increments. Optional
downsampling is disabled. The final values in Figure 4G are at 539 seconds.

Normalization divisions are guarded: a constant trajectory is assigned zero
after normalization. The CSV exports are created from the same prepared values
used by their corresponding plotting functions. The package contains no
p-values, significance tests, box summaries, image-provenance fields, source
image names, experiment identifiers, or replacement combined trajectory file.
