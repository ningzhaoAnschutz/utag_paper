"""Processing and plotting helpers for spot-distribution figure notebooks."""

from itertools import combinations
from pathlib import Path
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import gaussian_kde, mannwhitneyu


__all__ = [
    "build_audit_dataframe",
    "build_size_filter_audit",
    "extract_cell_arrays",
    "extract_particle_means_by_cell",
    "extract_particle_distribution",
    "extract_spots_per_cell",
    "load_tracking_records",
    "plot_cell_means_box_swarm",
    "plot_distributions",
    "summarize_spots_above_threshold",
    "run_distribution_analyses",
    "run_large_spot_analysis",
    "run_whisker_analyses",
    "save_pairwise_statistics",
]


_SIZE_FILTER_MODES = frozenset({"trajectory_mean", "detection"})


def _real_tracking_files(dataset_root: Path):
    """Return one tracking CSV from each results directory."""
    if not dataset_root.is_dir():
        raise NotADirectoryError(f"Dataset directory not found: {dataset_root}")

    result_dirs = sorted(
        path for path in dataset_root.iterdir()
        if path.is_dir() and "results_" in path.name
    )
    if not result_dirs:
        raise ValueError(f"No results directories found in {dataset_root}")

    tracking_files = []
    for result_dir in result_dirs:
        candidates = sorted(
            path for path in result_dir.iterdir()
            if path.is_file()
            and path.suffix.lower() == ".csv"
            and path.name.startswith("tracking_")
            and not path.name.startswith("._")
        )
        if len(candidates) != 1:
            raise ValueError(
                f"Expected one tracking CSV in {result_dir}; found {len(candidates)}"
            )
        tracking_files.append(candidates[0])
    return tracking_files


def _metadata_pixel_size_um(tracking_file: Path):
    """Read voxel_yx_nm from the metadata beside a tracking CSV."""
    metadata_files = sorted(
        path for path in tracking_file.parent.iterdir()
        if path.is_file()
        and path.name.startswith("Metadata_")
        and path.suffix.lower() == ".txt"
        and not path.name.startswith("._")
    )
    if len(metadata_files) != 1:
        raise ValueError(
            f"Expected one Metadata_*.txt beside {tracking_file.name}; "
            f"found {len(metadata_files)}"
        )

    metadata_text = metadata_files[0].read_text(errors="replace")
    match = re.search(
        r"^\s*voxel_yx_nm:\s*([0-9.eE+-]+)\s*$",
        metadata_text,
        re.MULTILINE,
    )
    if match is None:
        raise ValueError(f"voxel_yx_nm not found in {metadata_files[0]}")
    return float(match.group(1)) / 1000.0


def _particle_column(dataframe):
    """Return the available tracked-particle identifier column."""
    return "unique_particle" if "unique_particle" in dataframe.columns else "particle"


def _validate_size_filter_mode(size_filter_mode):
    """Return a validated size-filter mode."""
    if size_filter_mode not in _SIZE_FILTER_MODES:
        allowed = ", ".join(sorted(_SIZE_FILTER_MODES))
        raise ValueError(
            f"size_filter_mode must be one of {allowed}; "
            f"received {size_filter_mode!r}"
        )
    return size_filter_mode


def _prepare_size_filter_table(
    dataframe,
    *,
    size_column,
    pixel_size_um,
    threshold_um,
    size_filter_mode,
    snr_column=None,
    min_snr=None,
    selected_frame=None,
):
    """Annotate valid detections and apply the selected size-filter unit."""
    size_filter_mode = _validate_size_filter_mode(size_filter_mode)
    if threshold_um is not None and (
        not np.isfinite(threshold_um) or threshold_um <= 0
    ):
        raise ValueError("threshold_um must be positive and finite")
    if not np.isfinite(pixel_size_um) or pixel_size_um <= 0:
        raise ValueError("pixel_size_um must be positive and finite")

    prepared = dataframe.copy()
    if selected_frame is not None and "frame" in prepared.columns:
        prepared = prepared.loc[prepared["frame"] == selected_frame].copy()

    particle_column = _particle_column(prepared)
    if particle_column not in prepared.columns:
        raise KeyError("Tracked-particle identifier is missing")

    valid_mask = prepared[particle_column].notna().to_numpy()
    if snr_column is not None and snr_column in prepared.columns:
        snr_values = pd.to_numeric(
            prepared[snr_column], errors="coerce"
        ).to_numpy(dtype=float)
        valid_mask &= np.isfinite(snr_values)
        if min_snr is not None:
            valid_mask &= snr_values >= min_snr

    if threshold_um is None:
        prepared = prepared.loc[valid_mask].copy()
        prepared["_spot_size_um"] = np.nan
        prepared["_trajectory_mean_spot_size_um"] = np.nan
        prepared["_detection_above_size_threshold"] = False
        prepared["_trajectory_above_size_threshold"] = False
        prepared["_retained_by_size_filter"] = True
        return prepared

    if size_column not in prepared.columns:
        raise KeyError(f"{size_column} is missing")
    size_pixels = pd.to_numeric(
        prepared[size_column], errors="coerce"
    ).to_numpy(dtype=float)
    size_um = size_pixels * float(pixel_size_um)
    valid_mask &= np.isfinite(size_um) & (size_um > 0)
    prepared = prepared.loc[valid_mask].copy()
    prepared["_spot_size_um"] = size_um[valid_mask]

    if prepared.empty:
        prepared["_trajectory_mean_spot_size_um"] = pd.Series(dtype=float)
        prepared["_detection_above_size_threshold"] = pd.Series(dtype=bool)
        prepared["_trajectory_above_size_threshold"] = pd.Series(dtype=bool)
        prepared["_retained_by_size_filter"] = pd.Series(dtype=bool)
        return prepared

    prepared["_trajectory_mean_spot_size_um"] = prepared.groupby(
        particle_column, sort=False
    )["_spot_size_um"].transform("mean")
    prepared["_detection_above_size_threshold"] = (
        prepared["_spot_size_um"] > threshold_um
    )
    prepared["_trajectory_above_size_threshold"] = (
        prepared["_trajectory_mean_spot_size_um"] > threshold_um
    )
    if size_filter_mode == "trajectory_mean":
        retained = ~prepared["_trajectory_above_size_threshold"]
    else:
        retained = ~prepared["_detection_above_size_threshold"]
    prepared["_retained_by_size_filter"] = retained
    return prepared


def _size_filter_label(size_filter_mode, threshold_um):
    """Return a concise label describing the active size filter."""
    _validate_size_filter_mode(size_filter_mode)
    if size_filter_mode == "trajectory_mean":
        return rf"trajectory mean $\leq {threshold_um:g}\ \mu$m"
    return rf"detections $\leq {threshold_um:g}\ \mu$m"


def _size_filter_count_label(size_filter_mode, threshold_um):
    """Return the retention rule used by a filtered trajectory count."""
    _validate_size_filter_mode(size_filter_mode)
    if size_filter_mode == "trajectory_mean":
        return rf"mean size $\leq {threshold_um:g}\ \mu$m"
    return rf"$\geq 1$ detection $\leq {threshold_um:g}\ \mu$m"


def load_tracking_records(
    input_root,
    dataset_folders,
    use_metadata_pixel_size=False,
    fixed_pixel_xy_um=0.12989318982387477,
):
    """Load tracking tables and lateral pixel calibration by condition."""
    records = {key: [] for key in dataset_folders}
    for key, folder_name in dataset_folders.items():
        for tracking_file in _real_tracking_files(Path(input_root) / folder_name):
            pixel_size_um = (
                _metadata_pixel_size_um(tracking_file)
                if use_metadata_pixel_size
                else fixed_pixel_xy_um
            )
            records[key].append({
                "tracking_file": tracking_file,
                "cell_name": tracking_file.parent.name,
                "dataframe": pd.read_csv(tracking_file),
                "pixel_size_um": pixel_size_um,
            })
    return records


def build_audit_dataframe(tracking_records, condition_labels):
    """Summarize loaded cells, rows, particles, and pixel calibrations."""
    rows = []
    for key, records in tracking_records.items():
        rows.append({
            "condition": condition_labels[key],
            "cells/files": len(records),
            "dataframe rows": sum(len(record["dataframe"]) for record in records),
            "unique particles": sum(
                record["dataframe"]["particle"].nunique() for record in records
            ),
            "pixel sizes (um/px)": sorted({
                round(record["pixel_size_um"], 12) for record in records
            }),
        })
    return pd.DataFrame(rows)


def _save_figure(fig, output_dir, stem, save_formats, save_dpi):
    """Save a figure in each configured output format."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    for extension in save_formats:
        output_path = output_dir / f"{stem}.{extension}"
        fig.savefig(output_path, dpi=save_dpi, bbox_inches="tight")
        paths.append(output_path)
    print("Saved:", " | ".join(str(path) for path in paths))
    return paths


def extract_cell_arrays(
    tracking_records,
    condition_keys,
    selected_field,
    selected_frame=None,
    convert_pixels_to_um=False,
    max_spot_size_um=None,
    channel_index=0,
):
    """Return one array per cell for every condition, optionally filtered by spot size."""
    condition_arrays = []
    size_field = f"spot_size_ch_{channel_index}"
    for key in condition_keys:
        cell_arrays = []
        for record in tracking_records[key]:
            dataframe = record["dataframe"]
            if selected_field not in dataframe.columns:
                raise KeyError(f"{selected_field} missing from {record['tracking_file']}")
            if max_spot_size_um is not None and size_field not in dataframe.columns:
                raise KeyError(f"{size_field} missing from {record['tracking_file']}")

            sub_df = dataframe
            if selected_frame is not None and "frame" in sub_df.columns:
                sub_df = sub_df.loc[sub_df["frame"] == selected_frame]

            values = sub_df[selected_field].to_numpy(dtype=float)
            valid_mask = np.isfinite(values)

            if max_spot_size_um is not None:
                spot_size_um = (
                    sub_df[size_field].to_numpy(dtype=float)
                    * record["pixel_size_um"]
                )
                valid_mask = (
                    valid_mask
                    & np.isfinite(spot_size_um)
                    & (spot_size_um <= max_spot_size_um)
                )

            values = values[valid_mask]
            if convert_pixels_to_um:
                values = values * record["pixel_size_um"]
            cell_arrays.append(values)
        condition_arrays.append(cell_arrays)
    return condition_arrays


def extract_particle_means_by_cell(
    tracking_records,
    condition_keys,
    selected_field,
    max_particle_mean_size_um,
    selected_frame=None,
    convert_pixels_to_um=False,
    channel_index=0,
    min_snr=None,
    size_filter_mode="trajectory_mean",
):
    """Return equally weighted particle means per cell after a size filter.

    In ``trajectory_mean`` mode, an entire particle is retained when its mean
    size is at most the threshold. In ``detection`` mode, individual detections
    above the threshold are removed first. The returned arrays contain one mean
    selected-field value per remaining particle, preserving equal weighting.
    """
    size_filter_mode = _validate_size_filter_mode(size_filter_mode)
    condition_arrays = []
    size_field = f"spot_size_ch_{channel_index}"
    snr_field = f"snr_ch_{channel_index}"

    for key in condition_keys:
        cell_arrays = []
        for record in tracking_records[key]:
            dataframe = record["dataframe"]
            required_fields = {selected_field, size_field}
            missing_fields = required_fields.difference(dataframe.columns)
            if missing_fields:
                missing_text = ", ".join(sorted(missing_fields))
                raise KeyError(
                    f"{missing_text} missing from {record['tracking_file']}"
                )

            try:
                prepared = _prepare_size_filter_table(
                    dataframe,
                    size_column=size_field,
                    pixel_size_um=record["pixel_size_um"],
                    threshold_um=max_particle_mean_size_um,
                    size_filter_mode=size_filter_mode,
                    snr_column=snr_field,
                    min_snr=min_snr,
                    selected_frame=selected_frame,
                )
            except KeyError as exc:
                raise KeyError(
                    f"{exc.args[0]} in {record['tracking_file']}"
                ) from exc

            retained = prepared.loc[
                prepared["_retained_by_size_filter"]
            ].copy()
            particle_column = _particle_column(retained)
            selected_values = pd.to_numeric(
                retained[selected_field], errors="coerce"
            ).to_numpy(dtype=float)
            if convert_pixels_to_um:
                selected_values = selected_values * record["pixel_size_um"]

            valid_selected = np.isfinite(selected_values)
            particle_frame_values = pd.DataFrame({
                "particle": retained[particle_column].to_numpy()[valid_selected],
                "selected_value": selected_values[valid_selected],
            })
            eligible_values = (
                particle_frame_values.groupby("particle", sort=False)[
                    "selected_value"
                ]
                .mean()
                .to_numpy(dtype=float)
            )
            cell_arrays.append(eligible_values[np.isfinite(eligible_values)])
        condition_arrays.append(cell_arrays)
    return condition_arrays


def build_size_filter_audit(
    tracking_records,
    condition_keys,
    condition_labels,
    threshold_um,
    channel_index=0,
    min_snr=0.0,
    selected_frame=None,
    size_filter_mode="trajectory_mean",
):
    """Summarize detection and trajectory retention for every cell."""
    size_filter_mode = _validate_size_filter_mode(size_filter_mode)
    size_column = f"spot_size_ch_{channel_index}"
    snr_column = f"snr_ch_{channel_index}"
    rows = []

    for key in condition_keys:
        for record in tracking_records[key]:
            try:
                prepared = _prepare_size_filter_table(
                    record["dataframe"],
                    size_column=size_column,
                    pixel_size_um=record["pixel_size_um"],
                    threshold_um=threshold_um,
                    size_filter_mode=size_filter_mode,
                    snr_column=snr_column,
                    min_snr=min_snr,
                    selected_frame=selected_frame,
                )
            except KeyError as exc:
                raise KeyError(
                    f"{exc.args[0]} in {record['tracking_file']}"
                ) from exc

            particle_column = _particle_column(prepared)
            retained = prepared.loc[
                prepared["_retained_by_size_filter"]
            ]
            particles_before = int(prepared[particle_column].nunique())
            particles_retained = int(retained[particle_column].nunique())
            trajectory_flags = prepared.drop_duplicates(particle_column)
            rows.append({
                "condition_key": key,
                "condition": condition_labels[key],
                "cell_file": record["cell_name"],
                "tracking_csv": str(record["tracking_file"]),
                "pixel_size_um": record["pixel_size_um"],
                "threshold_um": threshold_um,
                "size_filter_mode": size_filter_mode,
                "valid_detections_before_filter": int(len(prepared)),
                "detections_above_threshold": int(
                    prepared["_detection_above_size_threshold"].sum()
                ),
                "detections_retained": int(len(retained)),
                "trajectories_before_filter": particles_before,
                "trajectories_with_mean_above_threshold": int(
                    trajectory_flags[
                        "_trajectory_above_size_threshold"
                    ].sum()
                ),
                "trajectories_retained": particles_retained,
                "trajectories_losing_all_detections": (
                    particles_before - particles_retained
                ),
            })
    return pd.DataFrame(rows)



def _significance_label(p_value):
    """Convert a p-value to its conventional significance label."""
    if not np.isfinite(p_value):
        return "na"
    if p_value < 0.0001:
        return "****"
    if p_value < 0.001:
        return "***"
    if p_value < 0.01:
        return "**"
    if p_value < 0.05:
        return "*"
    return "ns"


def plot_cell_means_box_swarm(
    conditions_data,
    condition_labels,
    output_dir,
    y_label,
    plot_stem,
    y_lim=None,
    show_stats=True,
    figsize=(6.5, 4.5),
    tick_size=18,
    max_percentile_significance=99.5,
    significance_top_padding=0.035,
    save_formats=("png", "svg"),
    save_dpi=600,
):
    """Plot one dot per cell using the mean of each supplied array."""
    rows = []
    for condition_label, repetitions in zip(condition_labels, conditions_data):
        for repetition in repetitions:
            values = np.asarray(repetition, dtype=float)
            rows.append({
                "Condition": condition_label,
                "Mean": np.nanmean(values) if values.size else np.nan,
            })
    plot_df = pd.DataFrame(rows).dropna(subset=["Mean"])

    fig, ax = plt.subplots(figsize=figsize, facecolor="white")
    sns.boxplot(
        data=plot_df,
        x="Condition",
        y="Mean",
        order=condition_labels,
        showfliers=False,
        boxprops={"facecolor": "white", "edgecolor": "black"},
        medianprops={"color": "red", "linewidth": 1.8},
        whiskerprops={"color": "black", "linewidth": 1.5},
        capprops={"color": "black", "linewidth": 1.5},
        linewidth=1.5,
        whis=[5, 95],
        width=0.5,
        ax=ax,
    )
    sns.swarmplot(
        data=plot_df,
        x="Condition",
        y="Mean",
        order=condition_labels,
        color="black",
        size=5,
        ax=ax,
        zorder=3,
    )

    ax.set_xlabel("")
    ax.set_ylabel(y_label, fontsize=tick_size + 4, color="black")
    ax.tick_params(
        axis="both",
        labelsize=tick_size,
        width=2,
        length=6,
        colors="black",
    )
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontname("Arial")
    for spine in ax.spines.values():
        spine.set_linewidth(2.0)
        spine.set_color("black")

    p_value_rows = []
    significant_pairs = []
    finite_means = plot_df["Mean"].to_numpy(dtype=float)
    global_min = np.nanmin(finite_means)
    global_max = np.nanpercentile(finite_means, max_percentile_significance)

    for i, j in combinations(range(len(condition_labels)), 2):
        group_1 = plot_df.loc[
            plot_df["Condition"] == condition_labels[i], "Mean"
        ].dropna()
        group_2 = plot_df.loc[
            plot_df["Condition"] == condition_labels[j], "Mean"
        ].dropna()
        p_value = (
            mannwhitneyu(group_1, group_2, alternative="two-sided").pvalue
            if len(group_1) and len(group_2)
            else np.nan
        )
        significance = _significance_label(p_value)
        p_value_rows.append({
            "comparison": f"{condition_labels[i]} vs {condition_labels[j]}",
            "p_value": p_value,
            "significance": significance,
        })
        if show_stats and significance not in {"ns", "na"}:
            significant_pairs.append((i, j, significance))

    significant_pairs = sorted(
        significant_pairs,
        key=lambda pair: (pair[1] - pair[0], pair[0]),
    )

    if y_lim is not None:
        lower, upper = map(float, y_lim)
        if not lower < upper:
            raise ValueError(f"y_lim must satisfy lower < upper; received {y_lim}")
        axis_span = upper - lower
        ax.set_ylim(lower, upper)

        if significant_pairs:
            visible_means = finite_means[
                (finite_means >= lower) & (finite_means <= upper)
            ]
            data_top = (
                np.nanpercentile(visible_means, max_percentile_significance)
                if visible_means.size
                else lower
            )
            cap_height = 0.008 * axis_span
            text_offset = 0.005 * axis_span
            lowest_line = data_top + 0.025 * axis_span
            highest_allowed = (
                upper
                - cap_height
                - text_offset
                - significance_top_padding * axis_span
            )
            number_of_bars = len(significant_pairs)
            level_spacing = 0.055 * axis_span
            if number_of_bars > 1:
                maximum_lowest_line = (
                    highest_allowed - level_spacing * (number_of_bars - 1)
                )
                lowest_line = min(lowest_line, maximum_lowest_line)
                minimum_lowest_line = lower + 0.10 * axis_span
                if lowest_line < minimum_lowest_line:
                    lowest_line = minimum_lowest_line
                    available = max(highest_allowed - lowest_line, 0)
                    level_spacing = available / (number_of_bars - 1)
            else:
                lowest_line = min(lowest_line, highest_allowed)
            bracket_levels = lowest_line + level_spacing * np.arange(number_of_bars)
        else:
            bracket_levels = []
    elif significant_pairs:
        data_range = global_max - global_min
        if not np.isfinite(data_range) or data_range <= 0:
            data_range = max(abs(global_max) * 0.1, 1.0)
        cap_height = 0.015 * data_range
        text_offset = 0.012 * data_range
        level_spacing = 0.10 * data_range
        first_level = global_max + 0.06 * data_range
        bracket_levels = first_level + level_spacing * np.arange(
            len(significant_pairs)
        )
        highest = bracket_levels[-1] + cap_height + text_offset
        lower = min(0, global_min) if global_min >= 0 else global_min
        ax.set_ylim(lower, highest + 0.10 * max(highest - lower, data_range))
    else:
        bracket_levels = []

    for (i, j, significance), y_line in zip(significant_pairs, bracket_levels):
        x_left = i + 0.04
        x_right = j - 0.04
        ax.plot(
            [x_left, x_left, x_right, x_right],
            [y_line, y_line + cap_height, y_line + cap_height, y_line],
            lw=1.25,
            color="black",
            solid_capstyle="round",
            solid_joinstyle="round",
            clip_on=True,
            zorder=4,
        )
        ax.text(
            (i + j) / 2,
            y_line + cap_height + text_offset,
            significance,
            ha="center",
            va="bottom",
            fontsize=max(tick_size - 5, 10),
            fontweight="normal",
            color="black",
            bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.25},
            clip_on=True,
            zorder=5,
        )

    if y_lim is not None:
        ax.set_ylim(y_lim)
    fig.tight_layout()
    _save_figure(fig, output_dir, plot_stem, save_formats, save_dpi)
    plt.show()
    plt.close(fig)

    summary_df = (
        plot_df.groupby("Condition", sort=False)["Mean"]
        .agg(["count", "mean", "std", "median"])
        .reset_index()
    )
    return plot_df, summary_df, pd.DataFrame(p_value_rows)


def extract_particle_distribution(
    tracking_records,
    condition_keys,
    selected_field_prefix,
    channel_index=0,
    min_snr=0.0,
    timepoint_frame=None,
    convert_pixels_to_um=False,
    max_spot_size_um=None,
    size_filter_mode="trajectory_mean",
):
    """Extract existing CSV values in three modes without pixel recalculation."""
    size_filter_mode = _validate_size_filter_mode(size_filter_mode)
    value_column = f"{selected_field_prefix}{channel_index}"
    snr_column = f"snr_ch_{channel_index}"
    size_column = f"spot_size_ch_{channel_index}"
    results = []

    for key in condition_keys:
        mean_per_particle = []
        all_timepoints = []
        at_timepoint = []
        valid_frames = []
        prepared_records = []

        for record in tracking_records[key]:
            dataframe = record["dataframe"].copy()
            if value_column not in dataframe.columns:
                raise KeyError(f"{value_column} missing from {record['tracking_file']}")
            try:
                prepared = _prepare_size_filter_table(
                    dataframe,
                    size_column=size_column,
                    pixel_size_um=record["pixel_size_um"],
                    threshold_um=max_spot_size_um,
                    size_filter_mode=size_filter_mode,
                    snr_column=snr_column,
                    min_snr=min_snr,
                )
            except KeyError as exc:
                raise KeyError(
                    f"{exc.args[0]} in {record['tracking_file']}"
                ) from exc

            prepared = prepared.loc[
                prepared["_retained_by_size_filter"]
            ].copy()
            selected_values = pd.to_numeric(
                prepared[value_column], errors="coerce"
            ).to_numpy(dtype=float)
            valid_selected = np.isfinite(selected_values) & (selected_values > 0)
            prepared = prepared.loc[valid_selected].copy()
            selected_values = selected_values[valid_selected]
            if convert_pixels_to_um:
                selected_values = selected_values * record["pixel_size_um"]
            prepared[value_column] = selected_values
            if prepared.empty:
                continue

            particle_column = _particle_column(prepared)
            means = prepared.groupby(particle_column)[value_column].mean()
            mean_per_particle.extend(means.to_list())
            all_timepoints.extend(prepared[value_column].to_list())
            if "frame" in prepared.columns:
                valid_frames.extend(prepared["frame"].to_list())
            prepared_records.append(prepared)

        selected_timepoint = (
            int(np.median(valid_frames))
            if timepoint_frame is None and valid_frames
            else timepoint_frame
        )
        if selected_timepoint is not None:
            for dataframe in prepared_records:
                if "frame" in dataframe.columns:
                    at_timepoint.extend(
                        dataframe.loc[
                            dataframe["frame"] == selected_timepoint,
                            value_column,
                        ].to_list()
                    )

        results.append({
            "condition_key": key,
            "mean_per_particle": np.asarray(mean_per_particle, dtype=float),
            "all_timepoints": np.asarray(all_timepoints, dtype=float),
            "at_timepoint": np.asarray(at_timepoint, dtype=float),
            "timepoint_frame_used": selected_timepoint,
        })
    return results


def _distribution_values(result, mode):
    """Return finite positive distribution values for the requested mode."""
    mode_map = {1: "mean_per_particle", 2: "all_timepoints", 3: "at_timepoint"}
    if mode not in mode_map:
        raise ValueError(f"mode must be 1, 2, or 3; received {mode}")
    values = np.asarray(result[mode_map[mode]], dtype=float)
    return values[np.isfinite(values) & (values > 0)]


def plot_distributions(
    distribution_results,
    condition_labels,
    condition_colors,
    output_dir,
    x_label,
    plot_stem,
    mode=1,
    xlim=None,
    figsize=(6.6, 5.1),
    show_inset=False,
    inset_xlim=(0.45, 1.0),
    inset_bounds=(0.43, 0.46, 0.50, 0.34),
    show_density_fill=True,
    legend_inside=False,
    show_tail_background=True,
    show_inset_fill=True,
    show_inset_title=True,
    inset_title_size=11.5,
    inset_tick_size=11.5,
    tick_size=19,
    spine_linewidth=2.0,
    bins=60,
    save_formats=("png", "svg"),
    save_dpi=600,
):
    """Plot KDE curves and return their observations and coordinates."""
    fig, ax = plt.subplots(figsize=figsize, facecolor="white")
    curves = []
    observation_rows = []
    curve_rows = []

    for result in distribution_results:
        key = result["condition_key"]
        values = _distribution_values(result, mode)
        if values.size == 0:
            continue

        observation_rows.extend({
            "condition_key": key,
            "value": value,
        } for value in values)

        color = condition_colors[key]
        label = condition_labels[key]
        if values.size > 3 and np.nanstd(values) > 0:
            x_values = np.linspace(values.min(), values.max(), 900)
            try:
                y_values = gaussian_kde(values, bw_method="scott")(x_values)
                if show_density_fill:
                    ax.fill_between(
                        x_values,
                        0,
                        y_values,
                        color=color,
                        alpha=0.075,
                        linewidth=0,
                        zorder=1,
                    )
                ax.plot(
                    x_values,
                    y_values,
                    color=color,
                    linewidth=2.5,
                    solid_capstyle="round",
                    label=label,
                    zorder=3,
                )

                median_value = float(np.median(values))
                median_density = float(np.interp(median_value, x_values, y_values))
                ax.vlines(
                    median_value,
                    0,
                    median_density,
                    color=color,
                    linewidth=1.15,
                    linestyles=(0, (3, 2)),
                    alpha=0.62,
                    zorder=2,
                )
                curves.append((x_values, y_values, color))
                curve_rows.extend({
                    "condition_key": key,
                    "x_value": x_value,
                    "probability_density": y_value,
                    "n_observations": int(values.size),
                } for x_value, y_value in zip(x_values, y_values))
            except np.linalg.LinAlgError:
                ax.hist(
                    values,
                    bins=bins,
                    density=True,
                    histtype="step",
                    linewidth=2.0,
                    color=color,
                    label=label,
                )
        else:
            ax.hist(
                values,
                bins=bins,
                density=True,
                histtype="step",
                linewidth=2.0,
                color=color,
                label=label,
            )

    ax.set_xlabel(x_label, fontsize=22, labelpad=6, color="black")
    ax.set_ylabel(
        "Probability Density", fontsize=22, labelpad=6, color="black"
    )
    ax.tick_params(
        axis="both",
        which="major",
        labelsize=tick_size,
        width=spine_linewidth,
        length=6,
        direction="out",
        colors="black",
    )
    ax.set_ylim(bottom=0)
    if xlim is not None:
        ax.set_xlim(xlim)
    for spine in ax.spines.values():
        spine.set_linewidth(spine_linewidth)
        spine.set_color("black")

    if legend_inside:
        ax.legend(
            loc="lower right",
            bbox_to_anchor=(0.985, 0.08),
            ncol=1,
            frameon=False,
            fontsize=13.0,
            handlelength=2.25,
            handletextpad=0.65,
            labelspacing=0.42,
            borderpad=0,
            borderaxespad=0,
            labelcolor="black",
        )
    else:
        ax.legend(
            loc="lower center",
            bbox_to_anchor=(0.5, 1.025),
            ncol=len(condition_labels),
            frameon=False,
            fontsize=14.5,
            handlelength=2.5,
            handletextpad=0.7,
            columnspacing=1.45,
            borderpad=0,
            borderaxespad=0,
            labelcolor="black",
        )

    if show_inset and curves:
        if show_tail_background:
            ax.axvspan(
                inset_xlim[0],
                inset_xlim[1],
                color="#6B7280",
                alpha=0.025,
                linewidth=0,
                zorder=0,
            )
        inset_ax = ax.inset_axes(inset_bounds)
        inset_ax.set_facecolor("white")
        inset_y_values = []
        for x_values, y_values, color in curves:
            if show_inset_fill:
                inset_ax.fill_between(
                    x_values,
                    0,
                    y_values,
                    color=color,
                    alpha=0.055,
                    linewidth=0,
                )
            inset_ax.plot(
                x_values,
                y_values,
                color=color,
                linewidth=1.9,
                solid_capstyle="round",
            )
            region = (x_values >= inset_xlim[0]) & (x_values <= inset_xlim[1])
            if np.any(region):
                inset_y_values.extend(y_values[region].tolist())

        inset_ax.set_xlim(inset_xlim)
        if inset_y_values:
            y_min = min(inset_y_values)
            y_max = max(inset_y_values)
            padding = 0.08 * max(y_max - y_min, y_max, 1e-9)
            inset_ax.set_ylim(y_min - padding, y_max + padding)
        if show_inset_title:
            inset_ax.set_title(
                rf"Tail: {inset_xlim[0]:g}–{inset_xlim[1]:g} $\mu$m",
                fontsize=inset_title_size,
                pad=4,
                color="black",
            )
        inset_ax.set_yticks([])
        inset_ax.set_xticks(np.linspace(inset_xlim[0], inset_xlim[1], 3))
        inset_ax.tick_params(
            axis="x",
            labelsize=inset_tick_size,
            width=max(1.0, spine_linewidth - 0.5),
            length=3.5,
            direction="out",
            colors="black",
        )
        for spine in inset_ax.spines.values():
            spine.set_visible(True)
            spine.set_linewidth(spine_linewidth)
            spine.set_color("black")

    fig.tight_layout(pad=1.15)
    _save_figure(
        fig,
        output_dir,
        plot_stem,
        save_formats,
        save_dpi,
    )
    plt.show()
    plt.close(fig)
    return pd.DataFrame(observation_rows), pd.DataFrame(curve_rows)


def summarize_spots_above_threshold(
    tracking_records,
    condition_keys,
    condition_labels,
    threshold_um,
    channel_index=0,
    min_snr=0.0,
    size_filter_mode="trajectory_mean",
):
    """Summarize per-cell trajectories or detections above a threshold."""
    size_filter_mode = _validate_size_filter_mode(size_filter_mode)
    size_column = f"spot_size_ch_{channel_index}"
    snr_column = f"snr_ch_{channel_index}"
    summary_unit = (
        "trajectory" if size_filter_mode == "trajectory_mean" else "detection"
    )
    rows = []

    for key in condition_keys:
        for record in tracking_records[key]:
            try:
                prepared = _prepare_size_filter_table(
                    record["dataframe"],
                    size_column=size_column,
                    pixel_size_um=record["pixel_size_um"],
                    threshold_um=threshold_um,
                    size_filter_mode=size_filter_mode,
                    snr_column=snr_column,
                    min_snr=min_snr,
                )
            except KeyError as exc:
                raise KeyError(
                    f"{exc.args[0]} in {record['tracking_file']}"
                ) from exc

            if size_filter_mode == "trajectory_mean":
                particle_column = _particle_column(prepared)
                units = prepared.drop_duplicates(particle_column)
                above_threshold = units[
                    "_trajectory_above_size_threshold"
                ]
            else:
                units = prepared
                above_threshold = units[
                    "_detection_above_size_threshold"
                ]
            count_above = int(above_threshold.sum())
            total_units = int(len(units))
            rows.append({
                "condition_key": key,
                "condition": condition_labels[key],
                "cell_file": record["cell_name"],
                "tracking_csv": str(record["tracking_file"]),
                "pixel_size_um": record["pixel_size_um"],
                "threshold_um": threshold_um,
                "size_filter_mode": size_filter_mode,
                "summary_unit": summary_unit,
                "units_above_threshold": count_above,
                "total_valid_units": total_units,
                "percentage_above_threshold": (
                    100 * count_above / total_units if total_units else np.nan
                ),
            })
    return pd.DataFrame(rows)


def extract_spots_per_cell(tracking_records, condition_keys):
    """Return one particle-count array per cell for every condition."""
    conditions = []
    for key in condition_keys:
        repetitions = []
        for record in tracking_records[key]:
            dataframe = record["dataframe"]
            particle_column = _particle_column(dataframe)
            count = (
                dataframe[particle_column].nunique()
                if particle_column in dataframe
                else len(dataframe)
            )
            repetitions.append(np.array([count], dtype=float))
        conditions.append(repetitions)
    return conditions


# Figure S10 orchestration and source-data export

def _holm_adjust(p_values):
    """Return Holm-adjusted p-values while preserving missing values."""
    p_values = np.asarray(p_values, dtype=float)
    adjusted = np.full(p_values.shape, np.nan, dtype=float)
    finite_indices = np.flatnonzero(np.isfinite(p_values))
    if not finite_indices.size:
        return adjusted

    finite_p_values = p_values[finite_indices]
    order = np.argsort(finite_p_values)
    sorted_p_values = finite_p_values[order]
    multipliers = len(sorted_p_values) - np.arange(len(sorted_p_values))
    adjusted_sorted = np.maximum.accumulate(multipliers * sorted_p_values)
    adjusted_sorted = np.minimum(adjusted_sorted, 1.0)
    adjusted[finite_indices[order]] = adjusted_sorted
    return adjusted


def _display_to_export_map(condition_keys, display_labels, export_labels):
    """Map plot labels to CSV-safe labels."""
    return {
        display_labels[key]: export_labels[key]
        for key in condition_keys
    }


def _export_cell_source(
    plot_df,
    value_column,
    output_path,
    display_to_export,
):
    """Export one row per cell from a box-and-whisker plot."""
    source_df = (
        plot_df.assign(
            condition=plot_df["Condition"].map(display_to_export)
        )
        .rename(columns={"Mean": value_column})
        [["condition", value_column]]
    )
    if source_df["condition"].isna().any():
        raise ValueError(f"Unmapped condition label while exporting {output_path}")
    source_df.to_csv(output_path, index=False)
    print(f"Saved: {output_path}")
    return source_df


def _format_pairwise_statistics(
    pvalue_df,
    plot_stem,
    display_to_export,
    show_stats,
):
    """Create an export table containing raw and Holm-adjusted p-values."""
    comparison_parts = pvalue_df["comparison"].str.split(
        " vs ", n=1, expand=True
    )
    raw_p_values = pvalue_df["p_value"].to_numpy(dtype=float)
    adjusted_p_values = _holm_adjust(raw_p_values)
    statistics_df = pd.DataFrame({
        "plot": plot_stem,
        "condition_1": comparison_parts[0].map(display_to_export),
        "condition_2": comparison_parts[1].map(display_to_export),
        "test": "Mann-Whitney U",
        "alternative": "two-sided",
        "p_value_raw": raw_p_values,
        "significance_raw": pvalue_df["significance"].to_numpy(),
        "p_value_holm": adjusted_p_values,
        "significance_holm": [
            _significance_label(value) for value in adjusted_p_values
        ],
    })
    if statistics_df[["condition_1", "condition_2"]].isna().any().any():
        raise ValueError(f"Unmapped condition label for {plot_stem}")
    statistics_df["annotation_displayed"] = (
        statistics_df["significance_raw"].isin(["*", "**", "***", "****"])
        if show_stats
        else False
    )
    statistics_df["annotation_basis"] = "raw p-value"
    return statistics_df


def _plot_whisker(
    condition_arrays,
    figure,
    condition_keys,
    display_labels,
    export_labels,
    figures_dir,
    raw_data_dir,
    style,
    source_column=None,
    max_percentile=99.0,
    significance_top_padding=0.035,
):
    """Create one box-and-whisker plot and its export tables."""
    plot_stem = figure["figure_name"]
    points, summary, pvalues = plot_cell_means_box_swarm(
        condition_arrays,
        condition_labels=[display_labels[key] for key in condition_keys],
        output_dir=figures_dir,
        y_label=figure["y_label"],
        plot_stem=plot_stem,
        y_lim=figure["y_limits"],
        show_stats=style["show_stats"],
        figsize=style["figsize"],
        tick_size=style["tick_size"],
        max_percentile_significance=max_percentile,
        significance_top_padding=significance_top_padding,
        save_formats=style["save_formats"],
        save_dpi=style["save_dpi"],
    )

    display_to_export = _display_to_export_map(
        condition_keys,
        display_labels,
        export_labels,
    )
    source_data = None
    if source_column is not None:
        source_path = Path(raw_data_dir) / figure["data_name"]
        source_data = _export_cell_source(
            points,
            source_column,
            source_path,
            display_to_export,
        )

    statistics = _format_pairwise_statistics(
        pvalues,
        plot_stem,
        display_to_export,
        style["show_stats"],
    )
    return {
        "points": points,
        "summary": summary,
        "pvalues": pvalues,
        "statistics": statistics,
        "source_data": source_data,
    }


def run_whisker_analyses(
    tracking_records,
    condition_keys,
    display_labels,
    export_labels,
    figures,
    figures_dir,
    raw_data_dir,
    channel_index,
    selected_frame,
    threshold_um,
    minimum_snr,
    style,
    size_filter_mode="trajectory_mean",
):
    """Run all count and cell-mean box-and-whisker analyses."""
    size_filter_mode = _validate_size_filter_mode(size_filter_mode)
    size_field = f"spot_size_ch_{channel_index}"
    filtered_size_arrays = extract_particle_means_by_cell(
        tracking_records,
        condition_keys,
        size_field,
        max_particle_mean_size_um=threshold_um,
        selected_frame=selected_frame,
        convert_pixels_to_um=True,
        channel_index=channel_index,
        min_snr=minimum_snr,
        size_filter_mode=size_filter_mode,
    )
    filtered_count_arrays = [
        [
            np.array([len(particle_means)], dtype=float)
            for particle_means in condition_cells
        ]
        for condition_cells in filtered_size_arrays
    ]

    condition_arrays = {
        "spots_per_cell": extract_spots_per_cell(
            tracking_records,
            condition_keys,
        ),
        "spots_per_cell_filtered": filtered_count_arrays,
        "spot_intensity": extract_cell_arrays(
            tracking_records,
            condition_keys,
            f"spot_int_ch_{channel_index}",
            selected_frame=selected_frame,
        ),
        "snr": extract_cell_arrays(
            tracking_records,
            condition_keys,
            f"snr_ch_{channel_index}",
            selected_frame=selected_frame,
        ),
        "spot_size": extract_cell_arrays(
            tracking_records,
            condition_keys,
            size_field,
            selected_frame=selected_frame,
            convert_pixels_to_um=True,
        ),
        "spot_intensity_filtered": extract_particle_means_by_cell(
            tracking_records,
            condition_keys,
            f"spot_int_ch_{channel_index}",
            max_particle_mean_size_um=threshold_um,
            selected_frame=selected_frame,
            channel_index=channel_index,
            min_snr=minimum_snr,
            size_filter_mode=size_filter_mode,
        ),
        "snr_filtered": extract_particle_means_by_cell(
            tracking_records,
            condition_keys,
            f"snr_ch_{channel_index}",
            max_particle_mean_size_um=threshold_um,
            selected_frame=selected_frame,
            channel_index=channel_index,
            min_snr=minimum_snr,
            size_filter_mode=size_filter_mode,
        ),
        "spot_size_filtered": filtered_size_arrays,
    }

    filter_source_token = (
        "trajectory_mean_size_filter"
        if size_filter_mode == "trajectory_mean"
        else "detection_size_filter"
    )
    source_columns = {
        "spots_per_cell": "total_unique_tracked_particles",
        "spots_per_cell_filtered": (
            f"eligible_particles_after_{filter_source_token}"
        ),
        "spot_intensity": "mean_spot_intensity_au",
        "snr": "mean_spot_snr",
        "spot_size": "mean_spot_size_um",
        "spot_intensity_filtered": (
            f"mean_particle_spot_intensity_au_after_{filter_source_token}"
        ),
        "snr_filtered": f"mean_particle_snr_after_{filter_source_token}",
        "spot_size_filtered": (
            f"mean_particle_spot_size_um_after_{filter_source_token}"
        ),
    }
    high_percentile_keys = {
        "spots_per_cell",
        "spots_per_cell_filtered",
        "spot_size",
        "spot_size_filtered",
    }
    results = {}
    statistics = {}
    for key, figure in figures.items():
        if key not in condition_arrays:
            raise KeyError(f"Unknown whisker analysis: {key}")
        plot_figure = dict(figure)
        plot_figure["y_label"] = figure["y_label"].format(
            threshold_um=threshold_um,
            filter_label=_size_filter_label(
                size_filter_mode, threshold_um
            ),
            filter_count_label=_size_filter_count_label(
                size_filter_mode, threshold_um
            ),
        )
        result = _plot_whisker(
            condition_arrays[key],
            plot_figure,
            condition_keys,
            display_labels,
            export_labels,
            figures_dir,
            raw_data_dir,
            style,
            source_column=source_columns[key],
            max_percentile=(99.9 if key in high_percentile_keys else 99.0),
        )
        results[key] = result
        statistics[key] = result["statistics"]
    return results, statistics


def _export_distribution_source(
    observations,
    curves,
    export_labels,
    raw_data_dir,
    observation_filename,
    curve_filename,
    observation_column,
    x_column,
    density_column,
):
    """Export particle observations and KDE coordinates."""
    observation_source = (
        observations.assign(
            condition=observations["condition_key"].map(export_labels)
        )
        .rename(columns={"value": observation_column})
        [["condition", observation_column]]
    )
    curve_source = (
        curves.assign(condition=curves["condition_key"].map(export_labels))
        .rename(columns={
            "x_value": x_column,
            "probability_density": density_column,
            "n_observations": "n_particle_means",
        })
        [["condition", x_column, density_column, "n_particle_means"]]
    )
    if (
        observation_source["condition"].isna().any()
        or curve_source["condition"].isna().any()
    ):
        raise ValueError(
            f"Unmapped condition while exporting {observation_filename}"
        )

    raw_data_dir = Path(raw_data_dir)
    observation_path = raw_data_dir / observation_filename
    curve_path = raw_data_dir / curve_filename
    observation_source.to_csv(observation_path, index=False)
    curve_source.to_csv(curve_path, index=False)
    print(f"Saved: {observation_path}")
    print(f"Saved: {curve_path}")
    return observation_source, curve_source


def run_distribution_analyses(
    tracking_records,
    condition_keys,
    display_labels,
    export_labels,
    colors,
    figures,
    figures_dir,
    raw_data_dir,
    channel_index,
    threshold_um,
    minimum_snr,
    style,
    size_filter_mode="trajectory_mean",
):
    """Run all configured particle-mean distribution analyses."""
    size_filter_mode = _validate_size_filter_mode(size_filter_mode)
    filter_source_token = (
        "trajectory_mean_size_filter"
        if size_filter_mode == "trajectory_mean"
        else "detection_size_filter"
    )
    specifications = {
        "intensity_distribution": {
            "field_prefix": "spot_int_ch_",
            "observation_column": "mean_particle_spot_intensity_au",
            "x_column": "spot_intensity_au",
            "density_column": "probability_density_per_au",
        },
        "size_distribution": {
            "field_prefix": "spot_size_ch_",
            "observation_column": "mean_particle_spot_size_um",
            "x_column": "spot_size_um",
            "density_column": "probability_density_per_um",
            "convert_pixels_to_um": True,
        },
        "intensity_distribution_filtered": {
            "field_prefix": "spot_int_ch_",
            "observation_column": (
                "mean_particle_spot_intensity_au_after_"
                f"{filter_source_token}"
            ),
            "x_column": "spot_intensity_au",
            "density_column": "probability_density_per_au",
            "max_spot_size_um": threshold_um,
        },
        "size_distribution_filtered": {
            "field_prefix": "spot_size_ch_",
            "observation_column": (
                f"mean_particle_spot_size_um_after_{filter_source_token}"
            ),
            "x_column": "spot_size_um",
            "density_column": "probability_density_per_um",
            "convert_pixels_to_um": True,
            "max_spot_size_um": threshold_um,
        },
    }

    results = {}
    for key, figure in figures.items():
        if key not in specifications:
            raise KeyError(f"Unknown distribution analysis: {key}")
        specification = specifications[key]
        distributions = extract_particle_distribution(
            tracking_records,
            condition_keys,
            selected_field_prefix=specification["field_prefix"],
            channel_index=channel_index,
            min_snr=minimum_snr,
            convert_pixels_to_um=specification.get(
                "convert_pixels_to_um", False
            ),
            max_spot_size_um=specification.get("max_spot_size_um"),
            size_filter_mode=size_filter_mode,
        )
        plot_options = dict(style)
        setting_names = {
            "figure_name",
            "data_name",
            "kde_data_name",
            "x_label",
            "x_limits",
        }
        plot_options.update({
            name: value
            for name, value in figure.items()
            if name not in setting_names
        })
        observations, curves = plot_distributions(
            distributions,
            condition_labels=display_labels,
            condition_colors=colors,
            output_dir=figures_dir,
            x_label=figure["x_label"].format(
                threshold_um=threshold_um,
                filter_label=_size_filter_label(
                    size_filter_mode, threshold_um
                ),
            ),
            plot_stem=figure["figure_name"],
            xlim=figure["x_limits"],
            **plot_options,
        )
        observation_source, curve_source = _export_distribution_source(
            observations,
            curves,
            export_labels,
            raw_data_dir,
            figure["data_name"],
            figure["kde_data_name"],
            specification["observation_column"],
            specification["x_column"],
            specification["density_column"],
        )
        results[key] = {
            "observations": observations,
            "curves": curves,
            "observation_source": observation_source,
            "curve_source": curve_source,
        }
    return results


def run_large_spot_analysis(
    tracking_records,
    condition_keys,
    display_labels,
    export_labels,
    figure,
    figures_dir,
    raw_data_dir,
    channel_index,
    threshold_um,
    minimum_snr,
    style,
    size_filter_mode="trajectory_mean",
):
    """Plot and export the per-cell percentage above the size threshold."""
    size_filter_mode = _validate_size_filter_mode(size_filter_mode)
    summary_unit = (
        "trajectory" if size_filter_mode == "trajectory_mean" else "detection"
    )
    unit_plural = f"{summary_unit}s"
    spot_summary = summarize_spots_above_threshold(
        tracking_records,
        condition_keys,
        display_labels,
        threshold_um=threshold_um,
        channel_index=channel_index,
        min_snr=minimum_snr,
        size_filter_mode=size_filter_mode,
    )
    threshold_token = f"{threshold_um:.2f}".replace(".", "p")
    count_column = f"{unit_plural}_above_{threshold_token}_um"
    total_column = f"total_valid_{unit_plural}"
    percentage_column = f"percentage_{unit_plural}_above_{threshold_token}_um"
    source_data = (
        spot_summary.assign(
            condition=spot_summary["condition_key"].map(export_labels)
        )
        .rename(columns={
            "units_above_threshold": count_column,
            "total_valid_units": total_column,
            "percentage_above_threshold": percentage_column,
        })
        [[
            "condition",
            "threshold_um",
            "size_filter_mode",
            "summary_unit",
            count_column,
            total_column,
            percentage_column,
        ]]
    )
    source_path = Path(raw_data_dir) / figure["data_name"]
    source_data.to_csv(source_path, index=False)
    print(f"Saved: {source_path}")

    condition_summary = (
        spot_summary.groupby("condition", sort=False)
        .agg(
            cells=("cell_file", "count"),
            units_above_threshold=("units_above_threshold", "sum"),
            total_valid_units=("total_valid_units", "sum"),
            mean_cell_percentage=("percentage_above_threshold", "mean"),
            median_cell_percentage=("percentage_above_threshold", "median"),
        )
        .reset_index()
    )
    condition_summary["pooled_percentage_above_threshold"] = np.where(
        condition_summary["total_valid_units"] > 0,
        100
        * condition_summary["units_above_threshold"]
        / condition_summary["total_valid_units"],
        np.nan,
    )
    condition_summary.insert(1, "size_filter_mode", size_filter_mode)
    condition_summary.insert(2, "summary_unit", summary_unit)

    percentage_arrays = []
    for key in condition_keys:
        values = spot_summary.loc[
            spot_summary["condition_key"] == key,
            "percentage_above_threshold",
        ].to_numpy(dtype=float)
        percentage_arrays.append([
            np.array([value]) for value in values if np.isfinite(value)
        ])

    plot_figure = dict(figure)
    plot_figure["y_label"] = figure["y_label"].format(
        threshold_um=threshold_um,
        filter_unit_plural=unit_plural.title(),
    )
    result = _plot_whisker(
        percentage_arrays,
        plot_figure,
        condition_keys,
        display_labels,
        export_labels,
        figures_dir,
        raw_data_dir,
        style,
        max_percentile=100,
        significance_top_padding=0.12,
    )
    result["spot_summary"] = spot_summary
    result["condition_summary"] = condition_summary
    result["source_data"] = source_data
    return result


def save_pairwise_statistics(statistics_by_plot, output_path):
    """Combine and save the pairwise statistics from all whisker plots."""
    if not statistics_by_plot:
        raise ValueError("No pairwise statistics were supplied")
    statistics = pd.concat(statistics_by_plot.values(), ignore_index=True)
    output_path = Path(output_path)
    statistics.to_csv(output_path, index=False)
    print(f"Saved: {output_path}")
    return statistics
