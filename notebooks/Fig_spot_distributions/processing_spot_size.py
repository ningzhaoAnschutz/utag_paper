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
    "extract_cell_arrays",
    "extract_particle_distribution",
    "extract_spots_per_cell",
    "load_tracking_records",
    "plot_cell_means_box_swarm",
    "plot_distributions",
    "summarize_spots_above_threshold",
]


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
):
    """Return one array per cell for every condition."""
    condition_arrays = []
    for key in condition_keys:
        cell_arrays = []
        for record in tracking_records[key]:
            dataframe = record["dataframe"]
            if selected_field not in dataframe.columns:
                raise KeyError(f"{selected_field} missing from {record['tracking_file']}")
            if selected_frame is not None and "frame" in dataframe.columns:
                values = dataframe.loc[
                    dataframe["frame"] == selected_frame,
                    selected_field,
                ].to_numpy(dtype=float)
            else:
                values = dataframe[selected_field].to_numpy(dtype=float)
            values = values[np.isfinite(values)]
            if convert_pixels_to_um:
                values = values * record["pixel_size_um"]
            cell_arrays.append(values)
        condition_arrays.append(cell_arrays)
    return condition_arrays


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
            highest_allowed = upper - cap_height - text_offset - 0.035 * axis_span
            number_of_bars = len(significant_pairs)
            level_spacing = 0.055 * axis_span
            if number_of_bars > 1:
                required_top = lowest_line + level_spacing * (number_of_bars - 1)
                if required_top > highest_allowed:
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
            color="#222222",
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
            color="#111111",
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
):
    """Extract existing CSV values in three modes without pixel recalculation."""
    value_column = f"{selected_field_prefix}{channel_index}"
    snr_column = f"snr_ch_{channel_index}"
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
            if snr_column in dataframe.columns:
                dataframe = dataframe.loc[dataframe[snr_column] >= min_snr].copy()
            dataframe = dataframe.loc[
                dataframe[value_column].notna() & (dataframe[value_column] > 0)
            ].copy()
            if dataframe.empty:
                continue

            if convert_pixels_to_um:
                dataframe[value_column] = (
                    dataframe[value_column].astype(float) * record["pixel_size_um"]
                )

            particle_column = _particle_column(dataframe)
            if particle_column in dataframe.columns:
                means = dataframe.groupby(particle_column)[value_column].mean()
                mean_per_particle.extend(means.to_list())
            else:
                mean_per_particle.append(dataframe[value_column].mean())

            all_timepoints.extend(dataframe[value_column].to_list())
            if "frame" in dataframe.columns:
                valid_frames.extend(dataframe["frame"].to_list())
            prepared_records.append(dataframe)

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
    bins=60,
    save_formats=("png", "svg"),
    save_dpi=600,
):
    """Plot KDE curves with an optional tail inset."""
    fig, ax = plt.subplots(figsize=figsize, facecolor="white")
    curves = []

    for result in distribution_results:
        key = result["condition_key"]
        values = _distribution_values(result, mode)
        if values.size == 0:
            continue

        color = condition_colors[key]
        label = condition_labels[key]
        if values.size > 3 and np.nanstd(values) > 0:
            x_values = np.linspace(values.min(), values.max(), 900)
            try:
                y_values = gaussian_kde(values, bw_method="scott")(x_values)
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

    ax.set_xlabel(x_label, fontsize=22, labelpad=6)
    ax.set_ylabel("Probability Density", fontsize=22, labelpad=6)
    ax.tick_params(
        axis="both",
        which="major",
        labelsize=17,
        width=1.5,
        length=6,
        direction="out",
        colors="#222222",
    )
    ax.set_ylim(bottom=0)
    if xlim is not None:
        ax.set_xlim(xlim)
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
        spine.set_color("#222222")

    ax.legend(
        loc="lower center",
        bbox_to_anchor=(0.5, 1.025),
        ncol=len(condition_labels),
        frameon=False,
        fontsize=14.5,
        handlelength=2.5,
        handletextpad=0.7,
        columnspacing=1.45,
        borderaxespad=0,
    )

    if show_inset and curves:
        ax.axvspan(
            inset_xlim[0],
            inset_xlim[1],
            color="#6B7280",
            alpha=0.025,
            linewidth=0,
            zorder=0,
        )
        inset_ax = ax.inset_axes(inset_bounds)
        inset_ax.set_facecolor("#FAFAFA")
        inset_y_values = []
        for x_values, y_values, color in curves:
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
            inset_ax.set_ylim(max(0, y_min - padding), y_max + padding)
        inset_ax.set_title(
            rf"Tail: {inset_xlim[0]:g}–{inset_xlim[1]:g} $\mu$m",
            fontsize=11.5,
            pad=4,
            color="#333333",
        )
        inset_ax.set_yticks([])
        inset_ax.set_xticks(np.linspace(inset_xlim[0], inset_xlim[1], 3))
        inset_ax.tick_params(
            axis="x",
            labelsize=11.5,
            width=1.0,
            length=3.5,
            direction="out",
            colors="#333333",
        )
        for spine in inset_ax.spines.values():
            spine.set_visible(True)
            spine.set_linewidth(0.9)
            spine.set_color("#A7A7A7")

    fig.tight_layout(pad=1.15)
    _save_figure(
        fig,
        output_dir,
        f"{plot_stem}_mode{mode}",
        save_formats,
        save_dpi,
    )
    plt.show()
    plt.close(fig)


def summarize_spots_above_threshold(
    tracking_records,
    condition_keys,
    condition_labels,
    threshold_um,
    channel_index=0,
    min_snr=0.0,
):
    """Summarize the per-cell fraction of particle means above a threshold."""
    size_column = f"spot_size_ch_{channel_index}"
    snr_column = f"snr_ch_{channel_index}"
    rows = []

    for key in condition_keys:
        for record in tracking_records[key]:
            dataframe = record["dataframe"].copy()
            if snr_column in dataframe.columns:
                dataframe = dataframe.loc[dataframe[snr_column] >= min_snr].copy()
            dataframe = dataframe.loc[
                dataframe[size_column].notna() & (dataframe[size_column] > 0)
            ].copy()

            particle_column = _particle_column(dataframe)
            particle_mean_pixels = dataframe.groupby(particle_column)[size_column].mean()
            particle_mean_um = particle_mean_pixels * record["pixel_size_um"]
            count_above = int((particle_mean_um > threshold_um).sum())
            total_particles = int(len(particle_mean_um))
            rows.append({
                "condition_key": key,
                "condition": condition_labels[key],
                "cell_file": record["cell_name"],
                "tracking_csv": str(record["tracking_file"]),
                "pixel_size_um": record["pixel_size_um"],
                "threshold_um": threshold_um,
                "spots_above_threshold": count_above,
                "total_valid_spots": total_particles,
                "percentage_above_threshold": (
                    100 * count_above / total_particles if total_particles else np.nan
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
